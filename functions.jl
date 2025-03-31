using Pkg
Pkg.activate(".")

using CSV
using DataFrames
using Random
using HDF5
#! ONE file has a different annotation
using ProgressBars


#! FIXME: there is no check on the checksum. This may lead to corrupted files
function create_tgca_fetch_data_file(sample_file,outpath)
    baseurl = "https://api.gdc.cancer.gov/data"
    outputfile = "tcga_fetch_data.sh"
    
    SAMPLES = CSV.read(sample_file, DataFrame, delim="\t")
    
    names(SAMPLES)
    files_uuid = SAMPLES[!,"File ID"]
    sample_id = SAMPLES[!,"Sample ID"]
    file_name = SAMPLES[!,"File Name"]

    f = open(outputfile, "w")
    
    for (idx,file) in enumerate(files_uuid)
        random_num = rand(Float64)
        cmd = "if ! [ -f \"$outpath/$(sample_id[idx]).tsv\" ]; then curl --remote-header-name $baseurl/$file -o $outpath/$(sample_id[idx]).tsv; sleep $random_num;fi\n"
        write(f, cmd)
    end

    close(f)
    return outputfile
end 

############################################### 

function merge_tcga_files_fixed(filedir, colname)
    dirnames = readdir(filedir)
    filter!(x -> occursin("TCGA",x), dirnames)

    n_files=0
    for d in dirnames
        n_files+=length(readdir("$filedir/$d"))
    end
    
    fullpath_ex = "$filedir/$(dirnames[1])/$(readdir("$filedir/$(dirnames[1])")[1])/latest/"
    example_file = "$fullpath_ex/$(readdir("$fullpath_ex")[1])"
    example_csv= CSV.read(example_file, DataFrame, header=2, skipto=7)
    n_genes = size(example_csv)[1]

    gene_ensg = example_csv[!,"gene_id"]
    gene_symbol = example_csv[!,"gene_name"]
    gene_type = example_csv[!,"gene_type"]

    GE_values = Matrix{Float32}(undef, n_files,n_genes)
    # Maybe symbol
    file_ids = Array{String}(undef, n_files)

    #make some info 
    curr_idx = 1
    for dir in dirnames
        type_dirs = readdir("$filedir/$dir")
        for f in ProgressBar(type_dirs)

            # Because some directories are not organised the same way
            halfpath = "$filedir/$dir/$f/"
            pathlists = readdir("$halfpath")

            if length(pathlists) > 1
                fullpath = "$filedir/$dir/$f/latest/"
            else
                fullpath = "$filedir/$dir/$f/$(pathlists[1])/"
            end

            tsv_file = "$fullpath/$(readdir("$fullpath")[1])"

            csvfile = CSV.read(tsv_file, DataFrame, header=2, skipto=7, ntasks=1)
            vals = csvfile[!,"$colname"]
            GE_values[curr_idx, :] = Float32.(vals)

            curr_idx+=1
        end
    end
    
    # This cna be put in the other loop when slowdown is figured out. 
    curr_idx = 1
    for (i,dir) in ProgressBar(enumerate(dirnames))
        type_dirs = readdir("$filedir/$dir")
        for (j, f) in ProgressBar(enumerate(type_dirs))
            file_ids[curr_idx] = f 
            curr_idx+=1
        end
    end

    #file_ids
    return GE_values, file_ids, gene_ensg, gene_symbol, gene_type
end


function build_tcga_h5_fixed(filedir, GE_values, file_ids, gene_ensg, gene_symbol, gene_type, outfile)

    #! Do the recursive search again
    
    dirnames = readdir(filedir)
    filter!(x -> occursin("TCGA",x), dirnames)
    

    cancer_types=[]
    tissue_types=[]

    tissue_types
    for dir in ProgressBar(dirnames)
        table_path = "$filedir/$dir/readcount.tsv"
        SAMPLES = CSV.read(table_path, DataFrame, delim="\t")

        sort!(SAMPLES, ["Project_ID", "Aliquot_ID"])

        cancer_types = vcat(cancer_types, SAMPLES[!,"Project_ID"])
        tissue_types = vcat(tissue_types, SAMPLES[!,"SampleType"])

    end

    fid = h5open(outfile, "w")

    fid["Sample ID"] = Array(file_ids)
    fid["Cancer Type"] = Array(String.(cancer_types)) #!
    fid["Tissue Type"] = Array(String.(tissue_types)) #!
    fid["Data Matrix"] = GE_values
    fid["Gene ENSG"] = Array(gene_ensg)
    fid["Gene Symbol"] = Array(gene_symbol)
    fid["Gene Type"] = Array(gene_type)

    close(fid)

end


function load_tcga_data_fixed(fname)

    fid = h5open(fname, "r")

    data = read(fid, "Data Matrix")
    samples = read(fid, "Sample ID")

    genes = read(fid, "Gene ENSG")
    genes_symbol = read(fid, "Gene Symbol")
    gene_type = read(fid, "Gene Type")
    
    cancer_types = read(fid, "Cancer Type")
    tissue_types = read(fid, "Tissue Type")

    close(fid)

    return data, samples, genes, genes_symbol, gene_type, cancer_types, tissue_types
end

function cancer_subset_fixed(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types, 
    target_type, outfile)
    cancer_idx = cancer_types .== target_type

    fid = h5open(outfile, "w")
    fid["Data Matrix"] = data[cancer_idx, : ]
    fid["Sample ID"] = samples[cancer_idx]
    fid["Cancer Type"] = Array(cancer_types)[cancer_idx]
    fid["Tissue Type"] = Array(tissue_types)[cancer_idx]
    fid["Gene ENSG"] = Array(genes)
    fid["Gene Symbol"] = genes_symbol
    fid["Gene Type"] = Array(gene_type)    
    close(fid)
end

###################################


#! TODO: finish this. But I can save the previous stuff into another file. 
function get_from_manifest(filedir)
    dirnames = readdir(filedir)
    filter!(x -> occursin("TCGA",x), dirnames)

    for dir in dirnames
        SAMPLES = CSV.read("$dir/readcount.tsv", DataFrame, delim="\t")
        sort!(SAMPLES, "")
        type_dirs = readdir("$filedir/$dir")
        for f in ProgressBar(type_dirs)
            SAMPLES = CSV.read("$filedir/$dir", DataFrame, delim="\t")
        end
    end
        
    SAMPLES = CSV.read(manifest, DataFrame, delim="\t")
end
filedir="/home/golem/GDC/data/tables"


# TODO: get the manifest file for each project... 
function build_tcga_h5(filedir, manifest, GE_values, gene_ensg, gene_symbol, gene_type, outfile)

    file_order = readdir(filedir)
    SAMPLES = CSV.read(manifest, DataFrame, delim="\t")
        
    sample_id = SAMPLES[!,"Sample ID"]

    order = [findfirst(==(id), sample_id) for id in [i[1] for i in split.(file_order, ".")]]

    sample_id = sample_id[order]

    files_uuid = SAMPLES[!,"File ID"][order]
    cancer_types = SAMPLES[!,"Project ID"][order]
    tissue_types = SAMPLES[!,"Tissue Type"][order]

    cancer_types

    
    fid = h5open(outfile, "w")

    fid["Sample ID"] = Array(sample_id)
    fid["File UUID"] = Array(files_uuid)
    fid["Cancer Type"] = Array(cancer_types)
    fid["Tissue Type"] = Array(tissue_types)
    fid["Data Matrix"] = GE_values
    fid["Gene ENSG"] = Array(gene_ensg)
    fid["Gene Symbol"] = Array(gene_symbol)
    fid["Gene Type"] = Array(gene_type)

    close(fid)

end


function merge_tcga_files(filedir, colname)

    fnames = readdir(filedir)
    
    n_samps = length(fnames)

    example_file = CSV.read("$filedir/$(fnames[1])", DataFrame, header=2, skipto=7)
    n_genes = size(example_file)[1]

    gene_ensg = example_file[!,"gene_id"]
    gene_symbol = example_file[!,"gene_name"]
    gene_type = example_file[!,"gene_type"]

    GE_values = Matrix{Float32}(undef, n_samps,n_genes)

    #make some info 
    for (i,f) in ProgressBar(enumerate(fnames))
        csvfile = CSV.read("$filedir/$f", DataFrame, header=2, skipto=7)
        vals = csvfile[!,"$colname"]
        GE_values[i, :] = Float32.(vals)
    end

    return GE_values, gene_ensg, gene_symbol, gene_type
end



# TODO: make this cleaner. 
function build_tcga_h5(filedir, manifest, GE_values, gene_ensg, gene_symbol, gene_type, outfile)

    file_order = readdir(filedir)
    SAMPLES = CSV.read(manifest, DataFrame, delim="\t")
        
    sample_id = SAMPLES[!,"Sample ID"]

    order = [findfirst(==(id), sample_id) for id in [i[1] for i in split.(file_order, ".")]]

    sample_id = sample_id[order]

    files_uuid = SAMPLES[!,"File ID"][order]
    cancer_types = SAMPLES[!,"Project ID"][order]
    tissue_types = SAMPLES[!,"Tissue Type"][order]

    cancer_types

    
    fid = h5open(outfile, "w")

    fid["Sample ID"] = Array(sample_id)
    fid["File UUID"] = Array(files_uuid)
    fid["Cancer Type"] = Array(cancer_types)
    fid["Tissue Type"] = Array(tissue_types)
    fid["Data Matrix"] = GE_values
    fid["Gene ENSG"] = Array(gene_ensg)
    fid["Gene Symbol"] = Array(gene_symbol)
    fid["Gene Type"] = Array(gene_type)

    close(fid)

end


function load_tcga_data(fname)

    fid = h5open(fname, "r")

    data = read(fid, "Data Matrix")
    samples = read(fid, "Sample ID")
    file_uuid = read(fid, "File UUID")

    genes = read(fid, "Gene ENSG")
    genes_symbol = read(fid, "Gene Symbol")
    gene_type = read(fid, "Gene Type")
    
    cancer_types = read(fid, "Cancer Type")
    tissue_types = read(fid, "Tissue Type")

    close(fid)

    return data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types
end


function cancer_subset(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types, 
    target_type, outfile)
    cancer_idx = cancer_types .== target_type

    fid = h5open(outfile, "w")
    fid["Data Matrix"] = data[cancer_idx, : ]
    fid["Sample ID"] = samples[cancer_idx]
    fid["File UUID"] = file_uuid[cancer_idx]
    fid["Cancer Type"] = Array(cancer_types)[cancer_idx]
    fid["Tissue Type"] = Array(tissue_types)[cancer_idx]
    fid["Gene ENSG"] = Array(genes)
    fid["Gene Symbol"] = genes_symbol
    fid["Gene Type"] = Array(gene_type)    
    close(fid)
end