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

function merge_tcga_files(filedir, colname)

    fnames = readdir(filedir)
    
    n_samps = length(fnames)

    example_file = CSV.read("$filedir/$(fnames[1])", DataFrame, header=2, skipto=7)
    n_genes = size(example_file)[1]

    gene_ensg = example_file[!,"gene_id"]
    gene_symbol = example_file[!,"gene_name"]
    gene_type = example_file[!,"gene_type"]

    GE_values = Matrix{Float64}(undef, n_samps,n_genes)

    #make some info 
    for (i,f) in ProgressBar(enumerate(fnames))
        csvfile = CSV.read("$filedir/$f", DataFrame, header=2, skipto=7)
        vals = csvfile[!,"$colname"]
        if length(vals) != n_genes
            println(f)
        end
        GE_values[i, :] = Float64.(vals)
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

    fid["Sample ID"] = sample_id
    fid["File UUID"] = files_uuid
    fid["Cancer Type"] = Array(cancer_types)
    fid["Tissue Type"] = Array(tissue_types)
    fid["Data Matrix"] = GE_values
    fid["Gene ENSG"] = gene_ensg
    fid["Gene Symbol"] = gene_symbol
    fid["Gene Type"] = Array(gene_type)

    close(fid)

end



using HDF5
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