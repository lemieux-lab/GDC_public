using CSV
using DataFrames
using Random

function create_tgca_fetch_data_file(sample_file, outpath)
    baseurl = "https://api.gdc.cancer.gov/data"
    outputfile = "tcga_fetch_data.sh"
    
    SAMPLES = CSV.read("gdc_sample_sheet.2025-03-21.tsv", DataFrame, delim="\t")
    
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

create_tgca_fetch_data_file("gdc_sample_sheet.2025-03-21.tsv", "GDC_raw")





using HDF5
function load_tcga_data(fname)

    fid = h5open(fname, "r")

    data = read(fid, "data_matrix")
    samples = read(fid, "case_id")
    genes = read(fid, "gene_id")
    gene_type = read(fid, "gene_type")
    cancer_types = read(fid, "cancer_type")

    close(fid)

    return data, samples, genes, gene_type, cancer_types
end

data,samples, genes, gene_type, cancer_types = load_tcga_data("full_TCGA_GE.h5")


data


cancer_idx = cancer_types .== "TCGA-BRCA"

function cancer_subset(data, samples, genes, gene_type, cancer_types, target_type, outfile)
    cancer_idx = cancer_types .== target_type

    fid = h5open(outfile, "w")

    fid["data_matrix"] = data[cancer_idx, : ]
    fid["samples"] = samples[cancer_idx]
    fid["gene_id"] = genes
    fid["gene_type"] = gene_type
    fid["cancer_type"] = cancer_types[cancer_idx]

    close(fid)
end

cancer_subset(data,samples, genes, gene_type, cancer_types, "TCGA-BRCA", "BRCA_TCGA_GE.h5")

#TODO: add this to the cancer_subset h5 for PAM50
pam50 = CSV.read("mmc2.csv", DataFrame)
pam50_labels = pam50[!,3]
count_map = Dict([(label, sum(pam50_labels .== label)) for label in unique(pam50_labels)])


length(pam50_labels)