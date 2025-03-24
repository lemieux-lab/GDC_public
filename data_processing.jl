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



