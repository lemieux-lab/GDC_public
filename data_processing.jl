
include("functions.jl")




filedir="/home/golem/GDC/data/rnaseq/readcount"
colname="stranded_second"
GE_values, file_ids, gene_ensg, gene_symbol, gene_type = merge_tcga_files_fixed(filedir, colname)

file_ids
GE_values



filedir="/home/golem/GDC/data/tables"
outfile="TCGA_counts_fixed.h5"
build_tcga_h5_fixed(filedir, GE_values, file_ids, gene_ensg, gene_symbol, gene_type, outfile)





data, samples, genes, genes_symbol, gene_type, cancer_types, 
    tissue_types = load_tcga_data_fixed("TCGA_counts_fixed.h5")


data
sum(tissue_types .=="Primary Tumor") #... there are more samplesbut fewer tumor?


# TODO: make TCGA


################################
create_tgca_fetch_data_file("gdc_sample_sheet.2025-03-21.tsv", "GDC_raw")


# ~10 mins
@time GE_values, gene_ensg, gene_symbol, gene_type = merge_tcga_files("GDC_raw", "stranded_second")
build_tcga_h5("GDC_raw", "gdc_sample_sheet.2025-03-21.tsv", GE_values, gene_ensg, gene_symbol, gene_type, "tcga_GE_counts.h5")


@time data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types = load_tcga_data("tcga_GE_counts.h5")



data
sum(tissue_types .== "Tumor")


cancer_subset(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types,
                "TCGA-BRCA", "BRCA_TCGA_counts.h5")

cancer_subset(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types,
                "TCGA-BRCA", "BRCA_TCGA_pam50_Thennavan.h5")


data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types = load_tcga_data("BRCA_TCGA_pam50_Thennavan.h5")


pam50 = CSV.read("mmc2.csv", DataFrame)
f = h5open("BRCA_TCGA_pam50_Thennavan.h5", "cw")
pam50_labels = sort(leftjoin(DataFrame(Dict(:CLID => samples)),pam50, on=:CLID), :CLID)[:,3]
pam50_labels[ismissing.(pam50_labels)] .= "NA"
f["PAM50"] = String.( pam50_labels)
close(f)


f = h5open("BRCA_TCGA_pam50_Thennavan.h5", "r")
unique(read(f, "PAM50"))
close(f)





####### TPM #############
@time GE_values, gene_ensg, gene_symbol, gene_type = merge_tcga_files("GDC_raw", "tpm_unstranded")
build_tcga_h5("GDC_raw", "gdc_sample_sheet.2025-03-21.tsv", GE_values, gene_ensg, gene_symbol, gene_type, "tcga_GE_tpm.h5")


data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types = load_tcga_data("tcga_GE_tpm.h5")


cancer_subset(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types,
                "TCGA-BRCA", "BRCA_TCGA_tpm.h5")

cancer_subset(data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types,
                "TCGA-BRCA", "BRCA_TCGA_tpm_pam50_Thennavan.h5")


data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types = load_tcga_data("BRCA_TCGA_tpm_pam50_Thennavan.h5")




pam50 = CSV.read("mmc2.csv", DataFrame)

f = h5open("BRCA_TCGA_tpm_pam50_Thennavan.h5", "cw")
pam50_labels = sort(leftjoin(DataFrame(Dict(:CLID => samples)),pam50, on=:CLID), :CLID)[:,3]
pam50_labels[ismissing.(pam50_labels)] .= "NA"
f["PAM50"] = String.( pam50_labels)
close(f)


f = h5open("BRCA_TCGA_tpm_pam50_Thennavan_surv.h5", "r")
read(f, "Sample ID")
close(f)




data, samples, file_uuid, genes, genes_symbol, gene_type, cancer_types, tissue_types = load_tcga_data("tcga_GE_counts.h5")