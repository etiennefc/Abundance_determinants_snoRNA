library(recount3)
library(recount)

# Find all available mouse projects
mouse_projects <- available_projects("mouse")

# Search for specific study of Shen et al. 2012 (mouse tissue RNA-Seq, project_id=SRP006787)
proj_info <- subset(mouse_projects, project == "SRP006787" &
                    project_type == "data_sources")

# Create a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene <- create_rse(proj_info, "gene")


# First, scale counts by total read coverage per sample (saved as counts matrix (must stay ass "counts", otherwise it doesn't work))
assay(rse_gene, "counts") <- transform_counts(rse_gene)

# Then compute TPM from the scaled counts
assay(rse_gene, "TPM") <- getTPM(rse_gene, length_var = "bp_length")

write.csv(assay(rse_gene, "TPM"), snakemake@output[["dataset"]])
#write.csv(assay(rse_gene, "TPM"), snakemake@output[["dataset"]], row.names = FALSE)
