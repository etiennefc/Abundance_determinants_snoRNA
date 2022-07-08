library(rlang)
sessionInfo()
library(branchpointer)
library(data.table)
library(tidyverse)
library(GenomicRanges)

# Get the unique (distinct) transcript_id of all host genes into a character vector 'transcript_ids'
#transcript_df <- read.table(snakemake@input[["transcript_id_df"]], header=TRUE, sep='\t')
transcript_df <- fread(snakemake@input[["transcript_id_df"]], header=TRUE, sep='\t')
transcript_ids <- transcript_df %>% distinct(transcript_id_host, .keep_all=TRUE)

#Remove problematic host genes for branchpointer (they only work when not querried with other host genes in a character vector; no explanation for this...)
#Otherwise, makeBranchpointWindowForExons doesn't work properly for these host genes
prob_hg <- c("ENST00000526015", "ENST00000379060", "ENST00000359074", "ENST00000413987", "ENST00000471759", "ENST00000554429", "ENST00000638012")
transcript_ids_2 <- transcript_ids[!(transcript_ids[["transcript_id_host"]] %in% prob_hg)]

# Create a gtf out of the exons of all host genes
exons <- gtfToExons(snakemake@input[["hg_gtf"]])
write.csv(exons, 'data/test_EXONS.csv', row.names=FALSE)

# Get bedtools bin directory (needed for predictBranchpoints below)
file_bedtools <- file(snakemake@input[["bedtools_dir"]],"r")
bedtools_dir <- readLines(file_bedtools,n=1)
close(file_bedtools)

# Create a window table for each problematic HG (the windows are 18 to 44 intronic nucleotides upstream of the exons in HG)
# and concat vertically (append) all these windows in a single Grange object called window_vec
# The windows are GRanges objects, not dataframes
window_vec <- GRanges()
for (i in prob_hg){
  temp_window <- makeBranchpointWindowForExons(i, idType="transcript_id", exons=exons)
  window_vec <- append(window_vec, temp_window)
}

# Create a window table for all other remaining HG that are not problematic for branchpointer and append it to window_vec
window_table <- makeBranchpointWindowForExons(transcript_ids_2[["transcript_id_host"]], idType="transcript_id", exons=exons)
window_vec <- append(window_vec, window_table)
write.csv(window_vec, snakemake@output[["bp_window_table"]], row.names = FALSE)

# Predict all possible branchpoints and their probability within each window (for all introns of snoRNA host genes)
# Cannot use locally the arguments useParallel or cores which make this script quite long to run (~1h) (did not test on a computer cluster) ...
bp <- predictBranchpoints(window_vec, queryType = "region", genome = snakemake@input[["genome"]], bedtoolsLocation=bedtools_dir)
write.csv(bp, snakemake@output[["bp_distance"]], row.names = FALSE)

##This is to visualize specific HG transcripts within a GRange object such as window_vec
#test <- c("ENST00000526015")
#print(window_vec[elementMetadata(window_vec)[, "transcript_id"] %in% test])
