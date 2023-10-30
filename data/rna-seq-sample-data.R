
# A small RNA-seq data set used for illustrations in the vignette and testing.
#

## Using trainomeMetaData
library(trainomeMetaData); library(tidyverse)

## A complete data se
# compl_data <- read_csv("ome-data/vol_transcript_rsem.csv")

# # Keep expressed transcripts
# compl_data <- compl_data[rowSums(compl_data[,-1]) > 100,]
#
# # Sample 100 rows
# rna_seq_sample <- compl_data[sample(1:nrow(compl_data), 1000), ]
#
#
# usethis::use_data(rna_seq_sample)
#
# ## Metadata
# rna_seq_metadata <- vol_samples
#
# usethis::use_data(rna_seq_metadata)
