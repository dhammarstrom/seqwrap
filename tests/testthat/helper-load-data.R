# Packages
library(testthat)
library(seqwrap)


# Use sample data for seqwrap tests
seqdata <- rna_seq_sample |>
  dplyr::mutate(across(-transcript_id, ~ as.integer(round(.x, 0))))

seqdatasubset <- seqdata[1:10, ]

metadata <- rna_seq_metadata |>
  dplyr::filter(!is.na(seq_sample_id))



