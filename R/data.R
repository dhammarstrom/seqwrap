#' A small (sample) RNA-seq data set used for illustrations
#' in the vignette and testing.
#'
#'
#' @format ## `rna_seq_sample`
#' A data frame with 1000 rows and 197 columns
#' \describe{
#'  \item{transcript_id}{Ensembl transcript id}
#'  \item{seq_sample_id}{Each remaining column represents sample id}
#' }
#'
#'

#'
#'
"rna_seq_sample"



#' Metadata to match the small RNA-seq data set used for illustrations
#' in the vignette and testing.
#'
#' @format A data frame with 268 rows and 10 variables:
#' \describe{
#'   \item{participant}{participant id}
#'   \item{sex}{mala/female}
#'   \item{include}{incl if included in analysis}
#'   \item{rnaseq_include}{incl if included in RNA-seq data set}
#'   \item{condition}{Multiple, 3 sets per exercise; Single, 1 set per exercise}
#'   \item{leg}{Right, R; Left, L}
#'   \item{time}{Study time-point: w0, w2pre, w2post, w12}
#'   \item{sample.weight}{Sample wet weight (mg)}
#'   \item{total.rna}{Total amount of RNA (ng)}
#'   \item{rqi}{RNA quality indicator}
#'   \item{seq_sample_id}{Sample ID in datasets downloaded with
#'   the function download_ome()}}
"rna_seq_metadata"
