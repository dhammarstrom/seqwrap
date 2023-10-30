#' @title Create lists of data frames from seqdata
#' @description This function takes a data frame of sequence data list of data frames based on the number of targets.
#' @param seqdata A data frame of sequence data with the first column being target IDs.
#' @return A list of data frames
create_lists <- function(seqdata) {

  ## Prepare data #####
  seqdata <- data.frame(seqdata)

  # Check if the data has a character vector for first column (indicating transcript id).
  if(!is.character(seqdata[,1])) {
    stop("The first column of the data is not character or factor, check if this column indicate target identifications.")
  }

  ## Quantity tables from the trainomeMetaData package comes with the first column
  # being trainscript id's. This splits the data into data frames based on number of
  # targets. Splits are saved as a list of data frames in the output list.

  output <- split(seqdata[,1:ncol(seqdata)], seqdata[,1])

  return(output)
}



