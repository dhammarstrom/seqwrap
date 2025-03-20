#' Check data frames
#'
#' This function checks data frames. In seqwrap certain patterns
#' are assumed, e.g. first column target identifications and following
#' columns of the same class.
#' @param d a data frame
#' @keywords internal
#'

df <- data.frame(id = c("aa", "bb", "cc"),
                 sample1 = as.integer(c(1, 2, 3)),
                 sample2 = as.integer(c(56, 58, 50)))

df_rn <- data.frame(row.names = c("aa", "bb", "cc"),
                    sample1 = as.integer(c(1, 2, 3)),
                    sample2 = as.integer(c(56, 58, 50)))


check_df <- function(df, expect.rownames = FALSE) {

  # expect.rownames is set from upper level function to
  # make it possible to characterize target id.
  if (expect.rownames) {



  }
  if (!expect.rownames) {
    out <- list(df.class = NULL,
                n.rows = NULL,
                n.cols = NULL,
                target.id = NULL)


  }



  sapply(df, class)



}



#' Data helper function.
#'
#' The function takes a named list of data frames
#' and combines them into target-wise data frames
#' for use in subsequent modelling.
#'
#' If the input is a data frame, the function returns a
#' a list of data frames with a single row (y in subsequent modelling)
#'
#'
#' @param dat Target-wise data, either a data frame or a named list.
#' @keywords internal
data_helper <- function(dat, rownames = FALSE) {
  # Check if input is a list or a data frame
  if (!is.list(dat) && !is.data.frame(dat)) {
    stop("Input must be either a list of data frames or a single data frame")
  }

  # If the input is a data frame
  if (is.data.frame(dat)) {




    if (rownames) {
      cli::cli_inform("Using row names as target identification")

      col_names <- colnames(dat)
      dat[["target_id"]] <- rownames(dat)
      # Move this column to be the first column
      dat <- dat[, c("target_id", col_names)]
      # Reset row names to numeric
      rownames(dat) <- NULL
    }

    # Split the data into a list of data frames for each target
    dfs <- split(dat[, seq_len(ncol(dat))[-1]], dat[, 1])

    dfs <- lapply(dfs, function(x) {
      rownames(x) <- "y"
      return(x)
    })

    return(dfs)
  }

  if (is.list(dat)) {
    if (rownames) {
      cli::cli_inform("Using row names as target identification")
      dat <- lapply(df_list_rn, function(x) {
        col_names <- colnames(x)
        x[["target_id"]] <- rownames(x)
        # Move this column to be the first column
        x <- x[, c("target_id", col_names)]
        # Reset row names to numeric
        rownames(x) <- NULL
        return(x)
      })
    }

    # Check if all elements are data frames
    if (!all(sapply(df_list, is.data.frame))) {
      stop("All elements in the list must be data frames")
    }

    # Get the names of the input data frames
    df_names <- names(df_list)
    if (is.null(df_names)) {
      stop("The list of data frames must be named")
    }

    # Find the number of rows in each data frame
    row_counts <- sapply(df_list, nrow)

    # Check if all data frames have the same number of rows
    if (length(unique(row_counts)) != 1) {
      stop("All data frames must have the same number of rows")
    }

    # Check if all data frames have at least one column
    if (any(sapply(df_list, ncol) < 1)) {
      stop("All data frames must have at least one column")
    }

    # Create a list to store the new data frames
    dfs <- list()

    # For each row index
    for (i in 1:row_counts[1]) {
      # Create a new data frame for this row
      new_df <- data.frame(row.names = df_names)

      # Get the name for this result element (from first column of first data frame)
      first_df <- df_list[[1]]
      first_col_name <- colnames(first_df)[1]
      result_name <- as.character(first_df[i, 1])

      # For each original data frame
      for (j in 1:length(df_list)) {
        df_name <- df_names[j]
        current_df <- df_list[[j]]

        # Extract the row from the current data frame
        row_data <- current_df[i, , drop = FALSE]

        # For each column in the current data frame (skip the first column)
        for (col_name in colnames(current_df)[-1]) {
          # Add the value to the new data frame
          new_df[df_name, col_name] <- row_data[[col_name]]
        }
      }

      # Add the new data frame to the result list with appropriate name
      dfs[[result_name]] <- new_df
    }

    return(dfs)
  }
}
