#' A flexible upper-level wrapper for iterative modelling using any available
#' fitting algorithm.
#'
#' @param fitting_fun A model fitting function like stats::lm,
#' glmmTMB::glmmTMB, lme4::lmer
#' @param arguments A list of arguments to be passed to the fitting function,
#' this should not contain data. Note that the formula must have y as the
#' dependent variable.
#' @param data A data frame with targets (i.e. genes, transcripts) as rows and
#' sample names as colums.
#' The first column is assumed to contain target names/identification
#' @param metadata A data frame with sample names (corresponding to column
#' names in the target matrix)
#' and design variables.
#' @param samplename A character value indicating the variable by which
#' metadata can merge with the target data. This defaults to "seq_sample_id"
#' as this is used in the trainomeMetaData package.
#' @param additional_vars A vector of additional variables that is contained
#' in the metadata data set that is needed to fit the model. By default the
#' metadata is reduced to variables contained in the slots
#' formula/model/fixed/random in additional arguments.
#' More variables may be needed for offsets, weights etc.
#' @param summary_fun A custom (user-created) function for
#' evaluating/summarizing models. If NULL, a list of fitted models are returned
#' @param eval_fun A custom (user-created) function for model
#' diagnostics/evaluation. If NULL, no evaluation/diagnostics of models are
#' returned
#' @param exported A list of functions, values etc. to be passed to
#' summary_fun and eval_fun. This list must contain any functions that
#' should be used in model summarise or evaluations.
#' @param return_models Logical, should models be returned as part of the
#' output? Save models during development on subsets of the data.
#' If used on large data sets, this will result in large memory usage.
#' @param save_models Logical, should models be saved? Models may be saved
#' on disk to save working memory.
#' @param model_path A character. The path to saved models.
#' @param subset A sequence, random samples or integers to indicate which
#' rows to keep in data. This is useful if you want to test the model in a
#' subset of targets. If keft to the default (NULL), all rows will be used.
#' @param cores An integer indicating the number of cores to be used in parallel
#'  computations. If NULL, a sequential for loop is used. If "max", all
#'  available cores are used.
#' @return A nested list with three upper levels slots: models, a list of
#' fitted objects; summaries, a list of summaries created from the summary_fun
#' function; evaluations, a list of diagnostics created from eval_fun.
#' @details This function provides a flexible wrapper to fit, summarize and
#' evaluate statistical models fitted to high dimensional omics-type data.
#' Models are fitted and passed to user defined functions to summarise and
#' evaluate models.
#' @export
seqwrap <- function(fitting_fun,
                    arguments,
                    data,
                    metadata,
                    samplename = "seq_sample_id",
                    additional_vars = NULL,
                    summary_fun = NULL,
                    eval_fun = NULL,
                    exported = list(),
                    return_models = TRUE,
                    save_models = FALSE,
                    model_path = NULL,
                    subset = NULL,
                    cores) {
  ## Prepare data #####
  data <- data.frame(data)
  metadata <- data.frame(metadata)

  ## Subset
  if (!is.null(subset)) data <- data[subset, ]

  ## Sanity checks

  # checks if arguments for the provided fitting function matches arguments
  stopifnot(
    "Arguments do not match named arguments of the selected \nmodel
            fitting function (fitting_fun)." =
      names(arguments) %in% names(formals(glmmTMB::glmmTMB))
  )

  # Check if the data has a character vector for first column
  # (indicating transcript id).
  stopifnot("The first column of the data is not character or
            factor,\ncheck if this column indicate target
            identifications." = is.character(data[, 1]))

  # Check if the sample name is present in the meta data
  stopifnot("The samplename does not exist in the metadata,\nno
            variable for matching metadata and
            data." = samplename %in% colnames(metadata))

  # Check that data is formatted correctly
  if (!all(metadata[, samplename] %in% names(data[, -1]))) {
    stop("The sample names in the metadata does
         not match the sample column names in the seqdata.")
  }


  # Split the data into a list of data frames for each target
  dfs <- split(data[, seq_len(ncol(data))], data[, 1])


  ### Fitting models in parallel #######

  # Determine the number of cores
  if (is.null(cores)) num_cores <- 1
  if (cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if (cores <= parallel::detectCores()) num_cores <- cores


  # Create a cluster using the number of cores specified
  cl <- parallel::makeCluster(num_cores)



  ## Export data to clusters
  parallel::clusterExport(cl,
    varlist = c(
      "metadata",
      "arguments",
      "fitting_fun",
      "samplename",
      "additional_vars",
      "fitting_fun",
      "save_models",
      "exported",
      "model_path",
      "return_models",
      "summary_fun",
      "eval_fun"
    ),
    envir = environment()
  )



  # Parallel execution of the fitting process
  cat("Transforming, merging and modelling data.\n")
  results <- pbapply::pblapply(
    cl = cl,
    X = dfs,
    FUN = seqwrap::merge_transform_fit,
    samp_name = samplename,
    metdat = metadata,
    arg_list = arguments,
    add_vars = additional_vars,
    mt_summary_fun = summary_fun,
    mt_eval_fun = eval_fun,
    ffun = fitting_fun,
    return_mod = return_models,
    save_mods = save_models,
    mod_path = model_path
  )

  parallel::stopCluster(cl)




  # Return results

  models <- NULL
  summaries <- NULL
  evaluations <- NULL
  errors <- NULL

  # Collect models
  if (return_models) models <- lapply(results, `[[`, "model")

  if (!is.null(summary_fun)) summaries <- lapply(results, `[[`, "summaries")
  if (!is.null(eval_fun)) evaluations <- lapply(results, `[[`, "evaluation")

  ## Create a data frame of all errors/warnings
  errors <- tibble::tibble(
    target = names(results),
    errors_fit = lapply(results, `[[`, "err"),
    warnings_fit = lapply(results, `[[`, "warn"),
    err_sum = lapply(results, `[[`, "err_sum"),
    warn_sum = lapply(results, `[[`, "warn_sum"),
    err_eval = lapply(results, `[[`, "err_eval"),
    warn_eval = lapply(results, `[[`, "warn_eval")
  )

  comb_results <- list(
    models = models,
    summaries = summaries,
    evaluations = evaluations,
    errors = errors
  )


  return(comb_results)
}
