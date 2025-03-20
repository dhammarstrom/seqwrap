# Define a custom property type that accepts either NULL or a list
null_or_list <- S7::new_property(
  # The validator function checks if the value is NULL or a list
  validator = function(value) {
    if (!(is.null(value) || is.list(value))) "must be NULL or a list"
  }
)


#' @export
seqwrapResults <- S7::new_class(
  name = "seqwrap_results",
  parent = S7::class_list,
  properties = list(
    models = null_or_list,
    summaries = null_or_list,
    evaluations = null_or_list,
    errors = S7::class_data.frame,
    n = S7::class_numeric,
    k = S7::class_numeric,
    call_arguments = S7::class_character,
    call_engine = S7::class_character
  )
)


#' Print method for objects of class seqwrapResults
#'
#' @description
#' Invoking the print method on seqwrapResults gives a summary
#' of the fitted objects.
#'
#' @param x A seqwrapResults object
#'
#' @return Invisibly returns the object
#'
#' @examples
#' results <- seqwrap(...)
#' print(results)
#'
#' @method print seqwrapResults
#' @name print.seqwrapResults
S7::method(print, seqwrapResults) <- function(x) {
  cli::cli_h1("seqwrap")
  cli::cli_inform(
    "A total of {x@n} sample{?s} and {x@k} target{?s} where
                  used in {.code {x@call_engine}} with arguments
                  {.code {x@call_arguments}}"
  )

  # Count non-null values in the error/warning data frame
  errors_sum <- sapply(x@errors, function(x) sum(!sapply(x, is.null)))

  if (any(errors_sum[-1] > 0)) {
    cli::cli_alert_info("Some targets had associated errors or warnings")

    cli::cli_inform(c(
      "*" = "Fitting algorithm (errors): n = {errors_sum[2]}
      ({100 * (errors_sum[2]/k)}%)",
      "*" = "Fitting algorithm (warnings): n = {errors_sum[3]}
      ({100 * (errors_sum[3]/k)}%)",
      "*" = "Summary function (errors): n = {errors_sum[4]}
      ({100 * (errors_sum[4]/k)}%)",
      "*" = "Summary function (warnings): n = {errors_sum[5]}
      ({100 * (errors_sum[5]/k)}%)",
      "*" = "Evaluation function (errors): n = {errors_sum[6]}
      ({100 * (errors_sum[6]/k)}%)",
      "*" = "Evaluation function (warnings): n = {errors_sum[7]}
      ({100 * (errors_sum[7]/k)}%)"
    ))
  } else cli::cli_alert_info("No targets had associated errors or warnings")

  invisible(x)
}


#' A flexible upper-level wrapper for iterative modelling using any available
#' fitting algorithm.
#'
#' @param fitting_fun A model fitting function like stats::lm,
#' glmmTMB::glmmTMB, lme4::lmer
#' @param arguments An alist or list of arguments to be passed to the fitting
#'  function, this should not contain data. Note that the formula must have
#'  y as the dependent variable.
#' @param data A data frame or a list of data frames with targets (e.g. genes,
#' transcripts) as rows and sample names as columns.
#' If rownames = FALSE (default), each data frame should have target
#' identifications as the first column in the data frame(s). If rownames = TRUE
#' row names will be converted to target identifications. If data is provided as
#' a list, each element of the list should be named. The corresponding names
#' be available as variables for the fitting function.
#' @param rownames should row names in data be used as target identifications?
#' Defaults to FALSE.
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
seqwrap <- function(
  fitting_fun,
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
  cores
) {
  ## Prepare data #####
  data <- data.frame(data)
  metadata <- data.frame(metadata)

  ## Subset the data
  if (!is.null(subset)) data <- data[subset, ]

  ## Sanity checks

  # checks if arguments for the provided fitting function matches arguments
  if (!all(names(arguments) %in% names(formals(fitting_fun)))) {
    cli::cli_abort(
      "Arguments do not match named arguments of the selected
                   \nmodel fitting function (fitting_fun).",
      call = call
    )
  }

  # Check if the data has a character vector for first column
  # (indicating transcript id).
  if (!is.character(data[, 1])) {
    cli::cli_abort(
      "The first column of the data is not character or
            factor,\ncheck if this column indicate target
            identifications.",
      call = call
    )
  }

  # Check if the sample name is present in the meta data
  if (!samplename %in% colnames(metadata)) {
    cli::cli_abort(
      "The samplename does not exist in the metadata,\nno
            variable for matching metadata and
            data.",
      call = call
    )
  }

  # Check that data is formatted correctly
  if (!all(metadata[, samplename] %in% names(data[, -1]))) {
    cli::cli_abort(
      "The sample names in the metadata does
         not match the sample column names in the seqdata.",
      call = call
    )
  }

  # data_helper function. Combine data into a list of data frames
  # containing variables, y in case of user-provided data frame;
  # names of variables from list names in case of user-provided
  # named list.
  dfs <- data_helper(data)

  ### Combine information for print method ###
  # Store number of targets
  k <- length(dfs)
  # Store the number of samples
  n <- length(metadata[, samplename])

  # Determine the number of cores
  if (is.null(cores)) num_cores <- 1
  if (cores >= parallel::detectCores()) num_cores <- parallel::detectCores()
  if (cores <= parallel::detectCores()) num_cores <- cores

  # Catch the function calls for printing
  funcall <- match.call()
  ## Arguments to string
  call_arguments <- deparse(funcall$arguments)
  call_arguments <- sub("^.*?\\((.*)\\).*$", "\\1", call_arguments)
  ## fitting_fun to string
  call_engine <- deparse(funcall$fitting_fun)

  # Print pre-fit information
  cli::cli_h1("seqwrap")
  cli::cli_inform(
    "Preparing to fit {n} sample{?s} and {k} target{?s} using
                    {call_engine} and {call_arguments}
                    with {num_cores} core{?s}."
  )

  ### Fitting models in parallel #######

  # Create a cluster using the number of cores specified
  cl <- parallel::makeCluster(num_cores)
  ## Export data to clusters
  parallel::clusterExport(
    cl,
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
  cli::cli_h2("Merging and modelling data...")
  results <- pbapply::pblapply(
    cl = cl,
    X = dfs,
    FUN = seqwrap_mtf,
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

  # Combine results

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

  # Count non-null values in the error/warning data frame
  errors_sum <- sapply(errors, function(x) sum(!sapply(x, is.null)))

  ## Evaluate errors for the resulting print function
  cli::cli_h2("Completed")

  if (any(errors_sum[-1] > 0)) {
    cli::cli_alert_info("Some targets had associated errors or warnings")

    cli::cli_inform(c(
      "*" = "Fitting algorithm (errors): n = {errors_sum[2]}
      ({100 * (errors_sum[2]/k)}%)",
      "*" = "Fitting algorithm (warnings): n = {errors_sum[3]}
      ({100 * (errors_sum[3]/k)}%)",
      "*" = "Summary function (errors): n = {errors_sum[4]}
      ({100 * (errors_sum[4]/k)}%)",
      "*" = "Summary function (warnings): n = {errors_sum[5]}
      ({100 * (errors_sum[5]/k)}%)",
      "*" = "Evaluation function (errors): n = {errors_sum[6]}
      ({100 * (errors_sum[6]/k)}%)",
      "*" = "Evaluation function (warnings): n = {errors_sum[7]}
      ({100 * (errors_sum[7]/k)}%)"
    ))
  }

  ## Combine the results into a seqwrapResults
  comb_results <- seqwrapResults(
    models = models,
    summaries = summaries,
    evaluations = evaluations,
    errors = errors,
    n = n,
    k = k,
    call_arguments = call_arguments,
    call_engine = call_engine
  )

  return(comb_results)
}
