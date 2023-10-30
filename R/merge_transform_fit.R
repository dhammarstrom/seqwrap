#' A function to fit models with a chosen fitting algorithm. To be used in seqwrap.
#' @param FUN Name of a fitting function, like glmmTMB::glmmTMB
#' @param ARG A list of arguments that can be evaluated by the fitting function
fit_fun <- function(FUN, ARG) {

  # Fit the model using the selected machinery and available arguments
  model <- do.call(FUN, ARG)

  # return the model
  return(model)

}






#' Transform, merge and fit models. The function is used inside seq_wrapper to combine
#' metadata with target-level data and perform the model fitting. The function is used in a
#' call to pbapply::pblapply.
#' @param x A data frame of target-specific quantities, in seqwrap a list from create_list is used iteratively.
#' @param samp_name Sample names from the upper level function
#' @param metadat Metadata from the upper level function
#' @param add_vars Additional variables to keep from the metadata
#' @param arg_list Arguments from the upper level function
#' @param mt_summary_fun Summary function from the upper level function
#' @param mt_eval_fun Evaluation function from the upper level function
#' @return_mod Logical, should the models be returned as part of the results?
#' @param save_mods Logical, should the models be saved?
#' @param mod_path Path to save the models
#' @param ffun the fitting function from the upper level function
merge_transform_fit <- function(x,
                            samp_name = samplename,
                            metdat = metadata,
                            arg_list = arguments,
                            add_vars = additional_vars,
                            mt_summary_fun = summary_fun,
                            mt_eval_fun = eval_fun,
                            ffun = fitting_fun,
                            return_mod = return_models,
                            save_mods = save_mods,
                            mod_path = mod_path) {


  transposed <- data.frame(seq_sample_id = rownames(t(x[,-1])),
                           y = as.vector(t(x[,-1])))

  colnames(transposed)[1] <- samp_name

  df <- merge(transposed, data.frame(metdat), by = samp_name)

  ## Keep only data needed for fitting
  if("formula" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$formula))
  if("model" %in% names(arg_list))  parsed <- all.vars(as.formula(arg_list$model))
  if(all(c("fixed", "random") %in% names(arg_list)))  parsed <- c(all.vars(as.formula(arg_list$fixed)), all.vars((arg_list$random[[1:length(arg_list$random)]])))

  ## Keep also additional variables that exists in the meta data data set
  if(!is.null(add_vars)) parsed <- c(parsed, add_vars)


  df <- df[,parsed, drop = FALSE]


  ## Remove attributes from the list of arguments (this solves an issue when using glmmTMB)
  arguments_final <- append(arg_list, list(data = df))

  for(i in 1:length(arguments_final)) {

    if(class(arguments_final[[i]]) == "formula") environment(arguments_final[[i]]) <- NULL

  }

  # Add warning/errors to outputs
  warn <- NULL
  err <- NULL
  mod <- NULL
  # Adding null values to outputs
  mod_sum <- NULL
  mod_eval <- NULL

  ## Fit the model
  tryCatch(
    mod <- seqwrap:::fit_fun(ffun, arguments_final),

    warning = function(w) warn <<- w,
    error = function(e) err <<- e
  )

  ## Do summarize function if it exists
  if(!is.null(mt_summary_fun)) {

    tryCatch(
      mod_sum <- do.call("mt_summary_fun", list(x = mod)),

      warning = function(w) warn_mod <<- w,
      error = function(e) err_mod <<- e
    )
  }

  ## Do evaluation function if it exists
  if(!is.null(mt_eval_fun)) {

    tryCatch(
      mod_sum <- do.call("mt_eval_fun", list(x = mod)),

      warning = function(w) warn_eval <<- w,
      error = function(e) err_eval <<- e
    )
  }

  mod_path <- if(!is.null(mod_path)) mod_path else paste0(getwd(), "/seqwrap-output")
  if(save_mods) saveRDS(mod, file = paste0(mod_path, "/", names(x), ".rds"))

  # Return the model if requested
  if(return_mod) {
    return(list(model = mod, summaries = mod_sum, evaluation = mod_eval, warn = warn, err = err))
  } else {
    return(list(summaries = mod_sum, evaluation = mod_eval, warn = warn, err = err))
  }





}


