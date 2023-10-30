#' Transform, merge and fit models. The function is used inside seq_wrapper to combine
#' metadata with target-level data and perform the model fitting. The function is used in a
#' call to pbapply::pblapply.
#' @param x A data frame of target-specific quantities
#' @param samp_name Sample names from the upper level function
#' @param arg_list Arguments from the upper level function
#' @param ffun the fitting function from the upper level function
transform_merge_fit <- function(x,
                                samp_name = samplename,
                                metdat = metadata,
                                arg_list = arguments,
                                add_vars = additional_vars,
                                tmf_summary_fun = summary_fun,
                                tmf_eval_fun = eval_fun,
                                ffun = fitting_fun,
                                save_mods = save_mods,
                                mod_path = mod_path) {


  transposed <- data.frame(seq_sample_id = rownames(t(x[,-1])), y = as.vector(t(x[,-1])))

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

  tryCatch(
    mod <- seqwrap:::fit_fun(ffun, arguments_final),
    warning = function(w) warn <<- w,
    error = function(e) err <<- e
  )


  if(!is.null(tmf_summary_fun)) mod_sum <- do.call(tmf_summary_fun, list(mod))



  if(!is.null(tmf_eval_fun)) mod_eval <- do.call(tmf_eval_fun, list(mod))





  mod_path <- if(!is.null(mod_path)) mod_path else paste0(getwd(), "/seqwrap-output")
  if(save_mods) saveRDS(mod, file = paste0(mod_path, "/", names(x), ".rds"))

  return(list(summaries = mod_sum, evaluation = mod_eval, warn = warn, err = err))




}
