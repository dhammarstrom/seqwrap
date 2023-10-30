

## Tests for merge_transform_fit

# Use sample data for tests
seqdata <- rna_seq_sample
metadata <- rna_seq_metadata |>
  dplyr::filter(!is.na(seq_sample_id))

# Use create_lists to create a list of data frames
dfs <- seqwrap::create_lists(seqdata)



## Create a fit manually ############
arguments <- list(formula = y ~ time + (1|participant),
                  family = glmmTMB::nbinom2)

# Transform and merge by hand
dat <- data.frame(seq_sample_id = rownames(t(x[,-1])),
                                y = as.vector(t(x[,-1])))

dat_merged <- merge(dat, data.frame(metadata), by = "seq_sample_id") |>
  mutate(y = as.integer(round(y, 0)))

testmod <- glmmTMB::glmmTMB(data = dat_merged, formula = y ~ time + (1|participant),
                 family = glmmTMB::nbinom2)

summary_fun(testmod)

coef(summary(testmod))$cond


summary_fun <- function(x){
  if(is.null(x$mod)) return(x$err) else {
    cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$cond))),
                               coef(summary(x))$cond,

                               row.names = NULL)

    return(cond_effects)
  }

}

# Extract one row from sample data and convert to expected format
x <- dfs[[1]] %>%
  mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

### Tests ######################################################
test_that("merge_transform_fit returns a list", {


  # Extract one row from sample data and convert to expected format
  x <- dfs[[1]] %>%
    mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

  expect_type(merge_transform_fit(x,
                                  samp_name = "seq_sample_id",
                                  metdat = metadata,
                                  arg_list = arguments,
                                  add_vars = NULL,
                                  mt_summary_fun =  summary_fun,
                                  mt_eval_fun = NULL,
                                  ffun = glmmTMB::glmmTMB,
                                  return_mod = FALSE,
                                  save_mods = FALSE,
                                  mod_path = NULL), "list")
})




test_that("merge_transform_fit returns no error when the fitting fails, but instead returns a message under $err", {


  # Extract one row from sample data and convert to expected format
  x <- dfs[[1]] %>%
    mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

  x_nonumbers <- x
  x_nonumbers[1,2:197] <- "nonumbers"


  expect_no_error(merge_transform_fit(x_nonumbers,
                                   samp_name = "seq_sample_id",
                                   metdat = metadata,
                                   arg_list = arguments,
                                   add_vars = NULL,
                                   mt_summary_fun =  summary_fun,
                                   mt_eval_fun = NULL,
                                   ffun = glmmTMB::glmmTMB,
                                   return_mod = FALSE,
                                   save_mods = FALSE,
                                   mod_path = NULL))


  results <- merge_transform_fit(x_nonumbers,
                      samp_name = "seq_sample_id",
                      metdat = metadata,
                      arg_list = arguments,
                      add_vars = NULL,
                      mt_summary_fun =  summary_fun,
                      mt_eval_fun = NULL,
                      ffun = glmmTMB::glmmTMB,
                      save_mods = FALSE,
                      mod_path = NULL)

  expect_type(results$err, "list")

})


test_that("merge_transform_fit can use multiple fitting algorithms to fit data provided", {







  # Extract one row from sample data and convert to expected format
  x <- dfs[[1]] %>%
    mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

  ## Summary function for glmmTMB ##############################
  summary_fun_glmmtmb <- function(x){
    if(is.null(x$mod)) return(x$err) else {
      cond_effects <- data.frame(cbind(data.frame(coef = rownames(coef(summary(x))$cond))),
                                 coef(summary(x))$cond,

                                 row.names = NULL)

      return(cond_effects)
    }

  }


#### Fitting glmmTMB ##########################################
  expect_no_error(merge_transform_fit(x,
                                      samp_name = "seq_sample_id",
                                      metdat = metadata,
                                      arg_list = list(formula = y ~ time + (1|participant),
                                                           family = glmmTMB::nbinom2),
                                      add_vars = NULL,
                                      mt_summary_fun =  summary_fun_glmmtmb,
                                      mt_eval_fun = NULL,
                                      ffun = glmmTMB::glmmTMB,
                                      save_mods = FALSE,
                                      mod_path = NULL))


  results <- merge_transform_fit(x,
                                 samp_name = "seq_sample_id",
                                 metdat = metadata,
                                 arg_list = list(formula = y ~ time + (1|participant),
                                                 family = glmmTMB::nbinom2),
                                 add_vars = NULL,
                                 mt_summary_fun =  summary_fun_glmmtmb,
                                 mt_eval_fun = NULL,
                                 ffun = glmmTMB::glmmTMB,
                                 return_mod = TRUE,
                                 save_mods = FALSE,
                                 mod_path = NULL)

  expect_s3_class(results$summaries, "data.frame")
  expect_s3_class(results$model, "glmmTMB")


  ### Summary function for glm.nb   ##############################
  summary_fun_glm.nb <- function(x){

      effects <- data.frame(cbind(data.frame(coef = rownames(summary(x)$coef))),
                            summary(x)$coef,

                                 row.names = NULL)

      return(effects)

  }



  # Extract one row from sample data and convert to expected format
  x <- dfs[[1]] %>%
    mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

  ## Fitting MASS::glm.nb ##########################################
  results_glm.nb <- merge_transform_fit(x,
                                 samp_name = "seq_sample_id",
                                 metdat = metadata,
                                 arg_list = list(formula = y ~ time),
                                 add_vars = NULL,
                                 mt_summary_fun =  summary_fun_glm.nb,
                                 mt_eval_fun = NULL,
                                 ffun = MASS::glm.nb,
                                 return_mod = TRUE,
                                 save_mods = FALSE,
                                 mod_path = NULL)



  expect_s3_class(results_glm.nb$summaries, "data.frame")
  expect_s3_class(results_glm.nb$model, "glm")


})

