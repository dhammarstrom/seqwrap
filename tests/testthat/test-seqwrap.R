

## Testing seqwrap

# Use sample data for tests
seqdata <- rna_seq_sample %>%
  mutate(across(-transcript_id, ~as.integer(round(.x, 0))))

seqdatasubset <- seqdata[1:10,]

metadata <- rna_seq_metadata |>
  dplyr::filter(!is.na(seq_sample_id))


test_that("seqwrap returns a list of models in the model slot when asked to return models", {

  test.glmmTMB <- seqwrap::seqwrap(fitting_fun = glmmTMB::glmmTMB,
                           arguments = list(formula = y ~ time + (1|participant),
                                            family = glmmTMB::nbinom2),
                           data = seqdatasubset,
                           metadata = metadata,
                           samplename = "seq_sample_id",
                           additional_vars = NULL,
                           summary_fun = NULL,
                           eval_fun = NULL,
                           exported = list(),
                           return_models = TRUE,
                           save_models = FALSE,
                           model_path = NULL,
                           subset = NULL,
                           cores = 1)

  expect_s3_class(test.glmmTMB$models[[1]], "glmmTMB")

  test.glm.nb <- seqwrap::seqwrap(fitting_fun = MASS::glm.nb,
                                   arguments = list(formula = y ~ time),
                                   data = seqdatasubset,
                                   metadata = metadata,
                                   samplename = "seq_sample_id",
                                   additional_vars = NULL,
                                   summary_fun = NULL,
                                   eval_fun = NULL,
                                   exported = list(),
                                   return_models = TRUE,
                                   save_models = FALSE,
                                   model_path = NULL,
                                   subset = NULL,
                                   cores = 1)


  expect_s3_class(test.glm.nb$models[[1]], c("glm", "lm", "negbin"))

  test.lm <- seqwrap::seqwrap(fitting_fun = stats::lm,
                              arguments = list(formula = y ~ time),
                              data = seqdatasubset,
                              metadata = metadata,
                              samplename = "seq_sample_id",
                              additional_vars = NULL,
                              summary_fun = NULL,
                              eval_fun = NULL,
                              exported = list(),
                              return_models = TRUE,
                              save_models = FALSE,
                              model_path = NULL,
                              subset = NULL,
                              cores = 1)

  expect_s3_class(test.lm$models[[1]], "lm")

  test.glm <- seqwrap::seqwrap(fitting_fun = stats::glm,
                              arguments = list(formula = y ~ time,
                                               family = poisson(link = "log")),
                              data = seqdatasubset,
                              metadata = metadata,
                              samplename = "seq_sample_id",
                              additional_vars = NULL,
                              summary_fun = NULL,
                              eval_fun = NULL,
                              exported = list(),
                              return_models = TRUE,
                              save_models = FALSE,
                              model_path = NULL,
                              subset = NULL,
                              cores = 1)

  expect_s3_class(test.glm$models[[1]], "glm")

  })


### Model evaluations
test_that("Model summaries and evaluations returns expected results", {

  ## Create a model summary function
  summary_fun_glmmTMB <- function(x) {


    # Extract conditional effects and store as a tibble
    results <-   tibble::as_tibble(coef(summary(x))$cond,
                                   rownames = "coef") |>
      dplyr::select(coef,
                    estimate = Estimate,
                    se = 'Std. Error',
                    z = 'z value',
                    p = 'Pr(>|z|)')


    ## Return results
    return(results)
    }


  eval_fun_glmmTMB <- function(x) {


    simresid <- DHARMa::simulateResiduals(fittedModel = x, plot = F, n = 1000)

    unif <- DHARMa::testUniformity(simresid, plot = F)
    disp <- DHARMa::testDispersion(simresid, plot = F)
    pdHess <- x$sdr$pdHess


    return(data.frame(unif = unif$p.value, disp = disp$p.value, pdHess = pdHess))


    }



  test.summary.glmmTMB <- seqwrap::seqwrap(fitting_fun = glmmTMB::glmmTMB,
                                   arguments = list(formula = y ~ time + (1|participant),
                                                    family = glmmTMB::nbinom2),
                                   data = seqdatasubset,
                                   metadata = metadata,
                                   samplename = "seq_sample_id",
                                   additional_vars = NULL,
                                   summary_fun = summary_fun_glmmTMB,
                                   eval_fun = eval_fun_glmmTMB,
                                   exported = list(),
                                   return_models = TRUE,
                                   save_models = FALSE,
                                   model_path = NULL,
                                   subset = NULL,
                                   cores = 1)

  ## Expected output from summary/evaluation functions
  expect_s3_class(test.summary.glmmTMB$summaries[[1]], "tbl")
  expect_s3_class(test.summary.glmmTMB$evaluations[[1]], "data.frame")

  ## Expect no errors in the error data
  expect_null(test.summary.glmmTMB$errors$err_sum[[1]])
  expect_null(test.summary.glmmTMB$errors$err_eval[[1]])

  ## Expect errors when a bad function is passed
  bad_eval_fun <- function(x) {
    as.data.frame(x)
  }

  bad_summary_fun <- function(x) {
    as.data.frame(x)
  }

  test.summary_bad.glmmTMB <- seqwrap::seqwrap(fitting_fun = glmmTMB::glmmTMB,
                                           arguments = list(formula = y ~ time + (1|participant),
                                                            family = glmmTMB::nbinom2),
                                           data = seqdatasubset,
                                           metadata = metadata,
                                           samplename = "seq_sample_id",
                                           additional_vars = NULL,
                                           summary_fun = bad_summary_fun,
                                           eval_fun = bad_eval_fun,
                                           exported = list(),
                                           return_models = TRUE,
                                           save_models = FALSE,
                                           model_path = NULL,
                                           subset = NULL,
                                           cores = 1)

  test.summary_bad.glmmTMB$errors$err_sum[[1]] |> expect_not_null()


})








