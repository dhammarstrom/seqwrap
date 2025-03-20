# Packages
library(testthat)
library(seqwrap)


test_that("seqwrap returns a list of models in the model
          slot when asked to return models", {
  test_glmmtmb <- seqwrap::seqwrap(
    fitting_fun = glmmTMB::glmmTMB,
    arguments = list(
      formula = y ~
        time +
          (1 | participant),
      family = glmmTMB::nbinom2
    ),
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
    cores = 1
  )

  expect_s3_class(test_glmmtmb@models[[1]], "glmmTMB")

  test_glmnb <- seqwrap::seqwrap(
    fitting_fun = MASS::glm.nb,
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
    cores = 1
  )

  expect_s3_class(test_glmnb@models[[1]], c("glm", "lm", "negbin"))

  test_lm <- seqwrap::seqwrap(
    fitting_fun = stats::lm,
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
    cores = 1
  )

  expect_s3_class(test_lm@models[[1]], "lm")

  test_glm <- seqwrap::seqwrap(
    fitting_fun = stats::glm,
    arguments = list(
      formula = y ~ time,
      family = poisson(link = "log")
    ),
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
    cores = 1
  )

  expect_s3_class(test_glm@models[[1]], "glm")
})


test_that("Model summaries and evaluations returns expected results", {
  ## Create a model summary function
  summaryfun_glmmtmb <- function(x) {
    # Extract conditional effects and store as a tibble
    results <- tibble::as_tibble(coef(summary(x))$cond, rownames = "coef") |>
      dplyr::select(
        coef,
        estimate = Estimate,
        se = "Std. Error",
        z = "z value",
        p = "Pr(>|z|)"
      )

    ## Return results
    return(results)
  }

  evalfun_glmmtmb <- function(x) {
    simresid <- DHARMa::simulateResiduals(
      fittedModel = x,
      plot = FALSE,
      n = 1000
    )

    unif <- DHARMa::testUniformity(simresid, plot = FALSE)
    disp <- DHARMa::testDispersion(simresid, plot = FALSE)
    pdhess <- x$sdr$pdHess

    return(data.frame(
      unif = unif$p.value,
      disp = disp$p.value,
      pdhess = pdhess
    ))
  }

  testsummary_glmmtmb <- seqwrap::seqwrap(
    fitting_fun = glmmTMB::glmmTMB,
    arguments = list(
      formula = y ~
        time +
          (1 | participant),
      family = glmmTMB::nbinom2
    ),
    data = seqdatasubset,
    metadata = metadata,
    samplename = "seq_sample_id",
    additional_vars = NULL,
    summary_fun = summaryfun_glmmtmb,
    eval_fun = evalfun_glmmtmb,
    exported = list(),
    return_models = TRUE,
    save_models = FALSE,
    model_path = NULL,
    subset = NULL,
    cores = 1
  )

  ## Expected output from summary/evaluation functions
  expect_s3_class(testsummary_glmmtmb@summaries[[1]], "tbl")
  expect_s3_class(testsummary_glmmtmb@evaluations[[1]], "data.frame")

  ## Expect no errors in the error data
  expect_null(testsummary_glmmtmb@errors$err_sum[[1]])
  expect_null(testsummary_glmmtmb@errors$err_eval[[1]])

  # Expect errors when a bad function is passed
  #
  # To be written
  #
  #
})
