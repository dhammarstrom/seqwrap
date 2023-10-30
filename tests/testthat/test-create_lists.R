

## Load required data
seqdata <- seqwrap::rna_seq_sample


test_that("create_lists returns a list", {
  expect_type(seqwrap::create_lists(seqdata), "list")
})

test_that("create_lists returns a list of data frames", {
  tmp <- create_lists(seqdata)
  expect_s3_class(tmp[[1]], "data.frame")
})


test_that("create_lists returns an error if the first column in seqdata is not character or factor", {
  seqdata <- seqdata[,-1]
  expect_error(create_lists(seqdata))
})
