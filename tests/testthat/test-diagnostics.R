# Tests for diagnostics.R

test_that("partial_corr returns expected structure", {
  data("bolus_inad")
  y <- bolus_inad$y
  
  result <- partial_corr(y)
  
  expect_type(result, "list")
  expect_true("partial_correlation" %in% names(result))
  expect_true(is.matrix(result$partial_correlation))
  expect_equal(nrow(result$partial_correlation), ncol(y))
  expect_equal(ncol(result$partial_correlation), ncol(y))
})

test_that("partial_corr with test = TRUE returns significance info", {
  data("bolus_inad")
  y <- bolus_inad$y
  
  result <- partial_corr(y, test = TRUE)
  
  expect_true("significant" %in% names(result))
  expect_true(is.matrix(result$significant))
})

test_that("partial_corr diagonal contains variances", {
  data("bolus_inad")
  y <- bolus_inad$y
  
  result <- partial_corr(y)
  
  diag_vals <- diag(result$partial_correlation)
  expect_true(all(diag_vals > 0))
})

test_that("partial_corr validates input", {
  # Test with single column matrix - should error
  expect_error(partial_corr(matrix(1:5, nrow = 5, ncol = 1)), "at least 2")
})

test_that("plot_prism runs without error", {
  skip_on_cran()
  data("bolus_inad")
  y <- bolus_inad$y[, 1:4]
  
  expect_silent(plot_prism(y))
})

test_that("plot_profile runs without error", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data("bolus_inad")
  y <- bolus_inad$y
  blocks <- bolus_inad$blocks
  
  p <- plot_profile(y, blocks = blocks)
  expect_s3_class(p, "ggplot")
})

test_that("plot_profile works without blocks", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  data("bolus_inad")
  y <- bolus_inad$y
  
  p <- plot_profile(y)
  expect_s3_class(p, "ggplot")
})
