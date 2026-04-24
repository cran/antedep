testthat::test_that("test_homogeneity_gau runs on complete data", {
  skip_on_cran()
  set.seed(1201)
  y <- simulate_gau(
    n_subjects = 60,
    n_time = 6,
    order = 1,
    phi = 0.35,
    sigma = 1
  )
  blocks <- rep(1:2, each = 30)

  res <- test_homogeneity_gau(y, blocks = blocks, p = 1, use_modified = FALSE)

  testthat::expect_s3_class(res, "gau_homogeneity_test")
  testthat::expect_true(is.finite(res$statistic))
  testthat::expect_equal(res$df, (2 - 1) * (2 * ncol(y) - 1) * (1 + 1) / 2)
  testthat::expect_true(res$p_value >= 0 && res$p_value <= 1)
})

testthat::test_that("test_homogeneity_gau matches textbook cochlear example with dropout handling", {
  skip_on_cran()
  data("cochlear_implant", package = "antedep")

  res <- test_homogeneity_gau(
    cochlear_implant$y,
    blocks = cochlear_implant$blocks,
    p = 1,
    use_modified = FALSE
  )

  # Book example (speech recognition / cochlear): stat ~ 3.56, p ~ 0.83.
  testthat::expect_s3_class(res, "gau_homogeneity_test")
  testthat::expect_equal(unname(res$statistic), 3.560827, tolerance = 1e-4)
  testthat::expect_equal(unname(res$df), 7)
  testthat::expect_equal(unname(res$p_value), 0.8287414, tolerance = 1e-4)
})

testthat::test_that("test_homogeneity_gau rejects non-monotone missingness", {
  skip_on_cran()
  y <- matrix(
    c(1, 2, 3, 4,
      2, NA, 3, 4,
      1, 1, 2, NA,
      2, 2, 2, 2),
    nrow = 4,
    byrow = TRUE
  )
  blocks <- c(1, 1, 2, 2)

  testthat::expect_error(
    test_homogeneity_gau(y, blocks = blocks, p = 1),
    "monotone dropout"
  )
})
