test_that("ci_gau returns Wald intervals for complete-data Gaussian AD fit", {
  skip_on_cran()
  set.seed(510)
  y <- simulate_gau(
    n_subjects = 90,
    n_time = 6,
    order = 1,
    phi = c(0, rep(0.4, 5))
  )
  fit <- fit_gau(y, order = 1)

  ci <- ci_gau(fit, level = 0.95, parameters = "all")

  expect_s3_class(ci, "gau_ci")
  expect_true(is.data.frame(ci$mu))
  expect_true(is.data.frame(ci$phi))
  expect_true(is.data.frame(ci$sigma))
  expect_equal(nrow(ci$mu), ncol(y))
  expect_equal(nrow(ci$sigma), ncol(y))
  expect_equal(nrow(ci$phi), ncol(y) - 1)
})

test_that("ci_gau rejects missing-data Gaussian AD fits", {
  skip_on_cran()
  set.seed(511)
  y <- simulate_gau(
    n_subjects = 60,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.35, 4))
  )
  y[sample(length(y), 12)] <- NA
  fit_em <- fit_gau(y, order = 1, na_action = "em", em_max_iter = 20)

  expect_error(
    ci_gau(fit_em),
    "complete-data fits only",
    fixed = FALSE
  )
})
