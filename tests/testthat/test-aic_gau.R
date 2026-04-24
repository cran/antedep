test_that("aic_gau matches manual formula for Gaussian AD fit", {
  set.seed(410)
  y <- simulate_gau(
    n_subjects = 30,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.35, 4))
  )

  fit <- fit_gau(y, order = 1, estimate_mu = TRUE)
  aic_pkg <- aic_gau(fit)

  N <- ncol(y)
  k <- N + N + (N - 1)
  aic_manual <- -2 * fit$log_l + 2 * k

  expect_equal(aic_pkg, aic_manual)
})
