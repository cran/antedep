test_that("aic_inad matches manual formula for INAD fit", {
  set.seed(411)
  y <- simulate_inad(
    n_subjects = 80,
    n_time = 6,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    alpha = 0.3,
    theta = 2
  )

  fit <- fit_inad(
    y = y,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    blocks = NULL
  )

  N <- ncol(y)
  k <- N + (N - 1)
  aic_manual <- -2 * fit$log_l + 2 * k

  expect_equal(aic_inad(fit), aic_manual)
})
