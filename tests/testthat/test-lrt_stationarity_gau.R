test_that("test_stationarity_gau returns expected df for order-1 constraints", {
  skip_on_cran()
  set.seed(1201)
  y <- simulate_gau(
    n_subjects = 40,
    n_time = 5,
    order = 1,
    mu = 0,
    phi = 0.45,
    sigma = 1
  )

  fit_uncon <- fit_gau(y, order = 1, na_action = "fail")

  res_phi <- test_stationarity_gau(
    y,
    order = 1,
    constrain = "phi",
    fit_unconstrained = fit_uncon,
    max_iter = 700
  )
  res_sigma <- test_stationarity_gau(
    y,
    order = 1,
    constrain = "sigma",
    fit_unconstrained = fit_uncon,
    max_iter = 700
  )
  res_both <- test_stationarity_gau(
    y,
    order = 1,
    constrain = "both",
    fit_unconstrained = fit_uncon,
    max_iter = 700
  )

  expect_s3_class(res_phi, "test_stationarity_gau")
  expect_s3_class(res_sigma, "test_stationarity_gau")
  expect_s3_class(res_both, "test_stationarity_gau")

  expect_equal(res_phi$df, ncol(y) - 2)
  expect_equal(res_sigma$df, ncol(y) - 1)
  expect_equal(res_both$df, (ncol(y) - 2) + (ncol(y) - 1))

  expect_true(is.finite(res_phi$lrt_stat))
  expect_true(is.finite(res_sigma$lrt_stat))
  expect_true(is.finite(res_both$lrt_stat))
  expect_true(is.finite(res_phi$p_value))
  expect_true(is.finite(res_sigma$p_value))
  expect_true(is.finite(res_both$p_value))
})


test_that("test_stationarity_gau detects strong order-1 non-stationarity", {
  skip_on_cran()
  set.seed(1202)
  y <- simulate_gau(
    n_subjects = 280,
    n_time = 6,
    order = 1,
    mu = 0,
    phi = c(0, 0.15, 0.85, -0.2, 0.8, -0.15),
    sigma = c(1.0, 1.8, 0.6, 1.7, 0.55, 1.6)
  )

  res <- test_stationarity_gau(
    y,
    order = 1,
    constrain = "both",
    max_iter = 1000
  )

  expect_s3_class(res, "test_stationarity_gau")
  expect_lt(res$p_value, 0.05)
})


test_that("run_stationarity_tests_gau returns expected order-2 battery", {
  skip_on_cran()
  set.seed(1203)
  phi <- matrix(0, nrow = 2, ncol = 6)
  phi[1, 2:6] <- 0.4
  phi[2, 3:6] <- 0.2

  y <- simulate_gau(
    n_subjects = 30,
    n_time = 6,
    order = 2,
    mu = 0,
    phi = phi,
    sigma = 1
  )

  out <- run_stationarity_tests_gau(y, order = 2, verbose = FALSE, max_iter = 700)

  expect_s3_class(out, "stationarity_tests_gau")
  expect_setequal(out$summary$constraint, c("phi1", "phi2", "phi", "sigma", "all"))
  expect_true(all(out$summary$df > 0))
  expect_equal(length(out$tests), 5)
})
