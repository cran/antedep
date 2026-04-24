test_that("em_gau matches fit_gau with na_action = 'em'", {
  skip_on_cran()
  set.seed(321)
  y <- simulate_gau(
    n_subjects = 60,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.4, 4))
  )
  y[sample(length(y), 15)] <- NA

  fit_direct <- fit_gau(
    y = y,
    order = 1,
    na_action = "em",
    em_max_iter = 20,
    em_tol = 1e-6,
    em_verbose = FALSE
  )
  fit_wrapper <- em_gau(
    y = y,
    order = 1,
    max_iter = 20,
    tol = 1e-6,
    verbose = FALSE
  )

  expect_equal(fit_wrapper$log_l, fit_direct$log_l)
  expect_equal(fit_wrapper$em_iterations, fit_direct$em_iterations)
  expect_equal(fit_wrapper$em_converged, fit_direct$em_converged)
  expect_equal(fit_wrapper$settings$order, fit_direct$settings$order)
})
