test_that("CAT missing-data inference support boundaries are explicit", {
  skip_on_cran()
  y_miss <- matrix(
    c(1, 2, 1,
      2, 1, NA,
      1, 2, 2,
      2, 1, 1),
    nrow = 4,
    byrow = TRUE
  )

  fit_missing_ci <- structure(
    list(
      settings = list(
        order = 1,
        na_action = "marginalize",
        na_action_effective = "marginalize"
      ),
      cell_counts = NULL
    ),
    class = "cat_fit"
  )
  expect_error(ci_cat(fit_missing_ci), "complete-data fits only")

  expect_s3_class(test_order_cat(y_miss, order_null = 0, order_alt = 1), "cat_lrt")
  expect_s3_class(
    test_homogeneity_cat(y_miss, blocks = c(1, 1, 2, 2), order = 1),
    "cat_lrt"
  )
  expect_true(is.list(run_order_tests_cat(y_miss)))

  expect_error(test_timeinvariance_cat(y_miss, order = 1), "complete data only")
  expect_error(test_stationarity_cat(y_miss, order = 1), "complete data only")
  expect_error(run_stationarity_tests_cat(y_miss, order = 1), "complete data only")

  fit_null <- structure(
    list(
      settings = list(order = 0, homogeneous = TRUE, n_blocks = 1,
                      na_action_effective = "marginalize"),
      cell_counts = NULL,
      log_l = -10,
      n_params = 2,
      aic = 24,
      bic = 25
    ),
    class = "cat_fit"
  )
  fit_alt <- structure(
    list(
      settings = list(order = 1, homogeneous = TRUE, n_blocks = 1,
                      na_action_effective = "marginalize"),
      cell_counts = NULL,
      log_l = -9,
      n_params = 4,
      aic = 26,
      bic = 27
    ),
    class = "cat_fit"
  )
  expect_s3_class(
    test_order_cat(fit_null = fit_null, fit_alt = fit_alt),
    "cat_lrt"
  )
})

test_that("INAD missing-data LRT support is enabled while CI remains complete-data only", {
  skip_on_cran()
  set.seed(2026)
  y_miss <- simulate_inad(
    n_subjects = 10,
    n_time = 3,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    alpha = 0.35,
    theta = 2.5
  )
  y_miss[sample(length(y_miss), 4)] <- NA
  blocks <- rep(1:2, each = 5)

  fit_stub <- list(
    settings = list(order = 1, thinning = "binom", innovation = "pois"),
    theta = rep(1, ncol(y_miss)),
    log_l = -10
  )

  expect_error(ci_inad(y_miss, fit_stub, blocks = blocks), "complete data only")
  ord_test <- suppressWarnings(
    test_order_inad(y_miss, order_null = 0, order_alt = 1, blocks = blocks)
  )
  expect_s3_class(ord_test, "test_order_inad")
  expect_s3_class(
    test_stationarity_inad(y_miss, order = 1, blocks = blocks),
    "test_stationarity_inad"
  )
  expect_s3_class(
    test_homogeneity_inad(y_miss, blocks = blocks, order = 1),
    "test_homogeneity_inad"
  )
  stat_tests <- suppressWarnings(
    run_stationarity_tests_inad(y_miss, order = 1, blocks = blocks, verbose = FALSE)
  )
  expect_s3_class(stat_tests, "stationarity_tests_inad")
  expect_s3_class(
    run_homogeneity_tests_inad(y_miss, blocks = blocks, order = 1),
    "homogeneity_tests_inad"
  )
})

test_that("AD likelihood-ratio utilities reject missing-data inference clearly", {
  skip_on_cran()
  y_miss <- matrix(
    c(1.0, 2.0, 3.0,
      2.0, 1.0, 0.0,
      0.5, NA, 1.5,
      3.0, 2.0, 1.0),
    nrow = 4,
    byrow = TRUE
  )
  blocks <- c(1, 1, 2, 2)

  expect_error(test_order_gau(y_miss, p = 0, q = 1), "complete data only")
  expect_error(test_homogeneity_gau(y_miss, blocks = blocks, p = 1), "monotone dropout")
  expect_error(
    test_one_sample_gau(y_miss, mu0 = c(0, 0, 0), p = 1),
    "complete data only"
  )
  expect_error(test_two_sample_gau(y_miss, blocks = blocks, p = 1), "complete data only")
  expect_error(
    test_contrast_gau(y_miss, C = matrix(c(1, -1, 0), nrow = 1), p = 1),
    "complete data only"
  )
  expect_error(test_stationarity_gau(y_miss, order = 1), "complete data only")
  expect_error(run_stationarity_tests_gau(y_miss, order = 1), "complete data only")
})
