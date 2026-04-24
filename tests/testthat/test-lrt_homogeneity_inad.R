# test-test_homogeneity_inad.R

testthat::test_that("test_homogeneity_inad validates inputs correctly", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rpois(100, 5), nrow = 20, ncol = 5)
  blocks <- c(rep(1, 10), rep(2, 10))
  
  # blocks required
  testthat::expect_error(
    test_homogeneity_inad(y, order = 1),
    "blocks must be provided"
  )
  
  # blocks length must match
  testthat::expect_error(
    test_homogeneity_inad(y, blocks = c(1, 2, 3), order = 1),
    "length equal to nrow"
  )
  
  # Need at least 2 groups
  testthat::expect_error(
    test_homogeneity_inad(y, blocks = rep(1, 20), order = 1),
    "at least 2 groups"
  )
  
  # y must be non-negative integers
  y_bad <- y
  y_bad[1, 1] <- -1
  testthat::expect_error(
    test_homogeneity_inad(y_bad, blocks, order = 1),
    "non-negative integers"
  )
})


testthat::test_that("test_homogeneity_inad runs all test types", {
  skip_on_cran()
  set.seed(456)
  
  # Small simulated data
  n1 <- 25
  n2 <- 25
  n_time <- 5
  
  # Simulate from same model (homogeneous)
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Test all three test types
  test_all <- test_homogeneity_inad(y, blocks, order = 1, test = "all")
  testthat::expect_s3_class(test_all, "test_homogeneity_inad")
  testthat::expect_true(test_all$lrt_stat >= 0)
  testthat::expect_true(test_all$df > 0)
  testthat::expect_true(test_all$p_value >= 0 && test_all$p_value <= 1)
  testthat::expect_equal(test_all$test, "all")
  
  test_mean <- test_homogeneity_inad(y, blocks, order = 1, test = "mean")
  testthat::expect_s3_class(test_mean, "test_homogeneity_inad")
  testthat::expect_equal(test_mean$test, "mean")
  
  test_dep <- test_homogeneity_inad(y, blocks, order = 1, test = "dependence")
  testthat::expect_s3_class(test_dep, "test_homogeneity_inad")
  testthat::expect_equal(test_dep$test, "dependence")
})


testthat::test_that("test_homogeneity_inad detects heterogeneous groups", {
  skip_on_cran()
  set.seed(789)

  n1 <- 20
  n2 <- 20
  n_time <- 5

  # Simulate from DIFFERENT models
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.2, theta = 2)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.6, theta = 6)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Test for overall heterogeneity - should detect differences
  test_all <- test_homogeneity_inad(y, blocks, order = 1, test = "all")
  
  # With very different parameters, we expect to reject H0
  # (p_value should be small)
  testthat::expect_true(test_all$p_value < 0.10)  # Relaxed threshold for randomness
})


testthat::test_that("test_homogeneity_inad detects mean differences", {
  skip_on_cran()
  set.seed(101)

  n1 <- 25
  n2 <- 25
  n_time <- 5

  # Same dependence, different means
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.4, theta = 2)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.4, theta = 6)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Test for mean differences - should detect
  test_mean <- test_homogeneity_inad(y, blocks, order = 1, test = "mean")
  testthat::expect_true(test_mean$p_value < 0.10)
})


testthat::test_that("run_homogeneity_tests_inad works correctly", {
  skip_on_cran()
  set.seed(202)
  
  n1 <- 30
  n2 <- 30
  n_time <- 5
  
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 4)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 4)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Run all tests at once
  tests <- run_homogeneity_tests_inad(y, blocks, order = 1)
  
  testthat::expect_s3_class(tests, "homogeneity_tests_inad")
  testthat::expect_true("test_all" %in% names(tests))
  testthat::expect_true("test_mean" %in% names(tests))
  testthat::expect_true("test_dependence" %in% names(tests))
  testthat::expect_true("summary_table" %in% names(tests))
  testthat::expect_true("model_table" %in% names(tests))
  testthat::expect_true("best_model_bic" %in% names(tests))
  
  # Check summary table structure
  testthat::expect_equal(nrow(tests$summary_table), 3)
  testthat::expect_true(all(c("test", "lrt_stat", "df", "p_value") %in% 
                             names(tests$summary_table)))
  
  # Check model table structure
  testthat::expect_equal(nrow(tests$model_table), 3)
  testthat::expect_true(all(c("model", "logLik", "n_params", "BIC") %in% 
                             names(tests$model_table)))
})


testthat::test_that("test_homogeneity_inad works with order 0", {
  skip_on_cran()
  set.seed(303)
  
  n1 <- 25
  n2 <- 25
  n_time <- 5
  
  # Order 0 means independence
  y1 <- matrix(rpois(n1 * n_time, 3), nrow = n1, ncol = n_time)
  y2 <- matrix(rpois(n2 * n_time, 5), nrow = n2, ncol = n_time)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  test <- test_homogeneity_inad(y, blocks, order = 0, test = "all")
  
  testthat::expect_s3_class(test, "test_homogeneity_inad")
  testthat::expect_true(test$lrt_stat >= 0)
})


testthat::test_that("test_homogeneity_inad works with more than 2 groups", {
  skip_on_cran()
  set.seed(404)
  
  n1 <- 20
  n2 <- 25
  n3 <- 20
  n_time <- 5
  
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y3 <- simulate_inad(n3, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y <- rbind(y1, y2, y3)
  blocks <- c(rep(1, n1), rep(2, n2), rep(3, n3))
  
  test <- test_homogeneity_inad(y, blocks, order = 1, test = "all")
  
  testthat::expect_s3_class(test, "test_homogeneity_inad")
  testthat::expect_equal(test$settings$n_blocks, 3)
})


testthat::test_that("test_homogeneity_inad BIC selects simpler model when homogeneous", {
  skip_on_cran()
  set.seed(505)

  n1 <- 25
  n2 <- 25
  n_time <- 5

  # Same parameters (homogeneous case) - use scalar values
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.4, theta = 4)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom",
                      innovation = "pois", alpha = 0.4, theta = 4)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  tests <- run_homogeneity_tests_inad(y, blocks, order = 1)
  
  # BIC should prefer simpler model (pooled or INADFE) over heterogeneous
  # when data is truly homogeneous
  testthat::expect_true(tests$best_model_bic %in% c("M1 (Pooled)", "M2 (INADFE)"))
})


testthat::test_that("test_homogeneity_inad df is correct for innovation='nbinom'", {
  skip_on_cran()
  set.seed(808)

  n_time <- 4
  n1 <- 20; n2 <- 20

  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom",
                      innovation = "nbinom", alpha = 0.3, theta = 3,
                      nb_inno_size = 2)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom",
                      innovation = "nbinom", alpha = 0.3, theta = 3,
                      nb_inno_size = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))

  # order=1: n_alpha = n_time-1 = 3, n_theta = 4, n_nb = 4 per group
  # M1 k = 3+4+4 = 11;  M3 k = 2*(3+4)+2*4 = 22  ->  df = 11
  test_all <- test_homogeneity_inad(y, blocks, order = 1, test = "all",
                                    innovation = "nbinom",
                                    thinning = "binom", max_iter = 10)
  expected_df <- (2 - 1) * (3 + n_time + n_time)   # (n_blocks-1)*(n_alpha+n_theta+n_nb)
  testthat::expect_equal(test_all$df, expected_df)
})


testthat::test_that("print methods work without error", {
  skip_on_cran()
  set.seed(606)
  
  n1 <- 20
  n2 <- 20
  n_time <- 4
  
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Single test print
  test <- test_homogeneity_inad(y, blocks, order = 1, test = "mean")
  testthat::expect_output(print(test), "Homogeneity")
  
  # All tests print
  tests <- run_homogeneity_tests_inad(y, blocks, order = 1)
  testthat::expect_output(print(tests), "Homogeneity Tests")
})


testthat::test_that("test_homogeneity_inad uses pre-fitted models", {
  skip_on_cran()
  set.seed(707)

  n1 <- 15
  n2 <- 15
  n_time <- 4
  
  y1 <- simulate_inad(n1, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y2 <- simulate_inad(n2, n_time, order = 1, thinning = "binom", 
                      innovation = "pois", alpha = 0.3, theta = 3)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Pre-fit models
  fit_pooled <- fit_inad(y, order = 1, thinning = "binom", 
                         innovation = "pois", blocks = NULL)
  fit_inadfe <- fit_inad(y, order = 1, thinning = "binom", 
                         innovation = "pois", blocks = blocks)
  
  # Use pre-fitted pooled and inadfe
  test_mean <- test_homogeneity_inad(y, blocks, order = 1, test = "mean",
                                     fit_pooled = fit_pooled, 
                                     fit_inadfe = fit_inadfe)
  
  testthat::expect_s3_class(test_mean, "test_homogeneity_inad")
  testthat::expect_equal(test_mean$fit_null$log_l, fit_pooled$log_l)
  testthat::expect_equal(test_mean$fit_alt$log_l, fit_inadfe$log_l)
})
