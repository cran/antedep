# test-lrt_cat.R - Tests for CAT hypothesis testing functions

# ============================================================
# Test test_order_cat
# ============================================================
test_that("test_order_cat works for AD(0) vs AD(1)", {
  skip_on_cran()
  set.seed(123)
  # Simulate AD(1) data - should reject AD(0)
  marg <- list(t1 = c(0.6, 0.4))
  trans <- list(
    t2 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE),
    t3 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE),
    t4 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
  )
  y <- simulate_cat(100, 4, order = 1, n_categories = 2,
                    marginal = marg, transition = trans)

  test <- test_order_cat(y, order_null = 0, order_alt = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_true(test$p_value >= 0 && test$p_value <= 1)
  expect_equal(test$order_null, 0)
  expect_equal(test$order_alt, 1)
  
  # For data with true AD(1), should reject AD(0)
  expect_true(test$p_value < 0.05)
})


test_that("test_order_cat works for AD(1) vs AD(2)", {
  skip_on_cran()
  set.seed(456)
  # Simulate AD(1) data - should NOT reject AD(1) in favor of AD(2)
  y <- simulate_cat(80, 5, order = 1, n_categories = 2)

  test <- test_order_cat(y, order_null = 1, order_alt = 2)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_equal(test$order_null, 1)
  expect_equal(test$order_alt, 2)
  
  # For data with true AD(1), should usually NOT reject AD(1)
  # (but with random data this isn't guaranteed, so just check it runs)
})

test_that("test_order_cat supports score and mlrt options", {
  skip_on_cran()
  set.seed(460)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)

  score_test <- test_order_cat(y, order_null = 0, order_alt = 1, test = "score")
  expect_s3_class(score_test, "cat_lrt")
  expect_equal(score_test$test, "score")
  expect_true(score_test$lrt_stat >= 0)
  expect_true(score_test$p_value >= 0 && score_test$p_value <= 1)

  mlrt_test <- test_order_cat(y, order_null = 0, order_alt = 1, test = "mlrt")
  expect_s3_class(mlrt_test, "cat_lrt")
  expect_equal(mlrt_test$test, "mlrt")
  expect_true(is.finite(mlrt_test$e_hat_mlrt))
  expect_true(mlrt_test$e_hat_mlrt > 0)
  expect_true(mlrt_test$lrt_stat >= 0)
  expect_true(mlrt_test$p_value >= 0 && mlrt_test$p_value <= 1)

  lrt_test <- test_order_cat(y, order_null = 0, order_alt = 1, test = "lrt")
  expect_false(isTRUE(all.equal(mlrt_test$lrt_stat, lrt_test$lrt_stat, tolerance = 1e-12)))
})

test_that("test_order_cat supports wald option", {
  skip_on_cran()
  set.seed(461)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)

  wald_test <- test_order_cat(y, order_null = 0, order_alt = 1, test = "wald")
  expect_s3_class(wald_test, "cat_lrt")
  expect_equal(wald_test$test, "wald")
  expect_true(is.finite(wald_test$lrt_stat))
  expect_true(wald_test$lrt_stat >= 0)
  expect_true(wald_test$p_value >= 0 && wald_test$p_value <= 1)
})


test_that("test_order_cat works with pre-fitted models", {
  skip_on_cran()
  set.seed(789)
  y <- simulate_cat(60, 4, order = 1, n_categories = 2)
  
  fit0 <- fit_cat(y, order = 0)
  fit1 <- fit_cat(y, order = 1)
  
  # Test with pre-fitted models
  test <- test_order_cat(fit_null = fit0, fit_alt = fit1)
  
  expect_s3_class(test, "cat_lrt")
  expect_equal(test$order_null, 0)
  expect_equal(test$order_alt, 1)
  
  # Test with y + one pre-fitted model
  test2 <- test_order_cat(y, order_null = 0, fit_alt = fit1)
  expect_equal(test$lrt_stat, test2$lrt_stat, tolerance = 1e-10)
})


test_that("test_order_cat validates inputs correctly", {
  skip_on_cran()
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  
  # Error: order_alt <= order_null
  expect_error(test_order_cat(y, order_null = 1, order_alt = 0))
  expect_error(test_order_cat(y, order_null = 1, order_alt = 1))
  
  # Error: negative order
  expect_error(test_order_cat(y, order_null = -1, order_alt = 1))
  
  # Error: order too high
  expect_error(test_order_cat(y, order_null = 2, order_alt = 3))
  
  # Error: no data and no pre-fitted models
  expect_error(test_order_cat(order_null = 0, order_alt = 1))
})


test_that("test_order_cat df calculation is correct", {
  skip_on_cran()
  set.seed(111)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)
  
  # AD(0) vs AD(1) with c=2, n=5
  # df = (c-1) * c^0 * (n-1) = 1 * 1 * 4 = 4
  test_01 <- test_order_cat(y, order_null = 0, order_alt = 1)
  expect_equal(test_01$df, 4)
  
  # AD(1) vs AD(2) with c=2, n=5
  # df = (c-1) * c^1 * (n-2) = 1 * 2 * 3 = 6
  test_12 <- test_order_cat(y, order_null = 1, order_alt = 2)
  expect_equal(test_12$df, 6)
})


test_that("run_order_tests_cat works correctly", {
  skip_on_cran()
  set.seed(222)
  # Simulate AD(1) data
  marg <- list(t1 = c(0.6, 0.4))
  trans <- list(
    t2 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t3 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t4 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t5 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE)
  )
  y <- simulate_cat(120, 5, order = 1, n_categories = 2,
                    marginal = marg, transition = trans)

  result <- run_order_tests_cat(y, max_order = 2)
  
  expect_true(is.list(result))
  expect_true("tests" %in% names(result))
  expect_true("table" %in% names(result))
  expect_true("fits" %in% names(result))
  expect_true("selected_order" %in% names(result))
  
  expect_equal(nrow(result$table), 2)  # Two tests: 0v1, 1v2
  expect_true(result$selected_order %in% 0:2)
  
  # For true AD(1), should select order 1
  # (0v1 significant, 1v2 not significant)
  expect_equal(result$selected_order, 1)
})


# ============================================================
# Test test_homogeneity_cat
# ============================================================
test_that("test_homogeneity_cat detects heterogeneity", {
  skip_on_cran()
  set.seed(333)
  # Create two groups with different transition probabilities
  marg1 <- list(t1 = c(0.8, 0.2))
  marg2 <- list(t1 = c(0.3, 0.7))
  trans1 <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE)
  )
  trans2 <- list(
    t2 = matrix(c(0.4, 0.6, 0.6, 0.4), 2, byrow = TRUE),
    t3 = matrix(c(0.4, 0.6, 0.6, 0.4), 2, byrow = TRUE)
  )
  
  y1 <- simulate_cat(80, 3, order = 1, n_categories = 2,
                     marginal = marg1, transition = trans1)
  y2 <- simulate_cat(80, 3, order = 1, n_categories = 2,
                     marginal = marg2, transition = trans2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 80), rep(2, 80))
  
  test <- test_homogeneity_cat(y, blocks, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$n_groups, 2)
  
  # Should reject homogeneity
  expect_true(test$p_value < 0.05)
})


test_that("test_homogeneity_cat accepts homogeneous data", {
  skip_on_cran()
  set.seed(444)
  # Create two groups with SAME parameters
  y1 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y2 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 50), rep(2, 50))
  
  test <- test_homogeneity_cat(y, blocks, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  
  # Should NOT reject homogeneity (p > 0.05 usually, but not guaranteed)
  # Just check that test runs and gives reasonable output
  expect_true(test$p_value >= 0 && test$p_value <= 1)
})

test_that("test_homogeneity_cat supports score option", {
  skip_on_cran()
  set.seed(445)
  y1 <- simulate_cat(50, 4, order = 1, n_categories = 2)
  y2 <- simulate_cat(50, 4, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 50), rep(2, 50))

  score_test <- test_homogeneity_cat(y, blocks, order = 1, test = "score")
  expect_s3_class(score_test, "cat_lrt")
  expect_equal(score_test$test, "score")
  expect_true(score_test$lrt_stat >= 0)
  expect_true(score_test$p_value >= 0 && score_test$p_value <= 1)
})

test_that("test_homogeneity_cat supports mlrt option", {
  skip_on_cran()
  old_opt <- options(antedep.cat_mlrt_nsim = 20L, antedep.cat_mlrt_seed = 101L)
  on.exit(options(old_opt), add = TRUE)

  set.seed(446)
  y1 <- simulate_cat(50, 4, order = 1, n_categories = 2)
  y2 <- simulate_cat(50, 4, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 50), rep(2, 50))

  mlrt_test <- test_homogeneity_cat(y, blocks, order = 1, test = "mlrt")
  expect_s3_class(mlrt_test, "cat_lrt")
  expect_equal(mlrt_test$test, "mlrt")
  expect_true(is.finite(mlrt_test$e_hat_mlrt))
  expect_true(mlrt_test$e_hat_mlrt > 0)
  expect_true(mlrt_test$lrt_stat >= 0)
})


test_that("test_homogeneity_cat works with pre-fitted models", {
  skip_on_cran()
  set.seed(555)
  y1 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y2 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 50), rep(2, 50))
  
  fit_homo <- fit_cat(y, order = 1, blocks = blocks, homogeneous = TRUE)
  fit_hetero <- fit_cat(y, order = 1, blocks = blocks, homogeneous = FALSE)
  
  test <- test_homogeneity_cat(fit_null = fit_homo, fit_alt = fit_hetero)
  
  expect_s3_class(test, "cat_lrt")
  expect_equal(test$n_groups, 2)
})


test_that("test_homogeneity_cat validates inputs", {
  skip_on_cran()
  y <- simulate_cat(50, 3, order = 1, n_categories = 2)
  
  # Error: no blocks provided
  expect_error(test_homogeneity_cat(y, blocks = NULL))
  
  # Error: no y and no pre-fitted models
  expect_error(test_homogeneity_cat(blocks = c(1, 1, 2, 2)))
})


# ============================================================
# Test test_timeinvariance_cat
# ============================================================
test_that("test_timeinvariance_cat works", {
  skip_on_cran()
  set.seed(666)
  # Simulate with time-invariant transitions (default)
  y <- simulate_cat(80, 5, order = 1, n_categories = 2)
  
  test <- test_timeinvariance_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$order, 1)
  
  # For time-invariant data, should NOT reject time-invariance
  # (but this is statistical, so just check it runs)
  expect_true(test$p_value >= 0 && test$p_value <= 1)
})


test_that("test_timeinvariance_cat detects time-varying transitions", {
  skip_on_cran()
  set.seed(777)
  # Simulate with time-VARYING transitions
  marg <- list(t1 = c(0.5, 0.5))
  trans_varying <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
    t4 = matrix(c(0.2, 0.8, 0.8, 0.2), 2, byrow = TRUE),
    t5 = matrix(c(0.7, 0.3, 0.3, 0.7), 2, byrow = TRUE)
  )
  
  y <- simulate_cat(200, 5, order = 1, n_categories = 2,
                    marginal = marg, transition = trans_varying)

  test <- test_timeinvariance_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  
  # Should reject time-invariance
  expect_true(test$p_value < 0.05)
})

test_that("test_timeinvariance_cat supports score option", {
  skip_on_cran()
  set.seed(778)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)
  score_test <- test_timeinvariance_cat(y, order = 1, test = "score")
  expect_s3_class(score_test, "cat_lrt")
  expect_equal(score_test$test, "score")
  expect_true(score_test$lrt_stat >= 0)
})

test_that("test_timeinvariance_cat supports mlrt option", {
  skip_on_cran()
  old_opt <- options(antedep.cat_mlrt_nsim = 20L, antedep.cat_mlrt_seed = 202L)
  on.exit(options(old_opt), add = TRUE)

  set.seed(779)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)
  mlrt_test <- test_timeinvariance_cat(y, order = 1, test = "mlrt")
  expect_s3_class(mlrt_test, "cat_lrt")
  expect_equal(mlrt_test$test, "mlrt")
  expect_true(is.finite(mlrt_test$e_hat_mlrt))
  expect_true(mlrt_test$e_hat_mlrt > 0)
  expect_true(mlrt_test$lrt_stat >= 0)
})


test_that("test_timeinvariance_cat validates inputs", {
  skip_on_cran()
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  
  # Error: order < 1
  expect_error(test_timeinvariance_cat(y, order = 0))
})


# ============================================================
# Test print methods
# ============================================================
test_that("print.cat_lrt works", {
  skip_on_cran()
  set.seed(888)
  y <- simulate_cat(100, 4, order = 1, n_categories = 2)
  
  test <- test_order_cat(y, order_null = 0, order_alt = 1)
  
  # Just verify it doesn't error
  expect_output(print(test), "Likelihood Ratio Test")
})

test_that("test_stationarity_cat supports score option (order 1)", {
  skip_on_cran()
  set.seed(889)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)
  score_test <- test_stationarity_cat(y, order = 1, test = "score")
  expect_s3_class(score_test, "cat_lrt")
  expect_equal(score_test$test, "score")
  expect_true(score_test$lrt_stat >= 0)
  expect_true(score_test$p_value >= 0 && score_test$p_value <= 1)
})

test_that("test_stationarity_cat supports score option (order 2)", {
  skip_on_cran()
  set.seed(890)
  y <- simulate_cat(100, 6, order = 2, n_categories = 2)
  score_test <- NULL
  expect_warning(
    score_test <- test_stationarity_cat(y, order = 2, test = "score"),
    "marginal-constancy plus time-invariant transitions"
  )
  expect_s3_class(score_test, "cat_lrt")
  expect_equal(score_test$test, "score")
  expect_true(score_test$lrt_stat >= 0)
  expect_true(score_test$p_value >= 0 && score_test$p_value <= 1)
})

test_that("test_stationarity_cat blocks lrt/mlrt at order 2", {
  skip_on_cran()
  set.seed(892)
  y <- simulate_cat(80, 6, order = 2, n_categories = 2)
  expect_error(
    test_stationarity_cat(y, order = 2, test = "lrt"),
    "not supported for order >= 2"
  )
  expect_error(
    test_stationarity_cat(y, order = 2, test = "mlrt"),
    "not supported for order >= 2"
  )
})

test_that("test_stationarity_cat supports mlrt option (order 1)", {
  skip_on_cran()
  old_opt <- options(antedep.cat_mlrt_nsim = 20L, antedep.cat_mlrt_seed = 303L)
  on.exit(options(old_opt), add = TRUE)

  set.seed(891)
  y <- simulate_cat(80, 5, order = 1, n_categories = 2)
  mlrt_test <- test_stationarity_cat(y, order = 1, test = "mlrt")
  expect_s3_class(mlrt_test, "cat_lrt")
  expect_equal(mlrt_test$test, "mlrt")
  expect_true(is.finite(mlrt_test$e_hat_mlrt))
  expect_true(mlrt_test$e_hat_mlrt > 0)
  expect_true(mlrt_test$lrt_stat >= 0)
  expect_true(mlrt_test$p_value >= 0 && mlrt_test$p_value <= 1)
})
