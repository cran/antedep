# test-stationarity_ci_cat.R - Tests for stationarity testing and CIs

# ============================================================
# Test test_stationarity_cat
# ============================================================
test_that("test_stationarity_cat works for order 0", {
  skip_on_cran()
  set.seed(100)
  # Simulate with constant marginal (stationary under independence)
  y <- simulate_cat(80, 5, order = 0, n_categories = 2)
  
  test <- test_stationarity_cat(y, order = 0)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$order, 0)
  expect_true(test$p_value >= 0 && test$p_value <= 1)
})


test_that("test_stationarity_cat works for order 1", {
  skip_on_cran()
  set.seed(101)
  # Simulate with time-invariant transitions (should be approximately stationary)
  y <- simulate_cat(80, 5, order = 1, n_categories = 2)
  
  test <- test_stationarity_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$order, 1)
})


test_that("test_stationarity_cat detects non-stationarity", {
  skip_on_cran()
  set.seed(102)
  # Simulate with time-varying marginals AND transitions
  marg <- list(t1 = c(0.9, 0.1))  # Very skewed initial
  trans_varying <- list(
    t2 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    t4 = matrix(c(0.2, 0.8, 0.8, 0.2), 2, byrow = TRUE),
    t5 = matrix(c(0.6, 0.4, 0.4, 0.6), 2, byrow = TRUE)
  )

  y <- simulate_cat(400, 5, order = 1, n_categories = 2,
                    marginal = marg, transition = trans_varying)
  
  test <- test_stationarity_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  # Should reject stationarity
  expect_true(test$p_value < 0.05)
})


test_that("test_stationarity_cat df is correct", {
  skip_on_cran()
  set.seed(103)
  y <- simulate_cat(50, 5, order = 1, n_categories = 2)
  
  test <- test_stationarity_cat(y, order = 1)
  
  # For stationary model:
  # - params: (c-1) for marginal + (c-1)*c^p for transition = 1 + 2 = 3
  # For non-stationary model:
  # - params from fit_cat
  # df = non-stationary params - stationary params
  expect_true(test$df > 0)
  expect_equal(test$fit_null$n_params, 3)  # (2-1) + (2-1)*2^1 = 1 + 2 = 3
})


test_that("run_stationarity_tests_cat works", {
  skip_on_cran()
  set.seed(104)
  y <- simulate_cat(60, 5, order = 1, n_categories = 2)
  
  result <- run_stationarity_tests_cat(y, order = 1)
  
  expect_true(is.list(result))
  expect_true("time_invariance" %in% names(result))
  expect_true("stationarity" %in% names(result))
  expect_true("table" %in% names(result))
  
  expect_equal(nrow(result$table), 2)
  expect_true(all(c("test", "lrt_stat", "df", "p_value") %in% names(result$table)))
})


# ============================================================
# Test ci_cat
# ============================================================
test_that("ci_cat works for order 0", {
  skip_on_cran()
  set.seed(200)
  y <- simulate_cat(80, 4, order = 0, n_categories = 2)
  fit <- fit_cat(y, order = 0)
  
  ci <- ci_cat(fit)
  
  expect_s3_class(ci, "cat_ci")
  expect_equal(ci$level, 0.95)
  expect_true(!is.null(ci$marginal))
  expect_true(is.null(ci$transition))  # No transitions for order 0
  
  # Check marginal CIs
  expect_true(is.data.frame(ci$marginal))
  expect_true(all(c("parameter", "estimate", "se", "lower", "upper") %in% names(ci$marginal)))
  
  # CIs should be within [0, 1]
  expect_true(all(ci$marginal$lower >= 0, na.rm = TRUE))
  expect_true(all(ci$marginal$upper <= 1, na.rm = TRUE))
  expect_true(all(ci$marginal$lower <= ci$marginal$upper, na.rm = TRUE))
})


test_that("ci_cat works for order 1", {
  skip_on_cran()
  set.seed(201)
  y <- simulate_cat(80, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  
  ci <- ci_cat(fit)
  
  expect_s3_class(ci, "cat_ci")
  expect_true(!is.null(ci$marginal))
  expect_true(!is.null(ci$transition))
  
  # Check transition CIs
  expect_true(is.list(ci$transition))
  expect_true(length(ci$transition) > 0)
  
  # Check first transition time point
  trans_df <- ci$transition[[1]]
  expect_true(is.data.frame(trans_df))
  expect_true(all(c("parameter", "estimate", "se", "lower", "upper") %in% names(trans_df)))
})


test_that("ci_cat works for order 2", {
  skip_on_cran()
  set.seed(202)
  y <- simulate_cat(100, 5, order = 2, n_categories = 2)
  fit <- fit_cat(y, order = 2)
  
  ci <- ci_cat(fit)
  
  expect_s3_class(ci, "cat_ci")
  expect_true(!is.null(ci$transition))
  
  # For order 2, transition CIs should have from1, from2, to columns
  trans_df <- ci$transition[[1]]
  expect_true(all(c("from1", "from2", "to") %in% names(trans_df)))
})


test_that("ci_cat respects level parameter", {
  skip_on_cran()
  set.seed(203)
  y <- simulate_cat(80, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  
  ci_95 <- ci_cat(fit, level = 0.95)
  ci_99 <- ci_cat(fit, level = 0.99)
  
  expect_equal(ci_95$level, 0.95)
  expect_equal(ci_99$level, 0.99)
  
  # 99% CIs should be wider than 95% CIs
  width_95 <- ci_95$marginal$upper - ci_95$marginal$lower
  width_99 <- ci_99$marginal$upper - ci_99$marginal$lower
  expect_true(all(width_99 >= width_95 - 1e-10))  # Allow tiny numerical tolerance
})


test_that("ci_cat parameters argument works", {
  skip_on_cran()
  set.seed(204)
  y <- simulate_cat(80, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  
  ci_marg <- ci_cat(fit, parameters = "marginal")
  ci_trans <- ci_cat(fit, parameters = "transition")
  ci_all <- ci_cat(fit, parameters = "all")
  
  expect_true(!is.null(ci_marg$marginal))
  expect_true(is.null(ci_marg$transition))
  
  expect_true(is.null(ci_trans$marginal))
  expect_true(!is.null(ci_trans$transition))
  
  expect_true(!is.null(ci_all$marginal))
  expect_true(!is.null(ci_all$transition))
})


test_that("ci_cat works with heterogeneous model", {
  skip_on_cran()
  set.seed(205)
  y1 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y2 <- simulate_cat(50, 3, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 50), rep(2, 50))
  
  fit <- fit_cat(y, order = 1, blocks = blocks, homogeneous = FALSE)
  ci <- ci_cat(fit)
  
  expect_s3_class(ci, "cat_ci")
  
  # Should have per-block CIs
  expect_true(is.list(ci$marginal))
  expect_true("block_1" %in% names(ci$marginal))
  expect_true("block_2" %in% names(ci$marginal))
})


test_that("ci_cat validates inputs", {
  skip_on_cran()
  set.seed(206)
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  
  # Invalid level
  expect_error(ci_cat(fit, level = 0))
  expect_error(ci_cat(fit, level = 1))
  expect_error(ci_cat(fit, level = 1.5))
  
  # Invalid fit
  expect_error(ci_cat("not a fit"))
})


test_that("ci_cat covers true values", {
  skip_on_cran()
  set.seed(207)
  # Simulate with known parameters
  true_marg <- c(0.6, 0.4)
  true_trans <- matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)

  marginal <- list(t1 = true_marg)
  transition <- list(
    t2 = true_trans,
    t3 = true_trans,
    t4 = true_trans
  )

  y <- simulate_cat(1000, 4, order = 1, n_categories = 2,
                    marginal = marginal, transition = transition)
  fit <- fit_cat(y, order = 1)
  ci <- ci_cat(fit, level = 0.95)
  
  # Check that true marginal values are covered by CIs
  marg_ci <- ci$marginal
  for (i in 1:2) {
    row_i <- marg_ci[marg_ci$category == i, ]
    expect_true(row_i$lower <= true_marg[i] && true_marg[i] <= row_i$upper,
                info = paste("Marginal category", i, "not covered"))
  }
  
  # Check that true transition values are covered (at least for t2)
  trans_ci <- ci$transition[["t2"]]
  for (i in 1:2) {
    for (j in 1:2) {
      row_ij <- trans_ci[trans_ci$from == i & trans_ci$to == j, ]
      expect_true(row_ij$lower <= true_trans[i, j] && true_trans[i, j] <= row_ij$upper,
                  info = paste("Transition", i, "->", j, "not covered"))
    }
  }
})


test_that("print.cat_ci works", {
  skip_on_cran()
  set.seed(208)
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  ci <- ci_cat(fit)
  
  # Just check it doesn't error
  expect_output(print(ci), "Confidence Intervals")
})


test_that("summary.cat_ci works", {
  skip_on_cran()
  set.seed(209)
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  ci <- ci_cat(fit)
  
  summ <- summary(ci)
  
  expect_true(is.data.frame(summ))
  expect_true(all(c("parameter", "type", "estimate", "lower", "upper") %in% names(summ)))
})
