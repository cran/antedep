# test-cat.R - Tests for categorical antedependence module

# ============================================================
# Test 1: Basic validation
# ============================================================
test_that("input validation works correctly", {
  skip_on_cran()
  # Valid data
  y_valid <- matrix(c(1, 2, 1, 2,
                      2, 1, 1, 1,
                      1, 1, 2, 2), nrow = 3, byrow = TRUE)
  
  validated <- antedep:::.validate_y_cat(y_valid)
  expect_equal(validated$n_categories, 2)
  expect_equal(dim(validated$y), c(3, 4))
  
  # Override n_categories
  validated2 <- antedep:::.validate_y_cat(y_valid, n_categories = 3)
  expect_equal(validated2$n_categories, 3)
  
  # Invalid: non-integer
  expect_error(antedep:::.validate_y_cat(matrix(c(1.5, 2), nrow = 1)))
  
  # Invalid: values < 1
  expect_error(antedep:::.validate_y_cat(matrix(c(0, 1, 2), nrow = 1)))
})


test_that("parameter counting is correct", {
  skip_on_cran()
  # Order 0: (c-1) * n = 1 * 4 = 4
  n_params_0 <- antedep:::.count_params_cat(order = 0, n_categories = 2, n_time = 4)
  expect_equal(n_params_0, 4)
  
  # Order 1: (c-1) * [1 + (n-1)*c] = 1 * [1 + 3*2] = 7
  n_params_1 <- antedep:::.count_params_cat(order = 1, n_categories = 2, n_time = 4)
  expect_equal(n_params_1, 7)
  
  # Order 2: (c-1) * [1 + c + (n-2)*c^2] = 1 * [1 + 2 + 2*4] = 11
  n_params_2 <- antedep:::.count_params_cat(order = 2, n_categories = 2, n_time = 4)
  expect_equal(n_params_2, 11)
  
  # With 3 categories, order 1, 5 time points
  # (c-1) * [1 + (n-1)*c] = 2 * [1 + 4*3] = 2 * 13 = 26
  n_params_3cat <- antedep:::.count_params_cat(order = 1, n_categories = 3, n_time = 5)
  expect_equal(n_params_3cat, 26)
})


# ============================================================
# Test 2: Cell counting
# ============================================================
test_that("cell counting works correctly", {
  skip_on_cran()
  y_test <- matrix(c(1, 1, 2,
                     1, 2, 1,
                     2, 1, 1,
                     1, 1, 1), nrow = 4, byrow = TRUE)
  
  # Count at time 1
  counts_t1 <- antedep:::.count_cells_table_cat(y_test, 1, n_categories = 2)
  expect_equal(counts_t1[1], 3)  # Three 1's
  expect_equal(counts_t1[2], 1)  # One 2
  
  # Count pairs (t1, t2)
  counts_t12 <- antedep:::.count_cells_table_cat(y_test, c(1, 2), n_categories = 2)
  # (1,1): subjects 1, 4 -> 2
  # (1,2): subject 2 -> 1
  # (2,1): subject 3 -> 1
  # (2,2): none -> 0
  expect_equal(counts_t12[1, 1], 2)
  expect_equal(counts_t12[1, 2], 1)
  expect_equal(counts_t12[2, 1], 1)
  expect_equal(counts_t12[2, 2], 0)
})


# ============================================================
# Test 3: Transition probability computation
# ============================================================
test_that("transition probability computation is correct", {
  skip_on_cran()
  y_test <- matrix(c(1, 1, 2,
                     1, 2, 1,
                     2, 1, 1,
                     1, 1, 1), nrow = 4, byrow = TRUE)
  
  counts_t12 <- antedep:::.count_cells_table_cat(y_test, c(1, 2), n_categories = 2)
  trans_probs <- antedep:::.counts_to_transition_probs(counts_t12)
  
  # P(Y2=1|Y1=1) = 2/3, P(Y2=2|Y1=1) = 1/3
  # P(Y2=1|Y1=2) = 1/1 = 1, P(Y2=2|Y1=2) = 0/1 = 0
  expect_equal(trans_probs[1, 1], 2/3, tolerance = 1e-10)
  expect_equal(trans_probs[1, 2], 1/3, tolerance = 1e-10)
  expect_equal(trans_probs[2, 1], 1, tolerance = 1e-10)
  expect_equal(trans_probs[2, 2], 0, tolerance = 1e-10)
})


# ============================================================
# Test 4: fit_cat basic functionality
# ============================================================
test_that("fit_cat works for order 0", {
  skip_on_cran()
  set.seed(42)
  y_sim <- matrix(sample(1:2, 100 * 5, replace = TRUE), nrow = 100, ncol = 5)
  
  fit0 <- fit_cat(y_sim, order = 0)
  expect_s3_class(fit0, "cat_fit")
  expect_true(is.finite(fit0$log_l))
  expect_true(fit0$log_l < 0)
  expect_equal(fit0$n_params, 5)  # (c-1)*n = 1*5 = 5
  expect_equal(fit0$n_obs, length(y_sim))
  expect_equal(fit0$n_missing, 0)
})


test_that("fit_cat works for order 1", {
  skip_on_cran()
  set.seed(42)
  y_sim <- matrix(sample(1:2, 100 * 5, replace = TRUE), nrow = 100, ncol = 5)
  
  fit0 <- fit_cat(y_sim, order = 0)
  fit1 <- fit_cat(y_sim, order = 1)
  
  expect_true(is.finite(fit1$log_l))
  expect_true(fit1$log_l >= fit0$log_l)  # Higher order fits at least as well
  expect_equal(fit1$n_params, 9)  # 1 + 4*2 = 9
})


test_that("fit_cat works for order 2", {
  skip_on_cran()
  set.seed(42)
  y_sim <- matrix(sample(1:2, 100 * 5, replace = TRUE), nrow = 100, ncol = 5)
  
  fit1 <- fit_cat(y_sim, order = 1)
  fit2 <- fit_cat(y_sim, order = 2)
  
  expect_true(is.finite(fit2$log_l))
  expect_true(fit2$log_l >= fit1$log_l)
})


test_that("fit_cat missing-data modes behave correctly", {
  skip_on_cran()
  set.seed(43)
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  y_miss <- y
  y_miss[1, 2] <- NA
  y_miss[4, 4] <- NA

  expect_error(fit_cat(y_miss, order = 1, na_action = "fail"), "Missing data")

  keep <- stats::complete.cases(y_miss)
  fit_complete_mode <- fit_cat(y_miss, order = 1, na_action = "complete")
  fit_complete_ref <- fit_cat(y_miss[keep, , drop = FALSE], order = 1)
  expect_equal(as.numeric(fit_complete_mode$log_l), as.numeric(fit_complete_ref$log_l), tolerance = 1e-10)
  expect_equal(fit_complete_mode$n_missing, 0)

  # With no missing values, marginalize should reduce to complete-data fit.
  fit_full <- fit_cat(y, order = 1)
  fit_marg_no_missing <- fit_cat(y, order = 1, na_action = "marginalize")
  expect_equal(as.numeric(fit_marg_no_missing$log_l), as.numeric(fit_full$log_l), tolerance = 1e-10)
  expect_equal(fit_marg_no_missing$settings$na_action_effective, "complete")
  expect_equal(fit_marg_no_missing$n_missing, 0)

  # EM entry point via fit_cat should be available for order 1.
  fit_em_mode <- fit_cat(y_miss, order = 1, na_action = "em", em_max_iter = 10)
  fit_em_direct <- em_cat(y_miss, order = 1, max_iter = 10)
  expect_equal(as.numeric(fit_em_mode$log_l), as.numeric(fit_em_direct$log_l), tolerance = 1e-8)
  expect_equal(fit_em_mode$settings$na_action, "em")
  expect_equal(fit_em_mode$settings$method, fit_em_direct$settings$method)
})


test_that("fit_cat marginalize handles missingness for order 2", {
  skip_on_cran()
  set.seed(44)
  y <- simulate_cat(60, 5, order = 2, n_categories = 2)
  
  # Create mixed monotone/intermittent missingness.
  y_miss <- y
  dropout_time <- sample(3:5, nrow(y_miss), replace = TRUE)
  dropout_flag <- runif(nrow(y_miss)) < 0.25
  for (i in seq_len(nrow(y_miss))) {
    if (dropout_flag[i]) {
      y_miss[i, dropout_time[i]:ncol(y_miss)] <- NA
    }
  }
  mask <- matrix(runif(length(y_miss)) < 0.08, nrow = nrow(y_miss))
  y_miss[mask] <- NA
  
  fit_miss <- fit_cat(y_miss, order = 2, na_action = "marginalize")
  ll_miss <- logL_cat(
    y_miss, order = 2, marginal = fit_miss$marginal, transition = fit_miss$transition,
    na_action = "marginalize"
  )
  
  expect_true(is.finite(fit_miss$log_l))
  expect_true(is.finite(ll_miss))

  expect_error(
    fit_cat(y_miss, order = 2, na_action = "em"),
    "currently supports order 0 and 1"
  )
})


# ============================================================
# Test 5: simulate_cat functionality
# ============================================================
test_that("simulate_cat produces valid output with defaults", {
  skip_on_cran()
  set.seed(123)
  y_uniform <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 3)
  
  expect_equal(nrow(y_uniform), 50)
  expect_equal(ncol(y_uniform), 4)
  expect_true(all(y_uniform >= 1 & y_uniform <= 3))
})


test_that("simulate_cat works with custom parameters", {
  skip_on_cran()
  marginal_custom <- list(t1 = c(0.7, 0.3))
  transition_custom <- list(
    t2 = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE),
    t4 = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
  )

  set.seed(456)
  y_custom <- simulate_cat(400, 4, order = 1, n_categories = 2,
                           marginal = marginal_custom,
                           transition = transition_custom)

  # Check that marginal of Y1 is approximately (0.7, 0.3)
  prop_y1 <- table(y_custom[, 1]) / nrow(y_custom)
  expect_equal(unname(prop_y1[1]), 0.7, tolerance = 0.08)
})


# ============================================================
# Test 6: logL_cat consistency with fit_cat
# ============================================================
test_that("logL_cat is consistent with fit_cat", {
  skip_on_cran()
  set.seed(42)
  y_sim <- matrix(sample(1:2, 100 * 5, replace = TRUE), nrow = 100, ncol = 5)
  
  fit_test <- fit_cat(y_sim, order = 1)
  ll_check <- logL_cat(y_sim, order = 1, 
                       marginal = fit_test$marginal, 
                       transition = fit_test$transition)
  
  # Use unname to ignore name attributes
  expect_equal(unname(ll_check), unname(fit_test$log_l), tolerance = 1e-8)
})


# ============================================================
# Test 7: BIC comparison
# ============================================================
test_that("bic_order_cat compares models correctly", {
  skip_on_cran()
  set.seed(42)
  y_sim <- matrix(sample(1:2, 100 * 5, replace = TRUE), nrow = 100, ncol = 5)
  
  result <- bic_order_cat(y_sim, max_order = 2)
  
  expect_equal(nrow(result$table), 3)  # Orders 0, 1, 2
  expect_true(all(c("order", "log_l", "n_params", "aic", "bic") %in% names(result$table)))
  expect_true(result$best_order %in% 0:2)
})


# ============================================================
# Test 8: Simulation-estimation consistency
# ============================================================
test_that("simulation-estimation is consistent", {
  skip_on_cran()
  set.seed(789)
  true_marginal <- list(t1 = c(0.6, 0.4))
  true_transition <- list(
    t2 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE),
    t4 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE),
    t5 = matrix(c(0.8, 0.2, 0.3, 0.7), nrow = 2, byrow = TRUE)
  )

  y_large <- simulate_cat(5000, 5, order = 1, n_categories = 2,
                          marginal = true_marginal,
                          transition = true_transition)

  fit_large <- fit_cat(y_large, order = 1)

  # Check marginal recovery - use unname to strip names
  est_marg <- unname(fit_large$marginal$t1)
  expect_equal(est_marg, true_marginal$t1, tolerance = 0.03)

  # Check transition recovery - use unname to strip dimnames
  est_trans <- unname(fit_large$transition$t2)
  expect_equal(est_trans, true_transition$t2, tolerance = 0.03)
})


# ============================================================
# Test 9: Blocks (heterogeneous)
# ============================================================
test_that("block handling works for homogeneous model", {
  skip_on_cran()
  set.seed(321)
  marginal_g1 <- list(t1 = c(0.8, 0.2))
  marginal_g2 <- list(t1 = c(0.3, 0.7))
  trans_g1 <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  )
  trans_g2 <- list(
    t2 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
  )

  y_g1 <- simulate_cat(60, 3, order = 1, n_categories = 2,
                       marginal = marginal_g1, transition = trans_g1)
  y_g2 <- simulate_cat(60, 3, order = 1, n_categories = 2,
                       marginal = marginal_g2, transition = trans_g2)

  y_combined <- rbind(y_g1, y_g2)
  blocks_vec <- c(rep(1, 60), rep(2, 60))

  # Fit homogeneous model
  fit_homo <- fit_cat(y_combined, order = 1, blocks = blocks_vec, homogeneous = TRUE)
  expect_true(is.finite(fit_homo$log_l))
})


test_that("block handling works for heterogeneous model", {
  skip_on_cran()
  set.seed(321)
  marginal_g1 <- list(t1 = c(0.8, 0.2))
  marginal_g2 <- list(t1 = c(0.3, 0.7))
  trans_g1 <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.1, 0.9), nrow = 2, byrow = TRUE)
  )
  trans_g2 <- list(
    t2 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
  )

  y_g1 <- simulate_cat(60, 3, order = 1, n_categories = 2,
                       marginal = marginal_g1, transition = trans_g1)
  y_g2 <- simulate_cat(60, 3, order = 1, n_categories = 2,
                       marginal = marginal_g2, transition = trans_g2)

  y_combined <- rbind(y_g1, y_g2)
  blocks_vec <- c(rep(1, 60), rep(2, 60))

  # Fit heterogeneous model
  fit_homo <- fit_cat(y_combined, order = 1, blocks = blocks_vec, homogeneous = TRUE)
  fit_hetero <- fit_cat(y_combined, order = 1, blocks = blocks_vec, homogeneous = FALSE)

  # Heterogeneous should fit better
  expect_true(fit_hetero$log_l > fit_homo$log_l)
  expect_equal(sort(fit_hetero$settings$block_levels), c("1", "2"))

  # Check parameter separation - block 1 has higher category-1 probability
  # than block 2 (true values 0.8 vs 0.3), direction should be correct
  est_g1_t1 <- unname(fit_hetero$marginal$block_1$t1)
  est_g2_t1 <- unname(fit_hetero$marginal$block_2$t1)
  expect_true(est_g1_t1[1] > est_g2_t1[1])
})

test_that("fit_cat preserves original block labels in settings", {
  skip_on_cran()
  set.seed(654)
  y <- simulate_cat(120, 4, order = 1, n_categories = 2)
  blocks <- rep(c(5, 2), each = 60)

  fit <- fit_cat(y, order = 1, blocks = blocks, homogeneous = FALSE)

  expect_equal(sort(fit$settings$block_levels), c("2", "5"))
  expect_equal(length(fit$marginal), 2)
})


# ============================================================
# Test 10: Order 2 model
# ============================================================
test_that("order 2 model works correctly", {
  skip_on_cran()
  set.seed(999)
  marginal_o2 <- list(
    t1 = c(0.5, 0.5),
    t2_given_1to1 = matrix(c(0.6, 0.4, 0.4, 0.6), nrow = 2, byrow = TRUE)
  )

  # For order 2, transition depends on (t-2, t-1)
  # Array dimensions: [y_{t-2}, y_{t-1}, y_t]
  trans_o2_single <- array(c(
    0.8, 0.2,  # y_{t-2}=1, y_{t-1}=1
    0.3, 0.7,  # y_{t-2}=1, y_{t-1}=2
    0.6, 0.4,  # y_{t-2}=2, y_{t-1}=1
    0.4, 0.6   # y_{t-2}=2, y_{t-1}=2
  ), dim = c(2, 2, 2))

  trans_o2 <- list(t3 = trans_o2_single, t4 = trans_o2_single, t5 = trans_o2_single)

  y_o2 <- simulate_cat(150, 5, order = 2, n_categories = 2,
                       marginal = marginal_o2, transition = trans_o2)

  fit_o2 <- fit_cat(y_o2, order = 2)
  expect_true(is.finite(fit_o2$log_l))

  # Compare with order 1 (order 2 should fit better on true AD(2) data)
  fit_o1_on_o2 <- fit_cat(y_o2, order = 1)
  expect_true(fit_o2$log_l > fit_o1_on_o2$log_l)
})
