# test-logL_cat.R - Tests for categorical AD log-likelihood

test_that("logL_cat validates inputs correctly", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(sample(1:3, 50, replace = TRUE), nrow = 10, ncol = 5)
  
  # Order must be non-negative
  expect_error(logL_cat(y, order = -1, marginal = list(), transition = list()),
               "non-negative")
  
  # Order must be less than n_time
  expect_error(logL_cat(y, order = 5, marginal = list(), transition = list()),
               "less than")
})


test_that("logL_cat missing-data modes behave correctly", {
  skip_on_cran()
  set.seed(124)
  y <- simulate_cat(50, 5, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)
  
  y_miss <- y
  y_miss[1, 3] <- NA
  y_miss[2, 5] <- NA
  
  expect_error(
    logL_cat(y_miss, order = 1, marginal = fit$marginal, transition = fit$transition,
             na_action = "fail"),
    "Missing data"
  )
  
  keep <- stats::complete.cases(y_miss)
  ll_complete_mode <- logL_cat(
    y_miss, order = 1, marginal = fit$marginal, transition = fit$transition,
    na_action = "complete"
  )
  ll_complete_ref <- logL_cat(
    y_miss[keep, , drop = FALSE], order = 1,
    marginal = fit$marginal, transition = fit$transition
  )
  expect_equal(as.numeric(ll_complete_mode), as.numeric(ll_complete_ref), tolerance = 1e-10)
  
  ll_no_missing <- logL_cat(y, order = 1, marginal = fit$marginal, transition = fit$transition)
  ll_no_missing_marg <- logL_cat(
    y, order = 1, marginal = fit$marginal, transition = fit$transition,
    na_action = "marginalize"
  )
  expect_equal(as.numeric(ll_no_missing_marg), as.numeric(ll_no_missing), tolerance = 1e-10)
})


test_that("logL_cat marginalize is finite for MAR-style missingness", {
  skip_on_cran()
  set.seed(125)
  y <- simulate_cat(60, 5, order = 1, n_categories = 2)
  
  # Monotone missingness (dropout)
  y_monotone <- y
  dropout_time <- sample(2:5, nrow(y_monotone), replace = TRUE)
  dropout_flag <- runif(nrow(y_monotone)) < 0.35
  for (i in seq_len(nrow(y_monotone))) {
    if (dropout_flag[i]) {
      y_monotone[i, dropout_time[i]:ncol(y_monotone)] <- NA
    }
  }
  
  fit_monotone <- fit_cat(y_monotone, order = 1, na_action = "marginalize")
  ll_monotone <- logL_cat(
    y_monotone, order = 1, marginal = fit_monotone$marginal,
    transition = fit_monotone$transition, na_action = "marginalize"
  )
  expect_true(is.finite(fit_monotone$log_l))
  expect_true(is.finite(ll_monotone))
  
  # Intermittent missingness
  y_intermittent <- y
  mask <- matrix(runif(length(y_intermittent)) < 0.15, nrow = nrow(y_intermittent))
  y_intermittent[mask] <- NA
  
  fit_intermittent <- fit_cat(y_intermittent, order = 1, na_action = "marginalize")
  ll_intermittent <- logL_cat(
    y_intermittent, order = 1, marginal = fit_intermittent$marginal,
    transition = fit_intermittent$transition, na_action = "marginalize"
  )
  expect_true(is.finite(fit_intermittent$log_l))
  expect_true(is.finite(ll_intermittent))
})


test_that("logL_cat computes correct log-likelihood for order 0", {
  skip_on_cran()
  set.seed(234)
  n_subjects <- 50
  n_time <- 4
  c <- 2
  
  # Create uniform marginals
  marginal <- lapply(1:n_time, function(k) c(0.5, 0.5))
  names(marginal) <- 1:n_time
  
  # Simulate uniform data
  y <- matrix(sample(1:c, n_subjects * n_time, replace = TRUE), 
              nrow = n_subjects, ncol = n_time)
  
  # Compute log-likelihood
  ll <- logL_cat(y, order = 0, marginal = marginal, n_categories = c)
  
  # Manual computation: each observation has prob 0.5
  ll_manual <- n_subjects * n_time * log(0.5)
  
  expect_equal(as.numeric(ll), ll_manual)
})


test_that("logL_cat computes correct log-likelihood for order 1", {
  skip_on_cran()
  set.seed(345)
  n_subjects <- 100
  n_time <- 3
  c <- 2
  
  # Known marginal
  marginal <- list(t1 = c(0.6, 0.4))
  
  # Known transitions
  transition <- list(
    t2 = matrix(c(0.8, 0.2,
                  0.3, 0.7), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.8, 0.2,
                  0.3, 0.7), nrow = 2, byrow = TRUE)
  )
  
  # Simulate data from this model
  y <- simulate_cat(n_subjects, n_time, order = 1, n_categories = c,
                    marginal = marginal, transition = transition)
  
  # Compute log-likelihood at true parameters
  ll_true <- logL_cat(y, order = 1, marginal = marginal, 
                      transition = transition, n_categories = c)
  
  # Compute at wrong parameters (uniform)
  marginal_wrong <- list(t1 = c(0.5, 0.5))
  transition_wrong <- list(
    t2 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2, byrow = TRUE)
  )
  ll_wrong <- logL_cat(y, order = 1, marginal = marginal_wrong, 
                       transition = transition_wrong, n_categories = c)
  
  # Log-likelihood should be higher at true parameters
  expect_true(as.numeric(ll_true) > as.numeric(ll_wrong))
})


test_that("logL_cat matches fit_cat log-likelihood", {
  skip_on_cran()
  set.seed(456)
  n_subjects <- 80
  n_time <- 5
  
  # Simulate data
  y <- simulate_cat(n_subjects, n_time, order = 1, n_categories = 2)
  
  # Fit model
  fit <- fit_cat(y, order = 1)
  
  # Compute log-likelihood at fitted parameters
  ll_check <- logL_cat(y, order = 1, 
                       marginal = fit$marginal, 
                       transition = fit$transition,
                       n_categories = 2)
  
  expect_equal(as.numeric(ll_check), as.numeric(fit$log_l), tolerance = 1e-10)
})


test_that("logL_cat works for order 2", {
  skip_on_cran()
  set.seed(567)
  n_subjects <- 100
  n_time <- 5
  c <- 2
  
  # Simulate order 2 data
  y <- simulate_cat(n_subjects, n_time, order = 2, n_categories = c)
  
  # Fit order 2 model
  fit <- fit_cat(y, order = 2)
  
  # Check that log-likelihood matches
  ll_check <- logL_cat(y, order = 2,
                       marginal = fit$marginal,
                       transition = fit$transition,
                       n_categories = c)
  
  expect_equal(as.numeric(ll_check), as.numeric(fit$log_l), tolerance = 1e-10)
})


test_that("logL_cat handles more than 2 categories", {
  skip_on_cran()
  set.seed(678)
  n_subjects <- 80
  n_time <- 4
  c <- 3
  
  # Simulate with 3 categories
  y <- simulate_cat(n_subjects, n_time, order = 1, n_categories = c)
  
  # Fit and verify
  fit <- fit_cat(y, order = 1, n_categories = c)
  
  ll_check <- logL_cat(y, order = 1,
                       marginal = fit$marginal,
                       transition = fit$transition,
                       n_categories = c)
  
  expect_equal(as.numeric(ll_check), as.numeric(fit$log_l), tolerance = 1e-10)
})


test_that("logL_cat works with heterogeneous blocks", {
  skip_on_cran()
  set.seed(789)
  n1 <- 40
  n2 <- 40
  n_time <- 4
  c <- 2
  
  # Simulate two groups with different parameters
  marg1 <- list(t1 = c(0.7, 0.3))
  trans1 <- list(
    t2 = matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE),
    t4 = matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
  )
  
  marg2 <- list(t1 = c(0.4, 0.6))
  trans2 <- list(
    t2 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
    t4 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE)
  )
  
  y1 <- simulate_cat(n1, n_time, order = 1, n_categories = c,
                     marginal = marg1, transition = trans1)
  y2 <- simulate_cat(n2, n_time, order = 1, n_categories = c,
                     marginal = marg2, transition = trans2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, n1), rep(2, n2))
  
  # Fit heterogeneous model
  fit <- fit_cat(y, order = 1, blocks = blocks, homogeneous = FALSE)
  
  # Check log-likelihood
  ll_check <- logL_cat(y, order = 1,
                       marginal = fit$marginal,
                       transition = fit$transition,
                       blocks = blocks,
                       homogeneous = FALSE,
                       n_categories = c)
  
  expect_equal(as.numeric(ll_check), as.numeric(fit$log_l), tolerance = 1e-10)
})


test_that("logL_cat is additive across independent subjects", {
  skip_on_cran()
  set.seed(901)
  n_time <- 4
  c <- 2
  
  # Create two subjects
  y <- matrix(sample(1:c, 2 * n_time, replace = TRUE), nrow = 2, ncol = n_time)
  
  # Fit model to get valid parameters
  fit <- fit_cat(y, order = 1, n_categories = c)
  
  # Full log-likelihood
  ll_full <- logL_cat(y, order = 1, marginal = fit$marginal,
                      transition = fit$transition, n_categories = c)
  
  # Individual log-likelihoods
  ll_1 <- logL_cat(y[1, , drop = FALSE], order = 1, marginal = fit$marginal,
                   transition = fit$transition, n_categories = c)
  ll_2 <- logL_cat(y[2, , drop = FALSE], order = 1, marginal = fit$marginal,
                   transition = fit$transition, n_categories = c)
  
  expect_equal(as.numeric(ll_full), as.numeric(ll_1) + as.numeric(ll_2), tolerance = 1e-10)
})


test_that("logL_cat handles n_categories override", {
  skip_on_cran()
  set.seed(112)
  n_subjects <- 50
  n_time <- 4
  
  # Generate data with only categories 1 and 2 observed
  y <- matrix(sample(1:2, n_subjects * n_time, replace = TRUE), 
              nrow = n_subjects, ncol = n_time)
  
  # Fit with 3 categories (even though category 3 not observed)
  fit <- fit_cat(y, order = 1, n_categories = 3)
  
  # Verify log-likelihood works with the override
  ll <- logL_cat(y, order = 1, marginal = fit$marginal,
                 transition = fit$transition, n_categories = 3)
  
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(as.numeric(ll), as.numeric(fit$log_l), tolerance = 1e-10)
})


test_that("logL_cat likelihood increases with correct model", {
  skip_on_cran()
  set.seed(223)
  n_subjects <- 60
  n_time <- 5
  c <- 2
  
  # Strong dependence structure
  marginal_true <- list(t1 = c(0.5, 0.5))
  transition_true <- list(
    t2 = matrix(c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE),
    t3 = matrix(c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE),
    t4 = matrix(c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE),
    t5 = matrix(c(0.95, 0.05, 0.05, 0.95), 2, byrow = TRUE)
  )
  
  # Simulate from true model
  y <- simulate_cat(n_subjects, n_time, order = 1, n_categories = c,
                    marginal = marginal_true, transition = transition_true)
  
  # Log-likelihood at true parameters
  ll_true <- logL_cat(y, order = 1, marginal = marginal_true,
                      transition = transition_true, n_categories = c)
  
  # Log-likelihood at MLE (should be >= true)
  fit <- fit_cat(y, order = 1, n_categories = c)
  ll_mle <- fit$log_l
  
  # MLE should be at least as good as true parameters
  expect_true(as.numeric(ll_mle) >= as.numeric(ll_true) - 1e-6)  # Small tolerance for numerical issues
})
