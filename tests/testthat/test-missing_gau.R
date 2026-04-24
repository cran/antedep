#' Test suite for missing data support
#'
#' Tests for missing_utils.R and AD missing data functionality

library(testthat)

# ==== Test missing_utils.R ====

test_that(".validate_missing works for complete data", {
  skip_on_cran()
  y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  
  result <- .validate_missing(y)
  
  expect_false(result$has_missing)
  expect_equal(result$n_complete, 10)
  expect_equal(result$n_intermittent, 0)
  expect_equal(result$pct_missing, 0)
  expect_true(all(result$patterns == "complete"))
})

test_that(".validate_missing identifies monotone dropout", {
  skip_on_cran()
  y <- matrix(c(
    1, 2, 3, 4, 5,
    1, 2, 3, NA, NA,
    1, 2, NA, NA, NA
  ), nrow = 3, byrow = TRUE)
  
  result <- .validate_missing(y)
  
  expect_true(result$has_missing)
  expect_equal(result$n_complete, 1)
  expect_equal(result$patterns[1], "complete")
  expect_equal(result$patterns[2], "dropout")
  expect_equal(result$patterns[3], "dropout")
})

test_that(".validate_missing identifies monotone drop-in", {
  skip_on_cran()
  y <- matrix(c(
    1, 2, 3, 4, 5,
    NA, NA, 3, 4, 5,
    NA, 2, 3, 4, 5
  ), nrow = 3, byrow = TRUE)
  
  result <- .validate_missing(y)
  
  expect_equal(result$patterns[2], "dropin")
  expect_equal(result$patterns[3], "dropin")
})

test_that(".validate_missing identifies intermittent missing", {
  skip_on_cran()
  y <- matrix(c(
    1, NA, 3, NA, 5,
    1, 2, NA, 4, 5
  ), nrow = 2, byrow = TRUE)
  
  result <- .validate_missing(y)
  
  expect_equal(result$n_intermittent, 2)
  expect_true(all(result$patterns == "intermittent"))
})

test_that(".get_truncation_bound works correctly", {
  skip_on_cran()
  y <- matrix(c(
    5, 10, 15, 20,
    3, 8, NA, 18,
    4, 9, 14, 19
  ), nrow = 3, byrow = TRUE)
  
  # Subject 2, time 3 is missing
  # Subject 2 max = 18, other subjects at time 3: 14
  # Expected: 1.2 * max(18, 14) = 1.2 * 18 = 21.6 -> 21
  bound <- .get_truncation_bound(y, subject_idx = 2, time_idx = 3, buffer = 1.2)
  
  expect_equal(bound, 21)
})

test_that(".get_truncation_bound handles edge cases", {
  skip_on_cran()
  # All missing in column
  y <- matrix(c(
    1, NA, 3,
    2, NA, 4,
    3, NA, 5
  ), nrow = 3, byrow = TRUE)
  
  bound <- .get_truncation_bound(y, subject_idx = 1, time_idx = 2, min_bound = 10)
  expect_gte(bound, 10)  # Should return at least min_bound
})

test_that(".safe_log handles zero correctly", {
  skip_on_cran()
  x <- c(-1, 0, 1, 2, exp(1))
  
  result <- .safe_log(x)
  
  expect_equal(result[1], -Inf)
  expect_equal(result[2], 0)  # 0 * log(0) = 0 convention
  expect_equal(result[3], log(1))
  expect_equal(result[4], log(2))
  expect_equal(result[5], 1)
})

test_that(".extract_complete_cases works", {
  skip_on_cran()
  y <- matrix(c(
    1, 2, 3,
    1, NA, 3,
    1, 2, 3,
    NA, NA, NA
  ), nrow = 4, byrow = TRUE)
  
  blocks <- c(1, 1, 2, 2)
  
  result <- .extract_complete_cases(y, blocks, warn = FALSE)
  
  expect_equal(nrow(result$y), 2)
  expect_equal(result$blocks, c(1, 2))
  expect_equal(result$n_removed, 2)
})

test_that(".extract_complete_cases warns when many removed", {
  skip_on_cran()
  y <- matrix(rnorm(100), nrow = 20, ncol = 5)
  y[1:15, 1] <- NA  # Make 75% incomplete
  
  expect_warning(
    result <- .extract_complete_cases(y, NULL),
    "Only 5/20"
  )
})

test_that(".extract_complete_cases errors when no complete cases", {
  skip_on_cran()
  y <- matrix(NA, nrow = 5, ncol = 3)
  
  expect_error(
    .extract_complete_cases(y, NULL),
    "No complete cases"
  )
})

# ==== Test logL_gau with missing data ====

test_that("logL_gau with na_action='fail' errors on missing data", {
  skip_on_cran()
  y <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2, byrow = TRUE)
  
  expect_error(
    logL_gau(y, order = 1, mu = rep(0, 3), phi = c(0.5, 0.5), 
            sigma = rep(1, 3), na_action = "fail"),
    "contains NA"
  )
})

test_that("logL_gau with na_action='complete' uses only complete cases", {
  skip_on_cran()
  set.seed(123)
  y_complete <- matrix(rnorm(30), nrow = 10, ncol = 3)
  y_missing <- y_complete
  y_missing[1:5, 1] <- NA  # Make first 5 subjects incomplete
  
  mu <- c(0, 0, 0)
  phi <- c(0.5, 0.5)
  sigma <- c(1, 1, 1)
  
  # Should only use last 5 subjects
  ll_complete <- logL_gau(y_missing, order = 1, mu = mu, phi = phi, sigma = sigma,
                         na_action = "complete")
  ll_last5 <- logL_gau(y_complete[6:10, ], order = 1, mu = mu, phi = phi, sigma = sigma,
                      na_action = "fail")
  
  expect_equal(ll_complete, ll_last5, tolerance = 1e-10)
})

test_that("logL_gau with na_action='marginalize' works for monotone dropout", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(15), nrow = 5, ncol = 3)
  y[1, 3] <- NA  # One missing at end
  
  mu <- c(0, 0, 0)
  phi <- c(0.5, 0.5)
  sigma <- c(1, 1, 1)
  
  # Should not error
  ll <- logL_gau(y, order = 1, mu = mu, phi = phi, sigma = sigma,
                na_action = "marginalize")
  
  expect_true(is.finite(ll))
})

test_that(".build_gau_covariance works for order 0", {
  skip_on_cran()
  sigma <- c(1, 2, 3)
  Sigma <- .build_gau_covariance(order = 0, phi = numeric(0), sigma = sigma, n_time = 3)
  
  expect_equal(diag(Sigma), sigma^2)
  expect_equal(Sigma[1, 2], 0)
  expect_equal(Sigma[1, 3], 0)
})

test_that(".build_gau_covariance works for order 1", {
  skip_on_cran()
  sigma <- c(1, 1, 1)
  phi <- c(0.5, 0.5)
  
  Sigma <- .build_gau_covariance(order = 1, phi = phi, sigma = sigma, n_time = 3)
  
  # Check variances increase (unless phi=0)
  expect_gte(Sigma[2, 2], Sigma[1, 1])
  expect_gte(Sigma[3, 3], Sigma[2, 2])
  
  # Check symmetry
  expect_equal(Sigma[1, 2], Sigma[2, 1])
})

test_that(".build_gau_covariance works for order 2", {
  skip_on_cran()
  set.seed(321)

  n_time <- 6
  sigma <- c(1.0, 1.1, 1.2, 1.0, 0.9, 1.1)
  phi <- matrix(0, nrow = 2, ncol = n_time)
  phi[1, 2:n_time] <- 0.45
  phi[2, 3:n_time] <- -0.20

  Sigma <- .build_gau_covariance(order = 2, phi = phi, sigma = sigma, n_time = n_time)

  expect_equal(dim(Sigma), c(n_time, n_time))
  expect_true(isSymmetric(Sigma, tol = 1e-10))
  expect_true(all(diag(Sigma) > 0))

  y <- simulate_gau(
    n_subjects = 5000,
    n_time = n_time,
    order = 2,
    mu = 0,
    phi = phi,
    sigma = sigma
  )
  Sigma_emp <- stats::cov(y)

  expect_equal(diag(Sigma_emp), diag(Sigma), tolerance = 0.12)
  expect_equal(Sigma_emp[3, 1], Sigma[3, 1], tolerance = 0.12)
  expect_equal(Sigma_emp[6, 4], Sigma[6, 4], tolerance = 0.12)
})

# ==== Test fit_gau with missing data ====

test_that("fit_gau with na_action='fail' errors on missing data", {
  skip_on_cran()
  y <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2, byrow = TRUE)
  
  expect_error(
    fit_gau(y, order = 1, na_action = "fail"),
    "contains NA"
  )
})

test_that("fit_gau with na_action='complete' removes incomplete subjects", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  y[1:3, 1] <- NA
  
  expect_warning(
    fit <- fit_gau(y, order = 1, na_action = "complete"),
    "Only 7/10"
  )
  
  expect_equal(fit$settings$n_subjects, 7)
  expect_equal(fit$n_obs, 7 * 5)
  expect_equal(fit$n_missing, 0)
})

test_that("fit_gau with complete data gives same result regardless of na_action", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  
  fit1 <- fit_gau(y, order = 1, na_action = "fail")
  fit2 <- fit_gau(y, order = 1, na_action = "complete")
  
  expect_equal(fit1$log_l, fit2$log_l, tolerance = 1e-6)
  expect_equal(fit1$mu, fit2$mu, tolerance = 1e-6)
})

# ==== Test EM algorithm ====

test_that(".initialize_gau_em works with enough complete cases", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(60), nrow = 20, ncol = 3)
  y[1:5, 1] <- NA  # 75% complete
  
  init <- .initialize_gau_em(y, order = 1, blocks = NULL, estimate_mu = TRUE)
  
  expect_length(init$mu, 3)
  expect_length(init$sigma, 3)
  expect_true(length(init$phi) %in% c(2, 3))
  if (length(init$phi) == 3) {
    expect_equal(init$phi[1], 0, tolerance = 1e-10)
  }
  expect_true(all(is.finite(init$mu)))
  expect_true(all(init$sigma > 0))
})

test_that(".initialize_gau_em falls back to marginal when few complete cases", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(30), nrow = 10, ncol = 3)
  y[1:9, 1] <- NA  # Only 1 complete case
  
  init <- .initialize_gau_em(y, order = 1, blocks = NULL, estimate_mu = TRUE)
  
  # Should still initialize
  expect_length(init$mu, 3)
  expect_true(all(is.finite(init$mu)))
})

test_that(".em_e_step_gau computes sufficient statistics", {
  skip_on_cran()
  set.seed(123)
  y <- matrix(rnorm(15), nrow = 5, ncol = 3)
  y[1, 3] <- NA
  
  mu <- c(0, 0, 0)
  phi <- c(0.5, 0.5)
  sigma <- c(1, 1, 1)
  
  suff_stats <- .em_e_step_gau(y, order = 1, mu, phi, sigma, blocks = NULL, tau = 0)
  
  expect_length(suff_stats$S1, 3)
  expect_equal(dim(suff_stats$S2), c(3, 3))
  expect_equal(suff_stats$n, 5)
  expect_true(all(is.finite(suff_stats$S1)))
})

test_that(".em_m_step_gau updates parameters", {
  skip_on_cran()
  # Simple sufficient statistics
  S1 <- c(5, 10, 15)
  S2 <- matrix(c(
    5, 7, 9,
    7, 20, 25,
    9, 25, 50
  ), nrow = 3, byrow = TRUE)
  suff_stats <- list(S1 = S1, S2 = S2, n = 10)
  
  params <- .em_m_step_gau(suff_stats, order = 1, blocks = NULL,
                          estimate_mu = TRUE, n_subjects = 10, n_time = 3)
  
  expect_length(params$mu, 3)
  expect_length(params$phi, 2)
  expect_length(params$sigma, 3)
  expect_true(all(params$sigma > 0))
})

test_that("fit_gau with EM converges on monotone missing data", {
  skip_on_cran()
  
  set.seed(123)
  y <- matrix(rnorm(30), nrow = 10, ncol = 3)
  y[1:3, 3] <- NA  # Monotone dropout
  
  fit <- fit_gau(y, order = 1, na_action = "em", em_max_iter = 50,
                em_verbose = FALSE)
  
  expect_true(fit$em_converged)
  expect_lte(fit$em_iterations, 50)
  expect_true(is.finite(fit$log_l))
  expect_gt(fit$n_missing, 0)
  expect_gt(fit$pct_missing, 0)
})

test_that("fit_gau EM handles blocks", {
  skip_on_cran()
  
  set.seed(123)
  y <- matrix(rnorm(40), nrow = 20, ncol = 2)
  y[1:5, 2] <- NA
  blocks <- c(rep(1, 10), rep(2, 10))
  
  fit <- fit_gau(y, order = 1, blocks = blocks, na_action = "em",
                em_max_iter = 30, em_verbose = FALSE)
  
  expect_true(fit$em_converged)
  expect_length(fit$tau, 2)
  expect_equal(fit$tau[1], 0)
})

test_that("EM log-likelihood is monotone increasing", {
  skip_on_cran()
  
  set.seed(123)
  y <- matrix(rnorm(25), nrow = 5, ncol = 5)
  y[sample(length(y), 5)] <- NA  # Random missing
  
  fit <- fit_gau(y, order = 1, na_action = "em", em_max_iter = 20,
                em_verbose = FALSE)
  
  ll_trace <- fit$em_ll_trace
  
  # Check monotonicity (allowing tiny numerical errors)
  diffs <- diff(ll_trace)
  expect_true(all(diffs >= -1e-6))
})

test_that("fit_gau EM order 2 converges and recovers dependence parameters", {
  skip_on_cran()

  set.seed(2026)

  n_subjects <- 400
  n_time <- 6
  mu_true <- seq(0.2, 0.7, length.out = n_time)
  sigma_true <- c(1.0, 1.1, 1.0, 0.9, 1.1, 1.0)
  phi_true <- matrix(0, nrow = 2, ncol = n_time)
  phi_true[1, 2:n_time] <- 0.50
  phi_true[2, 3:n_time] <- -0.25

  y_complete <- simulate_gau(
    n_subjects = n_subjects,
    n_time = n_time,
    order = 2,
    mu = mu_true,
    phi = phi_true,
    sigma = sigma_true
  )

  y_missing <- y_complete
  missing_idx <- sample(length(y_missing), size = round(0.12 * length(y_missing)))
  y_missing[missing_idx] <- NA

  fit_em <- NULL
  expect_warning(
    fit_em <- fit_gau(
      y_missing,
      order = 2,
      na_action = "em",
      em_max_iter = 120,
      em_tol = 1e-6,
      em_verbose = FALSE
    ),
    "provisional implementation"
  )
  fit_complete <- fit_gau(y_complete, order = 2, na_action = "fail")

  expect_true(fit_em$em_converged)
  expect_true(is.matrix(fit_em$phi))
  expect_equal(dim(fit_em$phi), c(2, n_time))
  expect_equal(fit_em$phi[1, 1], 0, tolerance = 1e-10)
  expect_equal(fit_em$phi[2, 1], 0, tolerance = 1e-10)
  expect_equal(fit_em$phi[2, 2], 0, tolerance = 1e-10)
  expect_true(all(diff(fit_em$em_ll_trace) >= -1e-6))

  # Agreement with truth and complete-data fit.
  expect_equal(fit_em$phi[1, 3:n_time], phi_true[1, 3:n_time], tolerance = 0.25)
  expect_equal(fit_em$phi[2, 3:n_time], phi_true[2, 3:n_time], tolerance = 0.25)
  expect_equal(fit_em$phi[1, 3:n_time], fit_complete$phi[1, 3:n_time], tolerance = 0.20)
  expect_equal(fit_em$phi[2, 3:n_time], fit_complete$phi[2, 3:n_time], tolerance = 0.20)
  expect_equal(fit_em$sigma[3:n_time], fit_complete$sigma[3:n_time], tolerance = 0.20)
})

# ==== Integration test ====

test_that("Complete workflow: simulate -> introduce missing -> fit with EM", {
  skip_on_cran()
  
  set.seed(42)
  
  # Simulate complete data (would use simulate_gau in practice)
  n <- 30
  p <- 4
  y_complete <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Introduce 10% MCAR missingness
  n_missing <- round(0.1 * n * p)
  missing_idx <- sample(n * p, n_missing)
  y_missing <- y_complete
  y_missing[missing_idx] <- NA
  
  # Fit with EM
  fit_em <- fit_gau(y_missing, order = 1, na_action = "em", 
                   em_max_iter = 50, em_verbose = FALSE)
  
  # Fit complete-case for comparison
  fit_cc <- suppressWarnings(
    fit_gau(y_missing, order = 1, na_action = "complete")
  )
  
  # EM should converge
  expect_true(fit_em$em_converged)
  
  # EM should have better (or equal) log-likelihood
  # (comparing on different data, so just check both are finite)
  expect_true(is.finite(fit_em$log_l))
  expect_true(is.finite(fit_cc$log_l))
  
  # EM uses more data
  expect_gt(fit_em$n_obs, fit_cc$n_obs)
})
