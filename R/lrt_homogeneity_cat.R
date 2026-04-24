# test_homogeneity_cat.R - Likelihood ratio test for homogeneity in categorical AD

#' Likelihood ratio test for homogeneity across groups (categorical AD data)
#'
#' Tests whether multiple groups share the same transition probability parameters
#' in a categorical antedependence model.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c. Can be NULL if both fit_null and
#'   fit_alt are provided.
#' @param blocks Integer vector of length n_subjects specifying group membership.
#'   Required unless pre-fitted models are provided.
#' @param order Antedependence order p. Default is 1.
#' @param n_categories Number of categories. If NULL, inferred from data.
#' @param fit_null Optional pre-fitted homogeneous model (class "cat_fit" with
#'   homogeneous = TRUE). If provided, y is not required for fitting under H0.
#' @param fit_alt Optional pre-fitted heterogeneous model (class "cat_fit" with
#'   homogeneous = FALSE). If provided, y is not required for fitting under H1.
#' @param test Type of test statistic. One of \code{"lrt"} (default),
#'   \code{"score"}, or \code{"mlrt"}.
#'
#' @return A list of class \code{"cat_lrt"} containing:
#' \describe{
#'   \item{method}{Inference method used: one of \code{"lrt"}, \code{"score"},
#'     \code{"mlrt"}, or \code{"wald"}.}
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{fit_null}{Fitted homogeneous model (H0)}
#'   \item{fit_alt}{Fitted heterogeneous model (H1)}
#'   \item{n_groups}{Number of groups}
#'   \item{table}{Summary data frame}
#' }
#'
#' @details
#' The null hypothesis is that all G groups share the same transition probability
#' parameters:
#' \deqn{H_0: \pi^{(1)} = \pi^{(2)} = \ldots = \pi^{(G)}}
#'
#' The alternative hypothesis allows each group to have its own parameters.
#'
#' The degrees of freedom are:
#' \deqn{df = (G-1) \times k}
#' where G is the number of groups and k is the number of free parameters per
#' population.
#'
#' @examples
#' \donttest{
#' # Simulate data with different transition probabilities for two groups
#' set.seed(123)
#' marg1 <- list(t1 = c(0.7, 0.3))
#' marg2 <- list(t1 = c(0.4, 0.6))
#' trans1 <- list(t2 = matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE),
#'                t3 = matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE))
#' trans2 <- list(t2 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
#'                t3 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE))
#'
#' y1 <- simulate_cat(100, 3, order = 1, n_categories = 2,
#'                    marginal = marg1, transition = trans1)
#' y2 <- simulate_cat(100, 3, order = 1, n_categories = 2,
#'                    marginal = marg2, transition = trans2)
#' y <- rbind(y1, y2)
#' blocks <- c(rep(1, 100), rep(2, 100))
#'
#' # Test homogeneity
#' test <- test_homogeneity_cat(y, blocks, order = 1)
#' print(test)
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @seealso \code{\link{fit_cat}}, \code{\link{test_order_cat}}
#'
#' @export
test_homogeneity_cat <- function(y = NULL, blocks = NULL, order = 1,
                                 n_categories = NULL,
                                 fit_null = NULL, fit_alt = NULL,
                                 test = c("lrt", "score", "mlrt")) {
  test <- match.arg(test)
  
  # Validate that we have either y+blocks or pre-fitted models
  if (is.null(y) && (is.null(fit_null) || is.null(fit_alt))) {
    stop("Either y and blocks must be provided, or both fit_null and fit_alt must be provided")
  }
  use_missing <- !is.null(y) && anyNA(y)
  na_action_fit <- if (use_missing) "marginalize" else "fail"
  
  if (!is.null(y) && is.null(blocks)) {
    stop("blocks must be provided when y is provided")
  }
  
  # Fit models if not provided
  if (is.null(fit_null)) {
    fit_null <- fit_cat(y, order = order, blocks = blocks,
                        homogeneous = TRUE, n_categories = n_categories,
                        na_action = na_action_fit)
  } else {
    if (!inherits(fit_null, "cat_fit")) {
      stop("fit_null must be a cat_fit object")
    }
    if (!fit_null$settings$homogeneous) {
      warning("fit_null should be a homogeneous model (homogeneous = TRUE)")
    }
  }
  
  if (is.null(fit_alt)) {
    fit_alt <- fit_cat(y, order = order, blocks = blocks,
                       homogeneous = FALSE, n_categories = n_categories,
                       na_action = na_action_fit)
  } else {
    if (!inherits(fit_alt, "cat_fit")) {
      stop("fit_alt must be a cat_fit object")
    }
    if (fit_alt$settings$homogeneous) {
      warning("fit_alt should be a heterogeneous model (homogeneous = FALSE)")
    }
  }
  
  # Verify they have same order
  if (fit_null$settings$order != fit_alt$settings$order) {
    stop("fit_null and fit_alt must have the same order")
  }
  
  # Extract info
  log_l_null <- fit_null$log_l
  log_l_alt <- fit_alt$log_l
  n_params_null <- fit_null$n_params
  n_params_alt <- fit_alt$n_params
  n_groups <- fit_alt$settings$n_blocks
  
  # Compute LRT statistic
  lrt_stat_raw <- -2 * (log_l_null - log_l_alt)
  
  # Degrees of freedom = (G-1) * k, where k is params per population
  df <- n_params_alt - n_params_null
  
  # Select test statistic
  stat_value <- lrt_stat_raw
  e_hat_mlrt <- NA_real_
  if (identical(test, "score")) {
    if (is.null(y)) {
      stop("y and blocks must be provided when test = 'score'")
    }
    if (anyNA(y) || .cat_fit_uses_missing_likelihood(fit_null) || .cat_fit_uses_missing_likelihood(fit_alt)) {
      .stop_cat_missing_inference("test_homogeneity_cat(test = 'score')")
    }
    stat_value <- .cat_score_homogeneity(
      y = y,
      blocks = blocks,
      p = fit_null$settings$order,
      c = fit_null$settings$n_categories,
      fit_null = fit_null
    )
  } else if (identical(test, "mlrt")) {
    p <- fit_null$settings$order
    c <- fit_null$settings$n_categories
    n_subjects <- fit_null$settings$n_subjects
    n_time <- fit_null$settings$n_time
    blocks_eff <- .cat_resolve_blocks(blocks, fit_alt)

    if (length(blocks_eff) != n_subjects) {
      stop("blocks length does not match n_subjects for modified homogeneity test")
    }

    simulate_fn <- function() {
      simulate_cat(
        n_subjects = n_subjects,
        n_time = n_time,
        order = p,
        n_categories = c,
        marginal = fit_null$marginal,
        transition = fit_null$transition,
        blocks = blocks_eff,
        homogeneous = TRUE
      )
    }
    lrt_raw_fn <- function(y_b) {
      fit0_b <- fit_cat(
        y_b,
        order = p,
        blocks = blocks_eff,
        homogeneous = TRUE,
        n_categories = c,
        na_action = "fail"
      )
      fit1_b <- fit_cat(
        y_b,
        order = p,
        blocks = blocks_eff,
        homogeneous = FALSE,
        n_categories = c,
        na_action = "fail"
      )
      -2 * (fit0_b$log_l - fit1_b$log_l)
    }

    e_hat_mlrt <- .cat_mlrt_expected_lrt_boot(simulate_fn, lrt_raw_fn)
    if (is.finite(e_hat_mlrt) && e_hat_mlrt > 0) {
      stat_value <- lrt_stat_raw * df / e_hat_mlrt
    } else {
      warning(
        "Modified LRT scaling factor could not be estimated; returning unmodified LRT statistic.",
        call. = FALSE
      )
    }
  }

  # P-value from chi-square
  p_value <- stats::pchisq(stat_value, df = df, lower.tail = FALSE)
  
  # Build summary table
  table_df <- data.frame(
    model = c("Homogeneous (H0)", "Heterogeneous (H1)"),
    log_l = c(log_l_null, log_l_alt),
    n_params = c(n_params_null, n_params_alt),
    aic = c(fit_null$aic, fit_alt$aic),
    bic = c(fit_null$bic, fit_alt$bic)
  )
  
  # Assemble output
  out <- list(
    method = test,
    lrt_stat = stat_value,
    statistic = stat_value,
    lrt_stat_raw = lrt_stat_raw,
    e_hat_mlrt = e_hat_mlrt,
    df = df,
    p_value = p_value,
    test = test,
    fit_null = fit_null,
    fit_alt = fit_alt,
    n_groups = n_groups,
    order = fit_null$settings$order,
    table = table_df
  )
  
  class(out) <- "cat_lrt"
  out
}


#' Likelihood ratio test for time-invariance (categorical data)
#'
#' Tests whether transition probabilities are constant over time in a
#' categorical antedependence model.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c.
#' @param order Antedependence order p. Default is 1.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership.
#' @param homogeneous Logical. If TRUE (default), parameters are shared across
#'   all groups.
#' @param n_categories Number of categories. If NULL, inferred from data.
#' @param test Type of test statistic. One of \code{"lrt"} (default),
#'   \code{"score"}, or \code{"mlrt"}.
#'
#' @return A list of class \code{"cat_lrt"} containing:
#' \describe{
#'   \item{method}{Inference method used: one of \code{"lrt"}, \code{"score"},
#'     \code{"mlrt"}, or \code{"wald"}.}
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{fit_null}{Fitted time-invariant model (H0)}
#'   \item{fit_alt}{Fitted time-varying model (H1)}
#'   \item{table}{Summary data frame}
#' }
#'
#' @details
#' The null hypothesis is that all transition probabilities (for k > p) are
#' equal across time:
#' \deqn{H_0: \pi_{y_k | y_{k-p}, \ldots, y_{k-1}} \text{ is constant for } k = p+1, \ldots, n}
#'
#' This reduces (n-p) separate transition matrices/arrays to a single one.
#'
#' The degrees of freedom are:
#' \deqn{df = (c-1) \times c^p \times (n - p - 1)}
#'
#' This function currently supports complete data only. If \code{y} contains
#' missing values, use model-fitting functions (for example \code{fit_cat})
#' directly with missing-data handling instead of this test wrapper.
#'
#' @examples
#' \donttest{
#' # Simulate data with time-invariant transitions
#' set.seed(123)
#' y <- simulate_cat(200, 6, order = 1, n_categories = 2)
#'
#' # Test time-invariance
#' test <- test_timeinvariance_cat(y, order = 1)
#' print(test)
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @seealso \code{\link{fit_cat}}, \code{\link{test_order_cat}}
#'
#' @export
test_timeinvariance_cat <- function(y, order = 1, blocks = NULL,
                                    homogeneous = TRUE, n_categories = NULL,
                                    test = c("lrt", "score", "mlrt")) {
  test <- match.arg(test)
  if (anyNA(y)) {
    .stop_cat_missing_inference("test_timeinvariance_cat")
  }
  
  # Validate data
  validated <- .validate_y_cat(y, n_categories)
  y <- validated$y
  c <- validated$n_categories
  n_time <- ncol(y)
  n_subjects <- nrow(y)
  p <- as.integer(order)
  
  if (p < 1) {
    stop("Time-invariance test requires order >= 1")
  }
  if (p >= n_time) {
    stop("order must be less than n_time")
  }
  
  # Fit unconstrained (time-varying) model - this is the alternative
  fit_alt <- fit_cat(y, order = p, blocks = blocks,
                     homogeneous = homogeneous, n_categories = c)
  
  # Fit constrained (time-invariant) model
  # This requires pooling all time points p+1 to n
  fit_null <- .fit_cat_timeinvariant(y, p, c, blocks, homogeneous)
  
  # Compute LRT statistic
  lrt_stat_raw <- -2 * (fit_null$log_l - fit_alt$log_l)
  
  # Degrees of freedom
  # Time-varying has (n-p) separate transition matrices, each with (c-1)*c^p params
  # Time-invariant has 1 transition matrix with (c-1)*c^p params
  # df = (n-p-1) * (c-1) * c^p
  n_trans_reduced <- (n_time - p - 1) * (c - 1) * c^p
  
  # Also need to account for marginal if p >= 1
  # For now, we only test time-invariance of transitions (not marginal)
  df <- n_trans_reduced
  
  # Adjust for heterogeneous case
  if (!homogeneous && !is.null(blocks)) {
    n_blocks <- length(unique(blocks))
    df <- df * n_blocks
  }
  
  # Select test statistic
  stat_value <- lrt_stat_raw
  e_hat_mlrt <- NA_real_
  if (identical(test, "score")) {
    pops <- .cat_split_populations(y, blocks = blocks, fit = fit_null)
    stat_value <- 0
    for (pop in pops) {
      pooled_transition <- pop$transition[["pooled"]]
      stat_value <- stat_value + .cat_score_timeinvariance_single(
        y = pop$y,
        p = p,
        c = c,
        pooled_transition = pooled_transition
      )
    }
  } else if (identical(test, "mlrt")) {
    n_subjects <- fit_alt$settings$n_subjects
    blocks_eff <- .cat_resolve_blocks(blocks, fit_alt)
    sim_params <- .cat_timeinvariant_sim_params(fit_null)

    simulate_fn <- function() {
      simulate_cat(
        n_subjects = n_subjects,
        n_time = n_time,
        order = p,
        n_categories = c,
        marginal = sim_params$marginal,
        transition = sim_params$transition,
        blocks = blocks_eff,
        homogeneous = homogeneous
      )
    }
    lrt_raw_fn <- function(y_b) {
      fit1_b <- fit_cat(
        y_b,
        order = p,
        blocks = blocks_eff,
        homogeneous = homogeneous,
        n_categories = c,
        na_action = "fail"
      )
      fit0_b <- .fit_cat_timeinvariant(
        y_b,
        p = p,
        c = c,
        blocks = blocks_eff,
        homogeneous = homogeneous
      )
      -2 * (fit0_b$log_l - fit1_b$log_l)
    }

    e_hat_mlrt <- .cat_mlrt_expected_lrt_boot(simulate_fn, lrt_raw_fn)
    if (is.finite(e_hat_mlrt) && e_hat_mlrt > 0) {
      stat_value <- lrt_stat_raw * df / e_hat_mlrt
    } else {
      warning(
        "Modified LRT scaling factor could not be estimated; returning unmodified LRT statistic.",
        call. = FALSE
      )
    }
  }

  # P-value from chi-square
  p_value <- stats::pchisq(stat_value, df = df, lower.tail = FALSE)
  
  # Build summary table
  table_df <- data.frame(
    model = c("Time-invariant (H0)", "Time-varying (H1)"),
    log_l = c(fit_null$log_l, fit_alt$log_l),
    n_params = c(fit_null$n_params, fit_alt$n_params),
    aic = c(fit_null$aic, fit_alt$aic),
    bic = c(fit_null$bic, fit_alt$bic)
  )
  
  # Assemble output
  out <- list(
    method = test,
    lrt_stat = stat_value,
    statistic = stat_value,
    lrt_stat_raw = lrt_stat_raw,
    e_hat_mlrt = e_hat_mlrt,
    df = df,
    p_value = p_value,
    test = test,
    fit_null = fit_null,
    fit_alt = fit_alt,
    order = p,
    table = table_df
  )
  
  class(out) <- "cat_lrt"
  out
}


#' Fit time-invariant categorical AD model
#'
#' Internal function to fit a model where transition probabilities are
#' constant across time (for k > p).
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param blocks Block vector (or NULL)
#' @param homogeneous Whether to pool across blocks
#'
#' @return A cat_fit-like object with time-invariant transitions
#'
#' @keywords internal
.fit_cat_timeinvariant <- function(y, p, c, blocks = NULL, homogeneous = TRUE) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Process blocks
  if (is.null(blocks)) {
    blocks <- rep(1L, n_subjects)
  }
  n_blocks <- max(blocks)
  
  if (homogeneous) {
    n_pops <- 1
  } else {
    n_pops <- n_blocks
  }
  
  # For each population, pool transition counts across all time points
  if (n_pops == 1) {
    result <- .fit_cat_timeinvariant_single(y, p, c, n_time, subject_mask = NULL)
    marginal <- result$marginal
    transition <- result$transition
    log_l <- result$log_l
    n_params <- result$n_params
  } else {
    marginal <- vector("list", n_pops)
    transition <- vector("list", n_pops)
    log_l <- 0
    n_params <- 0
    
    for (g in seq_len(n_pops)) {
      mask <- (blocks == g)
      result_g <- .fit_cat_timeinvariant_single(y, p, c, n_time, subject_mask = mask)
      marginal[[g]] <- result_g$marginal
      transition[[g]] <- result_g$transition
      log_l <- log_l + result_g$log_l
      n_params <- n_params + result_g$n_params
    }
    names(marginal) <- paste0("block_", seq_len(n_pops))
    names(transition) <- paste0("block_", seq_len(n_pops))
  }
  
  # Compute AIC and BIC
  aic <- -2 * log_l + 2 * n_params
  bic <- -2 * log_l + n_params * log(n_subjects)
  
  out <- list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    aic = aic,
    bic = bic,
    n_params = n_params,
    convergence = 0L,
    settings = list(
      order = p,
      n_categories = c,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = if (n_blocks > 1) blocks else NULL,
      homogeneous = homogeneous,
      n_blocks = n_blocks,
      time_invariant = TRUE
    )
  )
  
  class(out) <- "cat_fit"
  out
}


#' Fit time-invariant model for single population
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param subject_mask Logical mask for subjects (NULL = all)
#'
#' @return List with marginal, transition, log_l, n_params
#'
#' @keywords internal
.fit_cat_timeinvariant_single <- function(y, p, c, n_time, subject_mask = NULL) {
  
  # Apply mask
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  N <- nrow(y_sub)
  
  # Storage
  marginal <- list()
  log_l <- 0
  
  # Handle initial time points (same as regular fit)
  for (k in seq_len(min(p, n_time))) {
    time_indices <- seq_len(k)
    counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
    
    if (k == 1) {
      probs_k <- counts_k / N
      marginal[["t1"]] <- as.numeric(probs_k)
      names(marginal[["t1"]]) <- paste0("cat_", 1:c)
    } else {
      probs_k <- .counts_to_transition_probs(counts_k)
      marginal[[paste0("t", k, "_given_1to", k - 1)]] <- probs_k
    }
    
    log_l <- log_l + .loglik_contribution(counts_k, probs_k)
  }
  
  # Pool transition counts across all time points k = p+1 to n
  # Count all (Y_{k-p}, ..., Y_{k-1}, Y_k) patterns together
  pooled_counts <- array(0L, dim = rep(c, p + 1))
  
  for (k in (p + 1):n_time) {
    time_indices <- (k - p):k
    counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
    pooled_counts <- pooled_counts + counts_k
  }
  
  # Compute pooled transition probabilities
  pooled_probs <- .counts_to_transition_probs(pooled_counts)
  
  # Store single transition matrix
  transition <- list(pooled = pooled_probs)
  
  # Compute log-likelihood contribution from transitions
  # Need to evaluate at each time point using the pooled probs
  for (k in (p + 1):n_time) {
    time_indices <- (k - p):k
    counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
    log_l <- log_l + .loglik_contribution(counts_k, pooled_probs)
  }
  
  # Count parameters
  # Marginal: (c-1) for t1, plus (c-1)*c^{k-1} for each k=2,...,p
  # Transition: (c-1)*c^p (just one matrix, not n-p)
  n_params_marginal <- 0
  for (k in seq_len(min(p, n_time))) {
    if (k == 1) {
      n_params_marginal <- n_params_marginal + (c - 1)
    } else {
      n_params_marginal <- n_params_marginal + (c - 1) * c^(k - 1)
    }
  }
  n_params_trans <- (c - 1) * c^p
  n_params <- n_params_marginal + n_params_trans
  
  list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    n_params = n_params
  )
}
