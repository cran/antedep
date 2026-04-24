# loglik_cat.R - Log-likelihood computation for categorical antedependence models

#' Log-likelihood for categorical AD models (with missing data support)
#'
#' Evaluates the log-likelihood of an AD(p) model for categorical longitudinal
#' data at given parameter values.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c.
#' @param order Antedependence order p. Must be 0, 1, or 2.
#' @param marginal List of marginal/joint probabilities for initial time points.
#'   Structure depends on order (see Details).
#' @param transition List of transition probability arrays for time points
#'   k = p+1 to n. Each element should be an array of dimension c^p x c where
#'   the last dimension corresponds to the current time point.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership. Required if homogeneous = FALSE.
#' @param homogeneous Logical. If TRUE (default), same parameters used for all
#'   subjects. If FALSE, marginal and transition should be lists indexed by block.
#' @param n_categories Number of categories. If NULL, inferred from data.
#' @param na_action Handling of missing values in \code{y}. One of
#'   \code{"fail"} (default, error if any missing), \code{"complete"}
#'   (drop subjects with any missing values), or \code{"marginalize"}
#'   (integrate over missing categorical outcomes under the AD model).
#'
#' @return Scalar log-likelihood value.
#'
#' @details
#' The log-likelihood for AD(p) decomposes into contributions from initial
#' time points and transition time points.
#'
#' For order 0 (independence), the log-likelihood is the sum of log marginal
#' probabilities at each time point.
#'
#' Parameter structure for marginal:
#' \itemize{
#'   \item Order 0: List with elements t1, t2, ..., tn, each a vector of length c
#'   \item Order 1: List with element t1 (vector of length c)
#'   \item Order 2: List with t1 (vector), t2_given_1to1 (c x c matrix)
#' }
#'
#' Parameter structure for transition:
#' \itemize{
#'   \item Order 0: Not used (NULL or empty list)
#'   \item Order 1: List with elements t2, t3, ..., tn, each c x c matrix
#'   \item Order 2: List with elements t3, t4, ..., tn, each c x c x c array
#' }
#'
#' @examples
#' set.seed(1)
#' y <- simulate_cat(n_subjects = 40, n_time = 5, order = 1, n_categories = 3)
#' fit <- fit_cat(y, order = 1, n_categories = 3)
#' logL_cat(
#'   y = y,
#'   order = 1,
#'   marginal = fit$marginal,
#'   transition = fit$transition,
#'   n_categories = 3
#' )
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @export
logL_cat <- function(y, order, marginal, transition = NULL, blocks = NULL,
                     homogeneous = TRUE, n_categories = NULL,
                     na_action = c("fail", "complete", "marginalize")) {
  na_action <- match.arg(na_action)
  
  # Validate data
  validated <- .validate_y_cat(y, n_categories, allow_na = (na_action != "fail"))
  y <- validated$y
  c <- validated$n_categories
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  p <- as.integer(order)
  
  # Validate order
  if (p < 0) stop("order must be non-negative")
  if (p >= n_time) stop("order must be less than number of time points")
  if (na_action == "marginalize" && p > 2) {
    stop("na_action = 'marginalize' currently supports order 0, 1, and 2")
  }
  
  # Validate blocks
  blocks <- .validate_blocks_cat(blocks, n_subjects)
  
  # Handle missing values according to na_action
  if (na_action == "complete") {
    keep <- stats::complete.cases(y)
    y <- y[keep, , drop = FALSE]
    blocks <- blocks[keep]
    if (nrow(y) == 0) {
      stop("No complete subjects remain after na_action = 'complete'")
    }
    blocks <- .validate_blocks_cat(blocks, nrow(y))
  }
  
  n_subjects <- nrow(y)
  n_blocks <- max(blocks)
  
  # Determine number of populations
  if (homogeneous) {
    n_pops <- 1
  } else {
    n_pops <- n_blocks
    if (n_pops > 1 && !is.list(marginal[[1]])) {
      stop("For heterogeneous model, marginal and transition must be lists indexed by block")
    }
  }
  
  # Compute log-likelihood
  log_l <- 0
  
  for (s in seq_len(n_subjects)) {
    # Determine which population parameters to use
    if (homogeneous) {
      marg_s <- marginal
      trans_s <- transition
    } else {
      g <- blocks[s]
      marg_s <- marginal[[g]]
      trans_s <- transition[[g]]
    }
    
    # Add contribution from subject s
    log_l <- log_l + .logL_subject_cat_observed(y[s, ], p, c, marg_s, trans_s,
                                                na_action = na_action)
  }
  
  log_l
}


#' Compute log-likelihood contribution from one subject
#'
#' @param y_s Integer vector of length n_time (one subject's data)
#' @param p Order
#' @param c Number of categories
#' @param marginal Marginal parameters for this population
#' @param transition Transition parameters for this population
#'
#' @return Scalar log-likelihood contribution
#'
#' @keywords internal
.logL_subject_cat <- function(y_s, p, c, marginal, transition) {
  n_time <- length(y_s)
  ll <- 0
  
  if (p == 0) {
    # Independence model
    for (k in seq_len(n_time)) {
      prob <- marginal[[k]][y_s[k]]
      ll <- ll + .safe_log(prob)
    }
    return(ll)
  }
  
  # AD(p) model with p >= 1
  
  # Initial time points (k = 1 to p)
  # k = 1: P(Y_1)
  prob_1 <- marginal[["t1"]][y_s[1]]
  ll <- ll + .safe_log(prob_1)
  
  # k = 2 to min(p, n_time): need conditional on full history
  if (p >= 2 && n_time >= 2) {
    for (k in 2:min(p, n_time)) {
      # Get P(Y_k | Y_1, ..., Y_{k-1})
      trans_name <- paste0("t", k, "_given_1to", k - 1)
      trans_k <- marginal[[trans_name]]
      
      # Index into the array: [y_1, y_2, ..., y_{k-1}, y_k]
      idx <- c(y_s[1:(k-1)], y_s[k])
      prob_k <- trans_k[matrix(idx, nrow = 1)]
      ll <- ll + .safe_log(prob_k)
    }
  }
  
  # Time points k = p+1 to n_time: condition only on last p values
  if (n_time > p) {
    for (k in (p + 1):n_time) {
      trans_name <- paste0("t", k)
      trans_k <- transition[[trans_name]]
      
      # Index: [y_{k-p}, ..., y_{k-1}, y_k]
      idx <- y_s[(k - p):k]
      prob_k <- trans_k[matrix(idx, nrow = 1)]
      ll <- ll + .safe_log(prob_k)
    }
  }
  
  ll
}


#' Compute observed-data log-likelihood contribution from one subject
#'
#' @param y_s Integer vector of length n_time (one subject's data; may contain NA)
#' @param p Order
#' @param c Number of categories
#' @param marginal Marginal parameters for this population
#' @param transition Transition parameters for this population
#' @param na_action Missing-data handling mode
#'
#' @return Scalar log-likelihood contribution
#'
#' @keywords internal
.logL_subject_cat_observed <- function(y_s, p, c, marginal, transition, na_action) {
  has_na <- any(is.na(y_s))
  
  if (!has_na) {
    return(.logL_subject_cat(y_s, p, c, marginal, transition))
  }
  
  if (na_action == "fail") {
    stop("Missing data not supported in this call; use na_action = 'complete' or 'marginalize'")
  }
  
  if (na_action == "complete") {
    stop("Internal error: complete-case subject path should not contain missing values")
  }
  
  # na_action == "marginalize"
  if (p == 0) {
    ll <- 0
    n_time <- length(y_s)
    for (k in seq_len(n_time)) {
      if (!is.na(y_s[k])) {
        ll <- ll + log(marginal[[k]][y_s[k]])
      }
    }
    return(ll)
  }
  
  if (p == 1) {
    return(.logL_subject_cat_marginalize_p1(y_s, c, marginal, transition))
  }
  
  if (p == 2) {
    return(.logL_subject_cat_marginalize_p2(y_s, c, marginal, transition))
  }
  
  stop("na_action = 'marginalize' currently supports order 0, 1, and 2")
}


#' Stable log-sum-exp
#'
#' @param x Numeric vector
#'
#' @return Scalar log(sum(exp(x)))
#'
#' @keywords internal
.logsumexp_cat <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}


#' Observed-data log-likelihood for order-1 CAT model
#'
#' @param y_s Subject row (may contain NA)
#' @param c Number of categories
#' @param marginal Marginal list
#' @param transition Transition list
#'
#' @return Scalar log-likelihood contribution
#'
#' @keywords internal
.logL_subject_cat_marginalize_p1 <- function(y_s, c, marginal, transition) {
  n_time <- length(y_s)
  log_alpha <- rep(-Inf, c)
  
  # k = 1
  if (is.na(y_s[1])) {
    log_alpha <- log(marginal[["t1"]])
  } else {
    log_alpha[y_s[1]] <- log(marginal[["t1"]][y_s[1]])
  }
  
  if (n_time == 1) {
    return(.logsumexp_cat(log_alpha))
  }
  
  # k = 2,...,n
  for (k in 2:n_time) {
    trans_k <- transition[[paste0("t", k)]]
    log_trans_k <- log(trans_k)
    log_trans_k[trans_k <= 0] <- -Inf
    
    log_next <- rep(-Inf, c)
    if (is.na(y_s[k])) {
      for (j in seq_len(c)) {
        log_next[j] <- .logsumexp_cat(log_alpha + log_trans_k[, j])
      }
    } else {
      j_obs <- y_s[k]
      log_next[j_obs] <- .logsumexp_cat(log_alpha + log_trans_k[, j_obs])
    }
    log_alpha <- log_next
  }
  
  .logsumexp_cat(log_alpha)
}


#' Observed-data log-likelihood for order-2 CAT model
#'
#' @param y_s Subject row (may contain NA)
#' @param c Number of categories
#' @param marginal Marginal list
#' @param transition Transition list
#'
#' @return Scalar log-likelihood contribution
#'
#' @keywords internal
.logL_subject_cat_marginalize_p2 <- function(y_s, c, marginal, transition) {
  n_time <- length(y_s)
  
  # k = 1
  log_alpha1 <- rep(-Inf, c)
  if (is.na(y_s[1])) {
    log_alpha1 <- log(marginal[["t1"]])
  } else {
    log_alpha1[y_s[1]] <- log(marginal[["t1"]][y_s[1]])
  }
  
  if (n_time == 1) {
    return(.logsumexp_cat(log_alpha1))
  }
  
  # k = 2: state is (Y1, Y2)
  trans_2 <- marginal[["t2_given_1to1"]]
  log_trans_2 <- log(trans_2)
  log_trans_2[trans_2 <= 0] <- -Inf
  log_state <- matrix(-Inf, nrow = c, ncol = c)
  
  for (a in seq_len(c)) {
    if (!is.finite(log_alpha1[a])) next
    if (!is.na(y_s[1]) && y_s[1] != a) next
    for (b in seq_len(c)) {
      if (!is.na(y_s[2]) && y_s[2] != b) next
      log_state[a, b] <- log_alpha1[a] + log_trans_2[a, b]
    }
  }
  
  if (n_time == 2) {
    return(.logsumexp_cat(as.vector(log_state)))
  }
  
  # k = 3,...,n: state is (Y_{k-1}, Y_k)
  for (k in 3:n_time) {
    trans_k <- transition[[paste0("t", k)]]
    log_trans_k <- log(trans_k)
    log_trans_k[trans_k <= 0] <- -Inf
    log_next <- matrix(-Inf, nrow = c, ncol = c)
    
    for (b in seq_len(c)) {
      for (cat_k in seq_len(c)) {
        if (!is.na(y_s[k]) && y_s[k] != cat_k) next
        
        candidates <- rep(-Inf, c)
        for (a in seq_len(c)) {
          if (!is.finite(log_state[a, b])) next
          candidates[a] <- log_state[a, b] + log_trans_k[a, b, cat_k]
        }
        log_next[b, cat_k] <- .logsumexp_cat(candidates)
      }
    }
    
    log_state <- log_next
  }
  
  .logsumexp_cat(as.vector(log_state))
}


#' Compute log-likelihood from cell counts (faster for large datasets)
#'
#' Instead of looping over subjects, this version uses aggregated cell counts.
#' Equivalent to logL_cat but more efficient.
#'
#' @param y Data matrix
#' @param order AD order p
#' @param marginal Marginal parameters
#' @param transition Transition parameters
#' @param n_categories Number of categories
#'
#' @return Scalar log-likelihood
#'
#' @keywords internal
.logL_from_counts_cat <- function(y, order, marginal, transition, n_categories) {
  c <- n_categories
  n_time <- ncol(y)
  p <- order
  
  log_l <- 0
  
  if (p == 0) {
    # Independence: sum over each time point
    for (k in seq_len(n_time)) {
      counts_k <- .count_cells_table_cat(y, k, c)
      probs_k <- marginal[[k]]
      log_l <- log_l + sum(counts_k * .safe_log(probs_k))
    }
    return(log_l)
  }
  
  # AD(p) with p >= 1
  
  # Initial time points
  for (k in seq_len(min(p, n_time))) {
    time_indices <- seq_len(k)
    counts_k <- .count_cells_table_cat(y, time_indices, c)
    
    if (k == 1) {
      probs_k <- marginal[["t1"]]
      log_l <- log_l + sum(counts_k * .safe_log(probs_k))
    } else {
      trans_name <- paste0("t", k, "_given_1to", k - 1)
      probs_k <- marginal[[trans_name]]
      log_l <- log_l + sum(counts_k * .safe_log(probs_k))
    }
  }
  
  # Transition time points
  for (k in (p + 1):n_time) {
    time_indices <- (k - p):k
    counts_k <- .count_cells_table_cat(y, time_indices, c)
    trans_name <- paste0("t", k)
    probs_k <- transition[[trans_name]]
    log_l <- log_l + sum(counts_k * .safe_log(probs_k))
  }
  
  log_l
}
