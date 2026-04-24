# fit_cat.R - Maximum likelihood estimation for categorical antedependence models

#' Fit categorical antedependence model by maximum likelihood
#'
#' Computes maximum likelihood estimates for the parameters of an AD(p) model
#' for categorical longitudinal data. The model is parameterized by transition
#' probabilities, and MLEs are obtained in closed form.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c, where c is the number of categories.
#' @param order Antedependence order p. Must be 0, 1, or 2. Default is 1.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership. If NULL, all subjects are treated as one group.
#' @param homogeneous Logical. If TRUE (default), parameters are shared across
#'   all groups (blocks are ignored for estimation). If FALSE, separate
#'   transition probabilities are estimated for each group.
#' @param n_categories Number of categories. If NULL (default), inferred from
#'   the maximum value in y.
#' @param na_action Handling of missing values in \code{y}. One of
#'   \code{"fail"} (default, error if any missing), \code{"complete"}
#'   (drop subjects with any missing values), or \code{"marginalize"}
#'   (maximize observed-data likelihood by integrating over missing outcomes),
#'   or \code{"em"} (use \code{\link{em_cat}} for orders 0 and 1).
#' @param em_max_iter Maximum EM iterations used when \code{na_action = "em"}.
#' @param em_tol EM convergence tolerance used when \code{na_action = "em"}.
#' @param em_epsilon Numerical smoothing constant used when \code{na_action = "em"}.
#' @param em_safeguard Logical; if \code{TRUE}, use step-halving safeguard in
#'   \code{\link{em_cat}} when \code{na_action = "em"}.
#' @param em_verbose Logical; print EM progress when \code{na_action = "em"}.
#'
#' @return A list of class \code{"cat_fit"} containing:
#'   \item{marginal}{List of marginal/joint probabilities for initial time points}
#'   \item{transition}{List of transition probability arrays for k = p+1 to n}
#'   \item{log_l}{Log-likelihood at MLE}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{n_params}{Number of free parameters}
#'   \item{cell_counts}{List of cell counts: observed counts for closed-form fits
#'   (\code{na_action = "fail"/"complete"}), expected counts from the final E-step
#'   for EM fits (\code{na_action = "em"}), and \code{NULL} for
#'   \code{na_action = "marginalize"}}
#'   \item{convergence}{Optimizer convergence code (0 for closed-form solutions)}
#'   \item{settings}{List of model settings}
#'
#' @details
#' For AD(p), the model decomposes as:
#' \deqn{P(Y_1, \ldots, Y_n) = P(Y_1, \ldots, Y_p) \times \prod_{k=p+1}^{n} P(Y_k | Y_{k-p}, \ldots, Y_{k-1})}
#'
#' MLEs are computed as empirical proportions:
#' \itemize{
#'   \item Marginal/joint probabilities: count / N
#'   \item Transition probabilities: conditional count / marginal count
#' }
#'
#' Empty cells receive probability 0 (if denominator is also 0).
#'
#' When \code{na_action = "em"}, \code{fit_cat()} dispatches to
#' \code{\link{em_cat}}. In that case, \code{em_safeguard} controls step-halving
#' protection against likelihood-decreasing updates, and returned
#' \code{log_l}/AIC/BIC/\code{cell_counts} are synchronized via a final E-step
#' under the returned parameters.
#' For \code{order = 2}, \code{na_action = "em"} is not available and errors
#' explicitly; use \code{na_action = "marginalize"}.
#'
#' @examples
#' \donttest{
#' # Simulate binary AD(1) data
#' set.seed(123)
#' y <- simulate_cat(n_subjects = 100, n_time = 5, order = 1, n_categories = 2)
#'
#' # Fit model
#' fit <- fit_cat(y, order = 1)
#' print(fit)
#'
#' # Compare orders
#' fit0 <- fit_cat(y, order = 0)
#' fit1 <- fit_cat(y, order = 1)
#' fit2 <- fit_cat(y, order = 2)
#' c(AIC_0 = fit0$aic, AIC_1 = fit1$aic, AIC_2 = fit2$aic)
#'
#' # EM fit with missing data
#' y_miss <- y
#' y_miss[sample(length(y_miss), size = round(0.15 * length(y_miss)))] <- NA
#' fit_em <- fit_cat(
#'   y_miss,
#'   order = 1,
#'   na_action = "em",
#'   em_max_iter = 80,
#'   em_tol = 1e-6
#' )
#' fit_em$settings$n_iter
#' fit_em$settings$cell_counts_type
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @export
fit_cat <- function(y, order = 1, blocks = NULL, homogeneous = TRUE, 
                    n_categories = NULL,
                    na_action = c("fail", "complete", "marginalize", "em"),
                    em_max_iter = 100, em_tol = 1e-6, em_epsilon = 1e-8,
                    em_safeguard = TRUE, em_verbose = FALSE) {
  na_action <- match.arg(na_action)
  blocks_input <- blocks
  

  # Validate inputs
  validated <- .validate_y_cat(y, n_categories, allow_na = (na_action != "fail"))
  y <- validated$y
  c <- validated$n_categories
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  p <- as.integer(order)
  
  # Validate order
  if (p < 0) stop("order must be non-negative")
  if (p > 2) stop("order > 2 not currently supported")
  if (p >= n_time) stop("order must be less than number of time points")
  if (na_action == "marginalize" && p > 2) {
    stop("na_action = 'marginalize' currently supports order 0, 1, and 2")
  }
  
  # Validate and process blocks
  block_info <- .normalize_blocks(blocks, n_subjects)
  blocks <- block_info$blocks_id
  block_levels <- block_info$block_levels

  if (na_action == "em") {
    if (p == 2L) {
      stop("na_action = 'em' currently supports order 0 and 1; use na_action = 'marginalize' for order 2")
    }
    return(em_cat(
      y = y,
      order = p,
      blocks = blocks_input,
      homogeneous = homogeneous,
      n_categories = c,
      max_iter = em_max_iter,
      tol = em_tol,
      epsilon = em_epsilon,
      safeguard = em_safeguard,
      verbose = em_verbose
    ))
  }
  
  # Handle complete-case filtering if requested
  if (na_action == "complete") {
    keep <- stats::complete.cases(y)
    y <- y[keep, , drop = FALSE]
    if (!is.null(blocks_input)) {
      block_info <- .normalize_blocks(blocks_input[keep], nrow(y))
      blocks <- block_info$blocks_id
      block_levels <- block_info$block_levels
    } else {
      blocks <- rep(1L, nrow(y))
      block_levels <- "1"
    }
    if (nrow(y) == 0) {
      stop("No complete subjects remain after na_action = 'complete'")
    }
  }
  
  n_subjects <- nrow(y)
  effective_na_action <- na_action
  if (na_action == "marginalize" && !any(is.na(y))) {
    # Degenerate to complete-data closed-form fit when no missing values remain.
    effective_na_action <- "complete"
  }
  
  n_blocks <- max(blocks)
  
  # Determine number of populations to estimate
  if (homogeneous) {
    n_pops <- 1
  } else {
    n_pops <- n_blocks
  }
  
  # Count parameters
  n_params <- .count_params_cat(p, c, n_time, n_blocks, homogeneous)
  
  convergence <- 0L
  
  # Initialize storage for parameters
  if (n_pops == 1) {
    if (effective_na_action == "marginalize") {
      result <- .fit_cat_single_pop_marginalize(y, p, c, n_time)
      convergence <- as.integer(result$convergence)
    } else {
      result <- .fit_cat_single_pop(y, p, c, n_time, subject_mask = NULL)
    }
    marginal <- result$marginal
    transition <- result$transition
    cell_counts <- result$cell_counts
    log_l <- result$log_l
  } else {
    # Fit separately for each block
    marginal <- vector("list", n_pops)
    transition <- vector("list", n_pops)
    cell_counts <- vector("list", n_pops)
    log_l <- 0
    conv_codes <- integer(n_pops)
    
    for (g in seq_len(n_pops)) {
      mask <- (blocks == g)
      y_g <- y[mask, , drop = FALSE]
      if (effective_na_action == "marginalize") {
        result_g <- .fit_cat_single_pop_marginalize(y_g, p, c, n_time)
        conv_codes[g] <- as.integer(result_g$convergence)
      } else {
        result_g <- .fit_cat_single_pop(y, p, c, n_time, subject_mask = mask)
        conv_codes[g] <- 0L
      }
      marginal[[g]] <- result_g$marginal
      transition[[g]] <- result_g$transition
      # Preserve list length even when cell_counts is NULL under marginalization.
      cell_counts[g] <- list(result_g$cell_counts)
      log_l <- log_l + result_g$log_l
    }
    names(marginal) <- paste0("block_", seq_len(n_pops))
    names(transition) <- paste0("block_", seq_len(n_pops))
    names(cell_counts) <- paste0("block_", seq_len(n_pops))
    convergence <- as.integer(max(conv_codes))
  }
  
  # Compute AIC and BIC
  aic <- -2 * log_l + 2 * n_params
  bic <- -2 * log_l + n_params * log(n_subjects)
  
  # Assemble output
  out <- list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    aic = aic,
    bic = bic,
    n_params = n_params,
    n_obs = sum(!is.na(y)),
    n_missing = sum(is.na(y)),
    pct_missing = mean(is.na(y)) * 100,
    cell_counts = cell_counts,
    convergence = convergence,
    settings = list(
      order = p,
      n_categories = c,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = if (n_blocks > 1) blocks else NULL,
      block_levels = block_levels,
      homogeneous = homogeneous,
      n_blocks = n_blocks,
      na_action = na_action,
      na_action_effective = effective_na_action,
      cell_counts_type = if (identical(effective_na_action, "marginalize")) "none" else "observed"
    )
  )
  
  class(out) <- "cat_fit"
  out
}


#' Fit CAT model with missing data via observed-data likelihood optimization
#'
#' @param y Data matrix for one population (may contain NA)
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#'
#' @return List with marginal, transition, cell_counts, log_l, convergence
#'
#' @keywords internal
.fit_cat_single_pop_marginalize <- function(y, p, c, n_time) {
  # Initialize from complete cases if available; otherwise uniform probabilities.
  y_complete <- y[stats::complete.cases(y), , drop = FALSE]
  if (nrow(y_complete) > 0) {
    init_fit <- .fit_cat_single_pop(y_complete, p, c, n_time, subject_mask = NULL)
    init_marginal <- init_fit$marginal
    init_transition <- init_fit$transition
  } else {
    init <- .uniform_cat_params(p, c, n_time)
    init_marginal <- init$marginal
    init_transition <- init$transition
  }
  
  theta_start <- .pack_cat_params(init_marginal, init_transition, p, c, n_time)
  
  objective <- function(theta) {
    params <- .unpack_cat_params(theta, p, c, n_time)
    ll <- logL_cat(
      y = y,
      order = p,
      marginal = params$marginal,
      transition = params$transition,
      n_categories = c,
      na_action = "marginalize"
    )
    if (!is.finite(ll)) return(1e50)
    -as.numeric(ll)
  }
  
  if (length(theta_start) == 0) {
    # Degenerate single-category case
    params <- .unpack_cat_params(theta_start, p, c, n_time)
    ll <- -objective(theta_start)
    return(list(
      marginal = params$marginal,
      transition = params$transition,
      cell_counts = NULL,
      log_l = ll,
      convergence = 0L
    ))
  }
  
  opt <- stats::optim(
    par = theta_start,
    fn = objective,
    method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-10)
  )
  
  # Fallback if BFGS does not converge cleanly.
  if (opt$convergence != 0L || !is.finite(opt$value)) {
    opt_nm <- stats::optim(
      par = opt$par,
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 4000, reltol = 1e-10)
    )
    if (is.finite(opt_nm$value) && (!is.finite(opt$value) || opt_nm$value < opt$value)) {
      opt <- opt_nm
    }
  }
  
  params <- .unpack_cat_params(opt$par, p, c, n_time)
  
  list(
    marginal = params$marginal,
    transition = params$transition,
    cell_counts = NULL,
    log_l = -opt$value,
    convergence = as.integer(opt$convergence)
  )
}


#' Convert probability row to unconstrained logits
#'
#' @param prob Probability vector of length c
#' @param c Number of categories
#' @param eps Lower bound for numerical stability
#'
#' @return Numeric vector of length c-1
#'
#' @keywords internal
.cat_prob_to_theta <- function(prob, c, eps = 1e-8) {
  if (c <= 1) return(numeric(0))
  prob <- as.numeric(prob)
  if (length(prob) != c) {
    stop("Probability row has wrong length in CAT parameter packing")
  }
  prob <- pmax(prob, eps)
  prob <- prob / sum(prob)
  log(prob[seq_len(c - 1)] / prob[c])
}


#' Convert unconstrained logits to probability row
#'
#' @param theta_row Numeric vector of length c-1
#' @param c Number of categories
#'
#' @return Probability vector of length c
#'
#' @keywords internal
.cat_theta_to_prob <- function(theta_row, c) {
  if (c <= 1) return(1)
  if (length(theta_row) != (c - 1)) {
    stop("Theta row has wrong length in CAT parameter unpacking")
  }
  z <- c(theta_row, 0)
  z <- z - max(z)
  ez <- exp(z)
  ez / sum(ez)
}


#' Pack CAT parameters into unconstrained vector
#'
#' @param marginal Marginal parameter list
#' @param transition Transition parameter list
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#'
#' @return Numeric parameter vector
#'
#' @keywords internal
.pack_cat_params <- function(marginal, transition, p, c, n_time) {
  theta <- numeric(0)
  
  if (p == 0) {
    for (k in seq_len(n_time)) {
      theta <- c(theta, .cat_prob_to_theta(marginal[[k]], c))
    }
    return(theta)
  }
  
  theta <- c(theta, .cat_prob_to_theta(marginal[["t1"]], c))
  
  if (p == 1) {
    for (k in 2:n_time) {
      trans_k <- transition[[paste0("t", k)]]
      for (from in seq_len(c)) {
        theta <- c(theta, .cat_prob_to_theta(trans_k[from, ], c))
      }
    }
    return(theta)
  }
  
  # p == 2
  trans_t2 <- marginal[["t2_given_1to1"]]
  for (from in seq_len(c)) {
    theta <- c(theta, .cat_prob_to_theta(trans_t2[from, ], c))
  }
  for (k in 3:n_time) {
    trans_k <- transition[[paste0("t", k)]]
    for (from1 in seq_len(c)) {
      for (from2 in seq_len(c)) {
        theta <- c(theta, .cat_prob_to_theta(trans_k[from1, from2, ], c))
      }
    }
  }
  
  theta
}


#' Unpack unconstrained CAT parameter vector
#'
#' @param theta Numeric parameter vector
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#'
#' @return List with marginal and transition parameter lists
#'
#' @keywords internal
.unpack_cat_params <- function(theta, p, c, n_time) {
  idx <- 1L
  
  next_row <- function() {
    if (c <= 1) return(1)
    end_idx <- idx + c - 2
    if (end_idx > length(theta)) {
      stop("Theta vector ended early while unpacking CAT parameters")
    }
    row <- .cat_theta_to_prob(theta[idx:end_idx], c)
    idx <<- end_idx + 1L
    row
  }
  
  marginal <- list()
  transition <- list()
  
  if (p == 0) {
    for (k in seq_len(n_time)) {
      marginal[[k]] <- next_row()
      names(marginal[[k]]) <- paste0("cat_", seq_len(c))
    }
    names(marginal) <- paste0("t", seq_len(n_time))
  } else if (p == 1) {
    marginal[["t1"]] <- next_row()
    names(marginal[["t1"]]) <- paste0("cat_", seq_len(c))
    
    for (k in 2:n_time) {
      trans_k <- matrix(0, nrow = c, ncol = c)
      for (from in seq_len(c)) {
        trans_k[from, ] <- next_row()
      }
      transition[[paste0("t", k)]] <- trans_k
    }
  } else if (p == 2) {
    marginal[["t1"]] <- next_row()
    names(marginal[["t1"]]) <- paste0("cat_", seq_len(c))
    
    trans_t2 <- matrix(0, nrow = c, ncol = c)
    for (from in seq_len(c)) {
      trans_t2[from, ] <- next_row()
    }
    marginal[["t2_given_1to1"]] <- trans_t2
    
    for (k in 3:n_time) {
      trans_k <- array(0, dim = c(c, c, c))
      for (from1 in seq_len(c)) {
        for (from2 in seq_len(c)) {
          trans_k[from1, from2, ] <- next_row()
        }
      }
      transition[[paste0("t", k)]] <- trans_k
    }
  } else {
    stop("Only orders 0, 1, and 2 are supported")
  }
  
  if (idx <= length(theta)) {
    stop("Theta vector has extra elements after CAT parameter unpacking")
  }
  
  list(marginal = marginal, transition = transition)
}


#' Build uniform CAT parameter values
#'
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#'
#' @return List with marginal and transition
#'
#' @keywords internal
.uniform_cat_params <- function(p, c, n_time) {
  probs <- rep(1 / c, c)
  marginal <- list()
  transition <- list()
  
  if (p == 0) {
    marginal <- lapply(seq_len(n_time), function(k) {
      out <- probs
      names(out) <- paste0("cat_", seq_len(c))
      out
    })
    names(marginal) <- paste0("t", seq_len(n_time))
    return(list(marginal = marginal, transition = transition))
  }
  
  marginal[["t1"]] <- stats::setNames(probs, paste0("cat_", seq_len(c)))
  
  if (p == 1) {
    for (k in 2:n_time) {
      transition[[paste0("t", k)]] <- matrix(rep(probs, each = c), nrow = c, byrow = TRUE)
    }
  } else if (p == 2) {
    marginal[["t2_given_1to1"]] <- matrix(rep(probs, each = c), nrow = c, byrow = TRUE)
    for (k in 3:n_time) {
      transition[[paste0("t", k)]] <- array(
        rep(probs, each = c * c),
        dim = c(c, c, c)
      )
    }
  }
  
  list(marginal = marginal, transition = transition)
}


#' Fit CAT model for a single population
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param subject_mask Logical vector for subject selection (NULL = all)
#'
#' @return List with marginal, transition, cell_counts, log_l
#'
#' @keywords internal
.fit_cat_single_pop <- function(y, p, c, n_time, subject_mask = NULL) {
  
  # Apply mask
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  N <- nrow(y_sub)
  
  # Storage
  marginal <- list()
  transition <- list()
  cell_counts <- list()
  log_l <- 0
  
  if (p == 0) {
    # Independence model: just estimate marginal at each time point
    for (k in seq_len(n_time)) {
      counts_k <- .count_cells_table_cat(y_sub, k, c, subject_mask = NULL)
      probs_k <- counts_k / N
      
      marginal[[k]] <- as.numeric(probs_k)
      names(marginal[[k]]) <- paste0("cat_", 1:c)
      cell_counts[[paste0("t", k)]] <- counts_k
      
      # Log-likelihood contribution
      log_l <- log_l + .loglik_contribution(counts_k, probs_k)
    }
    names(marginal) <- paste0("t", seq_len(n_time))
    
  } else {
    # AD(p) model with p >= 1
    
    # Handle initial time points (k = 1 to p)
    # For these, we need the full joint distribution P(Y_1, ..., Y_k)
    
    for (k in seq_len(min(p, n_time))) {
      time_indices <- seq_len(k)
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      
      if (k == 1) {
        # Marginal of Y_1
        probs_k <- counts_k / N
        marginal[["t1"]] <- as.numeric(probs_k)
        names(marginal[["t1"]]) <- paste0("cat_", 1:c)
      } else {
        # P(Y_k | Y_1, ..., Y_{k-1})
        # counts_k has dimensions c^k
        # Need to compute conditional probabilities
        probs_k <- .counts_to_transition_probs(counts_k)
        marginal[[paste0("t", k, "_given_1to", k-1)]] <- probs_k
      }
      
      cell_counts[[paste0("t1_to_t", k)]] <- counts_k
      
      # Log-likelihood contribution for this time point
      # For k=1: sum N_y1 * log(pi_y1)
      # For k>1: sum N_{y1...yk} * log(pi_{yk|y1...y_{k-1}})
      if (k == 1) {
        log_l <- log_l + .loglik_contribution(counts_k, probs_k)
      } else {
        # Need joint counts and transition probs
        # LL contribution = sum over all (y1,...,yk) of N * log(pi_{yk|y1...y_{k-1}})
        log_l <- log_l + .loglik_contribution(counts_k, probs_k)
      }
    }
    
    # Handle time points k = p+1 to n (only condition on last p values)
    for (k in (p + 1):n_time) {
      # Count (Y_{k-p}, ..., Y_{k-1}, Y_k) combinations
      time_indices <- (k - p):k
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      
      # Compute transition probabilities P(Y_k | Y_{k-p}, ..., Y_{k-1})
      probs_k <- .counts_to_transition_probs(counts_k)
      
      transition[[paste0("t", k)]] <- probs_k
      cell_counts[[paste0("t", k-p, "_to_t", k)]] <- counts_k
      
      # Log-likelihood contribution
      log_l <- log_l + .loglik_contribution(counts_k, probs_k)
    }
  }
  
  list(
    marginal = marginal,
    transition = transition,
    cell_counts = cell_counts,
    log_l = log_l
  )
}
