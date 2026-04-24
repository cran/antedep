#' EM algorithm for AD with missing data
#'
#' Implements the EM algorithm for fitting AD(p) models with missing data
#' under MAR (Missing At Random) assumption.
#'
#' @param y Data matrix (n_subjects x n_time), may contain NA
#' @param order Antedependence order (0, 1, or 2)
#' @param blocks Block membership vector (optional)
#' @param estimate_mu Logical, whether to estimate mu (default TRUE)
#' @param max_iter Maximum number of EM iterations (default 100)
#' @param tol Convergence tolerance on log-likelihood (default 1e-6)
#' @param verbose Logical, print iteration info (default FALSE)
#'
#' @return List with fitted model (same structure as fit_gau, plus EM info)
#'
#' @keywords internal
.fit_gau_em <- function(y, order, blocks = NULL, estimate_mu = TRUE,
                       max_iter = 100, tol = 1e-6, verbose = FALSE, ...) {
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Validate missing data
  missing_info <- .validate_missing(y)
  
  if (verbose) {
    cat("\n=== EM Algorithm for AD with Missing Data ===\n")
    cat("Missing data: ", missing_info$pct_missing, "% (", 
        sum(is.na(y)), " values)\n", sep = "")
    cat("Complete subjects: ", missing_info$n_complete, "/", n_subjects, "\n")
    cat("Intermittent missing: ", missing_info$n_intermittent, " subjects\n\n")
  }
  
  # Initialize parameters from complete cases or simple estimates
  init <- .initialize_gau_em(y, order, blocks, estimate_mu)
  mu <- init$mu
  phi <- init$phi
  sigma <- init$sigma
  tau <- init$tau
  
  # Track log-likelihood
  ll_trace <- numeric(max_iter)
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    # E-step: Compute sufficient statistics
    suff_stats <- .em_e_step_gau(y, order, mu, phi, sigma, blocks, tau)
    
    # M-step: Update parameters
    params <- .em_m_step_gau(suff_stats, order, blocks, estimate_mu, n_subjects, n_time)
    mu <- params$mu
    phi <- params$phi
    sigma <- params$sigma
    tau <- params$tau
    
    # Compute observed-data log-likelihood
    ll <- .logL_gau_missing(y, order, mu, phi, sigma, blocks, tau)
    ll_trace[iter] <- ll
    
    # Check for invalid log-likelihood
    if (!is.finite(ll)) {
      warning("Non-finite log-likelihood at iteration ", iter, 
              ". EM terminated.", call. = FALSE)
      converged <- FALSE
      break
    }
    
    if (verbose && (iter %% 10 == 0 || iter == 1)) {
      cat("Iteration ", iter, ": log-lik = ", round(ll, 4), "\n", sep = "")
    }
    
    # Check convergence
    if (iter > 1) {
      if (abs(ll_trace[iter] - ll_trace[iter - 1]) < tol) {
        converged <- TRUE
        if (verbose) {
          cat("\nConverged at iteration ", iter, "\n")
        }
        break
      }
      
      # Check for non-monotonic increase (should not happen with EM)
      if (ll_trace[iter] < ll_trace[iter - 1] - 1e-6) {
        warning("EM log-likelihood decreased at iteration ", iter,
                ". This should not happen.", call. = FALSE)
      }
    }
  }
  
  if (!converged && verbose) {
    cat("\nWarning: EM did not converge in ", max_iter, " iterations\n")
  }
  
  # Compute final quantities
  n_blocks <- if (is.null(blocks)) 1 else length(unique(blocks))
  n_params <- .count_params_gau(order, n_time, n_blocks)
  
  # Return fitted model
  result <- list(
    mu = mu,
    phi = phi,
    sigma = sigma,
    tau = if (n_blocks > 1) tau else NULL,
    log_l = ll,
    aic = -2 * ll + 2 * n_params,
    bic = -2 * ll + n_params * log(n_subjects),
    n_params = as.integer(n_params),
    convergence = if (converged) 0L else 1L,
    message = if (converged) "EM converged" else "EM did not converge",
    
    # Missing data info
    n_obs = sum(!is.na(y)),
    n_missing = sum(is.na(y)),
    pct_missing = missing_info$pct_missing,
    missing_pattern = if (missing_info$n_intermittent > 0) "intermittent" 
                      else "monotone",
    
    # EM convergence info
    em_converged = converged,
    em_iterations = iter,
    em_ll_trace = ll_trace[1:iter],
    
    settings = list(
      order = order,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = blocks,
      estimate_mu = estimate_mu
    )
  )
  
  class(result) <- c("gau_fit", "list")
  result
}

#' Initialize parameters for AD EM algorithm
#'
#' @keywords internal
.initialize_gau_em <- function(y, order, blocks, estimate_mu) {
  # Try complete-case initialization
  cc <- .extract_complete_cases(y, blocks, warn = FALSE)
  
  if (nrow(cc$y) >= 10) {
    # Enough complete cases - initialize from the complete-case MLE.
    tryCatch({
      fit_cc <- fit_gau(
        y = cc$y,
        order = order,
        blocks = cc$blocks,
        na_action = "fail",
        estimate_mu = estimate_mu
      )

      tau <- fit_cc$tau
      if (is.null(tau)) tau <- 0

      list(
        mu = fit_cc$mu,
        phi = fit_cc$phi,
        sigma = fit_cc$sigma,
        tau = tau
      )
    }, error = function(e) {
      # Fallback to marginal estimates
      .initialize_gau_marginal(y, order, blocks)
    })
  } else {
    # Too few complete cases - use marginal estimates
    .initialize_gau_marginal(y, order, blocks)
  }
}

#' Initialize from marginal statistics (ignoring dependence)
#'
#' @keywords internal
.initialize_gau_marginal <- function(y, order, blocks) {
  n_time <- ncol(y)
  
  mu <- colMeans(y, na.rm = TRUE)
  sigma <- apply(y, 2, sd, na.rm = TRUE)
  sigma[sigma == 0 | is.na(sigma)] <- 1  # Avoid zero sigma
  
  # Default phi: moderate positive dependence
  if (order == 1) {
    phi <- rep(0.3, n_time - 1)
  } else if (order == 2) {
    phi <- matrix(0, nrow = 2, ncol = n_time)
    if (n_time >= 2) phi[1, 2:n_time] <- 0.25
    if (n_time >= 3) phi[2, 3:n_time] <- 0.10
  } else {
    phi <- numeric(0)
  }
  
  if (!is.null(blocks)) {
    n_blocks <- length(unique(blocks))
    tau <- numeric(n_blocks)
    tau[1] <- 0
    
    # Simple block means
    if (n_blocks > 1) {
      for (b in 2:n_blocks) {
        block_data <- y[blocks == b, , drop = FALSE]
        tau[b] <- mean(block_data, na.rm = TRUE) - mean(y, na.rm = TRUE)
      }
    }
  } else {
    tau <- 0
  }
  
  list(mu = mu, phi = phi, sigma = sigma, tau = tau)
}

#' E-step: Compute expected sufficient statistics
#'
#' Computes conditional expectations for missing values and second moments.
#'
#' @keywords internal
.em_e_step_gau <- function(y, order, mu, phi, sigma, blocks, tau) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Build full covariance matrix
  Sigma <- .build_gau_covariance(order, phi, sigma, n_time)
  
  # Adjust mean for blocks if present
  if (!is.null(blocks)) {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
    tau_matrix <- matrix(tau[blocks], nrow = n_subjects, ncol = n_time)
    mu_matrix <- mu_matrix + tau_matrix
  } else {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
  }
  
  # Initialize sufficient statistics
  # S1[t] = sum_s E[Y_t^(s)]
  # S2[t,t'] = sum_s E[Y_t^(s) * Y_t'^(s)]
  S1 <- numeric(n_time)
  S2 <- matrix(0, n_time, n_time)
  
  for (s in 1:n_subjects) {
    y_s <- y[s, ]
    mu_s <- mu_matrix[s, ]
    
    obs_idx <- which(!is.na(y_s))
    mis_idx <- which(is.na(y_s))
    
    if (length(mis_idx) == 0) {
      # Complete case
      y_full <- y_s
      S2_s <- tcrossprod(y_s)  # y * y'
      
    } else {
      # Has missing values
      y_obs <- y_s[obs_idx]
      mu_obs <- mu_s[obs_idx]
      mu_mis <- mu_s[mis_idx]
      
      Sigma_obs <- Sigma[obs_idx, obs_idx, drop = FALSE]
      Sigma_mis <- Sigma[mis_idx, mis_idx, drop = FALSE]
      Sigma_cross <- Sigma[mis_idx, obs_idx, drop = FALSE]
      
      # Conditional mean: E[Y_mis | Y_obs]
      Sigma_obs_inv <- tryCatch(
        solve(Sigma_obs),
        error = function(e) {
          solve(Sigma_obs + diag(1e-8, nrow(Sigma_obs)))
        }
      )
      E_mis <- mu_mis + Sigma_cross %*% Sigma_obs_inv %*% (y_obs - mu_obs)
      
      # Conditional covariance: Var[Y_mis | Y_obs]
      Var_mis <- Sigma_mis - Sigma_cross %*% Sigma_obs_inv %*% t(Sigma_cross)
      Var_mis <- 0.5 * (Var_mis + t(Var_mis))
      
      # Impute missing values
      y_full <- y_s
      y_full[mis_idx] <- E_mis
      
      # Second moment matrix
      S2_s <- tcrossprod(y_full)
      
      # Add conditional variance for missing-missing pairs
      S2_s[mis_idx, mis_idx] <- S2_s[mis_idx, mis_idx] + Var_mis
    }
    
    # Accumulate sufficient statistics
    S1 <- S1 + y_full
    S2 <- S2 + S2_s
  }
  
  list(S1 = S1, S2 = S2, n = n_subjects)
}

#' M-step: Update parameters from sufficient statistics
#'
#' @keywords internal
.em_m_step_gau <- function(suff_stats, order, blocks, estimate_mu, 
                          n_subjects, n_time) {
  S1 <- suff_stats$S1
  S2 <- suff_stats$S2
  n <- suff_stats$n
  
  # Update mu (marginal means)
  if (estimate_mu) {
    mu <- S1 / n
  } else {
    mu <- rep(0, n_time)  # Assume centered if not estimating
  }
  
  # Update block effects if present
  # (Simplified - full implementation would need block-specific sufficient stats)
  if (!is.null(blocks)) {
    n_blocks <- length(unique(blocks))
    tau <- numeric(n_blocks)
    tau[1] <- 0
    # TODO: Implement proper block effect estimation
  } else {
    tau <- 0
  }
  
  # Compute centered second moments
  S2_centered <- S2 / n - tcrossprod(mu)
  S2_centered <- 0.5 * (S2_centered + t(S2_centered))
  eps <- 1e-8
  
  # Update phi and sigma based on order
  if (order == 0) {
    phi <- numeric(0)
    sigma <- sqrt(pmax(diag(S2_centered), eps))
    
  } else if (order == 1) {
    phi <- numeric(n_time - 1)
    sigma <- numeric(n_time)
    sigma[1] <- sqrt(max(S2_centered[1, 1], eps))

    if (n_time >= 2) for (t in 2:n_time) {
      # phi[t] = Cov(Y_t, Y_{t-1}) / Var(Y_{t-1})
      denom <- max(S2_centered[t - 1, t - 1], eps)
      phi[t - 1] <- S2_centered[t, t - 1] / denom
      
      # sigma[t]^2 = Var(Y_t) - phi[t]^2 * Var(Y_{t-1})
      sigma2_t <- S2_centered[t, t] - phi[t - 1] * S2_centered[t, t - 1]
      sigma[t] <- sqrt(max(sigma2_t, eps))
    }
    
  } else if (order == 2) {
    # Exact order-2 conditional least-squares updates from expected moments.
    phi <- matrix(0, nrow = 2, ncol = n_time)
    sigma2 <- pmax(diag(S2_centered), eps)

    # Time 2 depends on time 1 only.
    if (n_time >= 2) {
      var1 <- max(S2_centered[1, 1], eps)
      cov21 <- S2_centered[2, 1]
      phi[1, 2] <- cov21 / var1
      sigma2[2] <- max(S2_centered[2, 2] - phi[1, 2] * cov21, eps)
    }

    # Times 3..n depend on the two previous time points.
    if (n_time >= 3) for (t in 3:n_time) {
      G <- matrix(
        c(
          S2_centered[t - 1, t - 1], S2_centered[t - 1, t - 2],
          S2_centered[t - 2, t - 1], S2_centered[t - 2, t - 2]
        ),
        nrow = 2,
        byrow = TRUE
      )
      rhs <- c(S2_centered[t, t - 1], S2_centered[t, t - 2])

      beta <- tryCatch(
        solve(G, rhs),
        error = function(e) solve(G + diag(eps, 2), rhs)
      )

      phi[, t] <- as.numeric(beta)
      sigma2[t] <- max(S2_centered[t, t] - sum(phi[, t] * rhs), eps)
    }

    sigma <- sqrt(sigma2)
  }
  
  list(mu = mu, phi = phi, sigma = sigma, tau = tau)
}
