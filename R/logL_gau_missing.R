#' Compute observed-data log-likelihood for AD with missing values
#'
#' Uses multivariate normal marginalization to compute the likelihood
#' of observed values, marginalizing over missing values.
#'
#' @param y Data matrix (n_subjects x n_time), may contain NA
#' @param order Antedependence order (0, 1, or 2)
#' @param mu Mean vector (length n_time)
#' @param phi Dependence parameter(s)
#' @param sigma Innovation standard deviations (length n_time)
#' @param blocks Block membership vector (optional)
#' @param tau Block effects, first element constrained to zero
#'
#' @return Observed-data log-likelihood (scalar)
#'
#' @keywords internal
.logL_gau_missing <- function(y, order, mu, phi, sigma, blocks = NULL, tau = 0) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Build full covariance matrix from AD parameters
  Sigma <- .build_gau_covariance(order, phi, sigma, n_time)
  
  # Adjust mean for blocks if present
  if (!is.null(blocks)) {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
    tau_matrix <- matrix(tau[blocks], nrow = n_subjects, ncol = n_time)
    mu_matrix <- mu_matrix + tau_matrix
  } else {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
  }
  
  # Compute log-likelihood by marginalizing over missing values
  log_lik <- 0
  
  for (s in 1:n_subjects) {
    y_s <- y[s, ]
    mu_s <- mu_matrix[s, ]
    
    obs_idx <- which(!is.na(y_s))
    
    if (length(obs_idx) == 0) {
      # All missing - skip
      next
    }
    
    # Extract observed values and parameters
    y_obs <- y_s[obs_idx]
    mu_obs <- mu_s[obs_idx]
    Sigma_obs <- Sigma[obs_idx, obs_idx, drop = FALSE]
    
    # Compute MVN density on observed values
    log_lik_s <- tryCatch({
      # Use Cholesky decomposition for numerical stability
      L <- chol(Sigma_obs)
      z <- backsolve(L, y_obs - mu_obs, transpose = TRUE)
      log_det <- 2 * sum(log(diag(L)))
      -0.5 * (length(obs_idx) * log(2 * pi) + log_det + sum(z^2))
    }, error = function(e) {
      NA_real_
    })

    if (!is.finite(log_lik_s)) return(-Inf)
    log_lik <- log_lik + log_lik_s
  }
  
  return(log_lik)
}

#' Build AD covariance matrix from parameters
#'
#' @param order Antedependence order
#' @param phi Dependence parameters
#' @param sigma Innovation standard deviations
#' @param n_time Number of time points
#' @return Covariance matrix (n_time x n_time)
#' @keywords internal
.build_gau_covariance <- function(order, phi, sigma, n_time) {
  Sigma <- matrix(0, nrow = n_time, ncol = n_time)

  if (length(sigma) != n_time || any(!is.finite(sigma)) || any(sigma <= 0)) {
    stop("sigma must be positive and have length n_time")
  }

  if (order == 0) {
    diag(Sigma) <- sigma^2
    return(Sigma)
  }

  if (order == 1) {
    if (length(phi) == n_time) {
      phi1 <- as.numeric(phi)
    } else if (length(phi) == n_time - 1) {
      phi1 <- c(0, as.numeric(phi))
    } else {
      stop("phi must have length n_time or n_time - 1 for order 1")
    }

    V <- numeric(n_time)
    V[1] <- sigma[1]^2
    for (t in 2:n_time) {
      V[t] <- phi1[t]^2 * V[t - 1] + sigma[t]^2
    }

    for (i in 1:n_time) {
      Sigma[i, i] <- V[i]
      if (i < n_time) {
        for (j in (i + 1):n_time) {
          prod_phi <- prod(phi1[(i + 1):j])
          Sigma[i, j] <- prod_phi * V[i]
          Sigma[j, i] <- Sigma[i, j]
        }
      }
    }
    return(Sigma)
  }

  if (order == 2) {
    if (is.matrix(phi)) {
      if (nrow(phi) == 2 && ncol(phi) == n_time) {
        phi2 <- phi
      } else if (nrow(phi) == n_time - 2 && ncol(phi) == 2) {
        phi2 <- matrix(0, nrow = 2, ncol = n_time)
        if (n_time >= 3) {
          phi2[1, 3:n_time] <- phi[, 1]
          phi2[2, 3:n_time] <- phi[, 2]
        }
      } else {
        stop("phi must be 2 x n_time or (n_time - 2) x 2 matrix for order 2")
      }
    } else {
      phi_vec <- as.numeric(phi)
      if (length(phi_vec) == 2 * n_time) {
        phi2 <- matrix(phi_vec, nrow = 2, byrow = TRUE)
      } else if (length(phi_vec) == 2 * (n_time - 2)) {
        phi2 <- matrix(0, nrow = 2, ncol = n_time)
        if (n_time >= 3) {
          phi_compact <- matrix(phi_vec, ncol = 2, byrow = TRUE)
          phi2[1, 3:n_time] <- phi_compact[, 1]
          phi2[2, 3:n_time] <- phi_compact[, 2]
        }
      } else {
        stop("phi must be compatible with order 2")
      }
    }

    Sigma[1, 1] <- sigma[1]^2

    if (n_time >= 2) {
      Sigma[2, 1] <- phi2[1, 2] * Sigma[1, 1]
      Sigma[1, 2] <- Sigma[2, 1]
      Sigma[2, 2] <- phi2[1, 2]^2 * Sigma[1, 1] + sigma[2]^2
    }

    if (n_time >= 3) {
      for (t in 3:n_time) {
        for (k in 1:(t - 1)) {
          Sigma[t, k] <- phi2[1, t] * Sigma[t - 1, k] + phi2[2, t] * Sigma[t - 2, k]
          Sigma[k, t] <- Sigma[t, k]
        }
        Sigma[t, t] <- phi2[1, t] * Sigma[t - 1, t] +
          phi2[2, t] * Sigma[t - 2, t] + sigma[t]^2
      }
    }

    return(Sigma)
  }

  stop("order must be 0, 1, or 2")
}
