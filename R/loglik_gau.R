#' Log-likelihood for Gaussian AD models (with missing data support)
#'
#' Computes the log-likelihood for Gaussian antedependence models of order 0, 1, or 2.
#' Supports missing data under MAR assumption via na_action parameter.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns. May contain NA.
#' @param order Antedependence order, one of 0, 1, or 2.
#' @param mu Mean vector (length n_time).
#' @param phi Dependence coefficient(s). For order 1: vector of length n_time-1.
#'   For order 2: matrix with 2 columns or vector of length 2*(n_time-2).
#' @param sigma Innovation standard deviations (length n_time).
#' @param blocks Integer vector of block membership (length n_subjects), or NULL.
#' @param tau Block effects, first element constrained to zero
#' @param na_action How to handle missing values:
#'   \itemize{
#'     \item \code{fail}: Error if any NA is present (default)
#'     \item \code{complete}: Use only complete cases
#'     \item \code{marginalize}: Compute observed-data likelihood
#'   }
#'
#' @return Scalar log-likelihood value.
#'
#' @details
#' For complete data (no NA), all three na_action options give the same result.
#'
#' For missing data:
#' \itemize{
#'   \item marginalize: Uses MVN marginalization to compute P(Y_obs). This is
#'     the correct observed-data likelihood for MAR missing data.
#'   \item complete: Removes subjects with any missing values. May lose information.
#'   \item fail: Stops with error. Useful to ensure no missing data present.
#' }
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 30, n_time = 5, order = 1, phi = 0.3)
#' fit <- fit_gau(y, order = 1)
#' logL_gau(y, order = 1, mu = fit$mu, phi = fit$phi, sigma = fit$sigma)
#'
#' @export
#' @importFrom stats dnorm
logL_gau <- function(y, order, mu = NULL, phi = NULL, sigma = NULL, blocks = NULL, tau = 0,
                    na_action = c("fail", "complete", "marginalize")) {
  na_action <- match.arg(na_action)

  if (!is.matrix(y)) y <- as.matrix(y)

  n_subjects <- nrow(y)
  n_time <- ncol(y)

  if (!order %in% c(0, 1, 2)) {
    stop("order must be 0, 1, or 2", call. = FALSE)
  }

  .expand_vec <- function(x, default, name) {
    if (is.null(x)) {
      out <- rep(default, n_time)
    } else if (length(x) == 1) {
      out <- rep(as.numeric(x), n_time)
    } else {
      out <- as.numeric(x)
    }
    if (length(out) != n_time) {
      stop(paste0(name, " must have length n_time"), call. = FALSE)
    }
    out
  }

  mu <- .expand_vec(mu, 0, "mu")
  sigma <- .expand_vec(sigma, 1, "sigma")

  if (any(!is.finite(mu))) return(-Inf)
  if (any(!is.finite(sigma)) || any(sigma <= 0)) return(-Inf)

  if (order == 0) {
    phi_norm <- numeric(0)
  } else if (order == 1) {
    if (is.null(phi)) {
      phi_norm <- rep(0, n_time)
    } else if (length(phi) == 1) {
      phi_norm <- rep(as.numeric(phi), n_time)
    } else {
      phi_vec <- as.numeric(phi)
      if (length(phi_vec) == n_time) {
        phi_norm <- phi_vec
      } else if (length(phi_vec) == n_time - 1) {
        phi_norm <- c(0, phi_vec)
      } else {
        stop("phi must have length n_time or n_time - 1 for order 1", call. = FALSE)
      }
    }
  } else {
    if (is.null(phi)) {
      phi_norm <- matrix(0, nrow = 2, ncol = n_time)
    } else if (is.matrix(phi)) {
      if (nrow(phi) == 2 && ncol(phi) == n_time) {
        phi_norm <- phi
      } else {
        stop("phi must be a 2 x n_time matrix for order 2", call. = FALSE)
      }
    } else {
      phi_vec <- as.numeric(phi)
      if (length(phi_vec) == 2 * n_time) {
        phi_norm <- matrix(phi_vec, nrow = 2, byrow = TRUE)
      } else if (length(phi_vec) == 2 * (n_time - 2)) {
        phi_norm <- matrix(0, nrow = 2, ncol = n_time)
        if (n_time >= 3) {
          phi_compact <- matrix(phi_vec, ncol = 2, byrow = TRUE)
          phi_norm[1, 3:n_time] <- phi_compact[, 1]
          phi_norm[2, 3:n_time] <- phi_compact[, 2]
        }
      } else {
        stop("phi must be compatible with order 2", call. = FALSE)
      }
    }
  }

  has_missing <- any(is.na(y))

  if (has_missing) {
    if (na_action == "fail") {
      stop("y contains NA values. Use na_action = 'marginalize' or 'complete'",
           call. = FALSE)
    }

    if (na_action == "complete") {
      cc <- .extract_complete_cases(y, blocks, warn = FALSE)
      y <- cc$y
      blocks <- cc$blocks
      n_subjects <- nrow(y)
      has_missing <- FALSE
    }

    if (na_action == "marginalize" && has_missing) {
      return(.logL_gau_missing(y, order, mu, phi_norm, sigma, blocks, tau))
    }
  }

  if (!is.null(blocks)) {
    if (length(blocks) != n_subjects) stop("blocks must have length n_subjects", call. = FALSE)
    n_blocks <- length(unique(blocks))
    if (length(tau) == 1 && n_blocks > 1) {
      tau <- c(0, rep(as.numeric(tau), n_blocks - 1))
    }
    if (length(tau) != n_blocks) stop("tau must have length equal to number of blocks", call. = FALSE)
    if (abs(tau[1]) > 1e-10) {
      warning("tau[1] should be 0 for identifiability. Setting tau[1] = 0.")
      tau[1] <- 0
    }
  } else {
    tau <- 0
  }

  log_lik <- 0

  for (s in 1:n_subjects) {
    y_s <- y[s, ]

    if (!is.null(blocks)) {
      mu_s <- mu + tau[blocks[s]]
    } else {
      mu_s <- mu
    }

    log_lik <- log_lik + dnorm(y_s[1], mean = mu_s[1], sd = sigma[1], log = TRUE)

    if (order == 0) {
      for (t in 2:n_time) {
        log_lik <- log_lik + dnorm(y_s[t], mean = mu_s[t], sd = sigma[t], log = TRUE)
      }

    } else if (order == 1) {
      for (t in 2:n_time) {
        cond_mean <- mu_s[t] + phi_norm[t] * (y_s[t - 1] - mu_s[t - 1])
        log_lik <- log_lik + dnorm(y_s[t], mean = cond_mean, sd = sigma[t], log = TRUE)
      }

    } else if (order == 2) {
      if (n_time >= 2) {
        cond_mean_2 <- mu_s[2] + phi_norm[1, 2] * (y_s[1] - mu_s[1])
        log_lik <- log_lik + dnorm(y_s[2], mean = cond_mean_2, sd = sigma[2], log = TRUE)
      }

      for (t in 3:n_time) {
        cond_mean <- mu_s[t] +
          phi_norm[1, t] * (y_s[t - 1] - mu_s[t - 1]) +
          phi_norm[2, t] * (y_s[t - 2] - mu_s[t - 2])
        log_lik <- log_lik + dnorm(y_s[t], mean = cond_mean, sd = sigma[t], log = TRUE)
      }
    }
  }

  log_lik
}
