# File: R/gau_utils.R
# Shared utility functions for Gaussian antedependence models

#' Compute residual sum of squares from AD regression
#'
#' For time point i, regress Y_i on its p_i predecessors (Y_(i-1), ..., Y_(i-p_i))
#' and return the residual sum of squares.
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param i Time index (1-based).
#' @param p Number of predecessors to include in regression (can be 0).
#' @param mu Optional mean vector. If provided, data is centered by mu before regression.
#'   If NULL, regression includes an intercept.
#' @param include_intercept Logical. If TRUE and mu is NULL, include intercept in regression.
#'
#' @return Residual sum of squares (scalar).
#'
#' @details
#' When mu is provided, the regression is:
#'   \code{(Y_i - mu_i) ~ (Y_(i-1) - mu_(i-1)) + ... + (Y_(i-p) - mu_(i-p))}
#'   with no intercept.
#'
#' When mu is NULL and include_intercept is TRUE:
#'   \code{Y_i ~ 1 + Y_(i-1) + ... + Y_(i-p)}
#'
#' For i = 1 or p = 0, the RSS is simply the sum of squared deviations from the mean
#' (or from mu_1 if provided).
#'
#' @keywords internal
.rss_gau <- function(y, i, p, mu = NULL, include_intercept = TRUE) {
    n <- nrow(y)

    if (i < 1 || i > ncol(y)) stop("i out of range")

    # Effective order for this time point
    p_i <- min(p, i - 1L)

    if (!is.null(mu)) {
        # Center data by provided mean
        response <- y[, i] - mu[i]

        if (p_i == 0L) {
            # No predecessors - RSS is sum of squared centered values
            return(sum(response^2))
        }

        # Build design matrix (centered predecessors, no intercept)
        X <- matrix(NA_real_, nrow = n, ncol = p_i)
        for (k in 1:p_i) {
            X[, k] <- y[, i - k] - mu[i - k]
        }

        # Regression without intercept
        fit <- stats::lm.fit(X, response)
        rss <- sum(fit$residuals^2)

    } else {
        # No centering - use intercept
        response <- y[, i]

        if (p_i == 0L) {
            # No predecessors - RSS from mean
            if (include_intercept) {
                return(sum((response - mean(response))^2))
            } else {
                return(sum(response^2))
            }
        }

        # Build design matrix
        if (include_intercept) {
            X <- cbind(1, matrix(NA_real_, nrow = n, ncol = p_i))
            for (k in 1:p_i) {
                X[, k + 1] <- y[, i - k]
            }
        } else {
            X <- matrix(NA_real_, nrow = n, ncol = p_i)
            for (k in 1:p_i) {
                X[, k] <- y[, i - k]
            }
        }

        fit <- stats::lm.fit(X, response)
        rss <- sum(fit$residuals^2)
    }

    rss
}


#' Compute RSS vector for all time points under AD(p)
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param p Antedependence order.
#' @param mu Optional mean vector for centering.
#'
#' @return Numeric vector of length n_time with RSS for each time point.
#'
#' @keywords internal
.rss_vector_gau <- function(y, p, mu = NULL) {
    N <- ncol(y)
    rss <- numeric(N)

    for (i in 1:N) {
        rss[i] <- .rss_gau(y, i, p, mu = mu)
    }

    rss
}


#' Compute intervenor-adjusted sample partial correlation
#'
#' Computes the sample partial correlation between Y_i and Y_j, adjusted for
#' the intervenors Y_(j+1), ..., Y_(i-1), from the residuals of the fitted
#' mean structure.
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param i Row/column index (larger time index).
#' @param j Row/column index (smaller time index), j < i.
#' @param mu Optional mean vector for centering.
#'

#' @return Sample partial correlation (scalar).
#'
#' @details
#' The intervenor-adjusted partial correlation r_(i,j|(j+1:i-1)) is computed
#' from residuals of regressing Y_i and Y_j on the intervenors.
#'
#' Under AD(p), the partial correlation r_(i,i-k|(i-k+1:i-1)) should be
#' approximately zero for k > p.
#'
#' @keywords internal
.partial_corr_gau <- function(y, i, j, mu = NULL) {
    if (j >= i) stop("j must be less than i")

    n <- nrow(y)

    # Center data if mu provided
    if (!is.null(mu)) {
        y_i <- y[, i] - mu[i]
        y_j <- y[, j] - mu[j]
    } else {
        y_i <- y[, i] - mean(y[, i])
        y_j <- y[, j] - mean(y[, j])
    }

    # If adjacent (i = j + 1), no intervenors - return ordinary correlation
    if (i == j + 1L) {
        return(stats::cor(y_i, y_j))
    }

    # Build intervenor matrix: Y_{j+1}, ..., Y_{i-1}
    n_intervenors <- i - j - 1L
    Z <- matrix(NA_real_, nrow = n, ncol = n_intervenors)
    for (k in 1:n_intervenors) {
        idx <- j + k
        if (!is.null(mu)) {
            Z[, k] <- y[, idx] - mu[idx]
        } else {
            Z[, k] <- y[, idx] - mean(y[, idx])
        }
    }

    # Regress y_i and y_j on intervenors.
    # When mu is provided the data is already centered, so no intercept is needed.
    # When mu is NULL the data is centered by column means, so an intercept is included.
    if (!is.null(mu)) {
        fit_i <- stats::lm.fit(Z, y_i)
        fit_j <- stats::lm.fit(Z, y_j)
    } else {
        fit_i <- stats::lm.fit(cbind(1, Z), y_i)
        fit_j <- stats::lm.fit(cbind(1, Z), y_j)
    }

    # Partial correlation is correlation of residuals
    stats::cor(fit_i$residuals, fit_j$residuals)
}


#' Compute matrix of intervenor-adjusted partial correlations
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param mu Optional mean vector for centering.
#'
#' @return Matrix where entry (i,j) for i > j is r_(i,j|(j+1:i-1)).
#'
#' @keywords internal
.partial_corr_matrix_gau <- function(y, mu = NULL) {
    N <- ncol(y)
    R <- matrix(0, nrow = N, ncol = N)

    for (i in 2:N) {
        for (j in 1:(i-1)) {
            R[i, j] <- .partial_corr_gau(y, i, j, mu = mu)
            R[j, i] <- R[i, j]  # Symmetric
        }
    }

    R
}


#' Compute RSS for two-sample case (pooled and separate)
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param blocks Integer vector indicating group membership (1 or 2).
#' @param i Time index.
#' @param p Antedependence order.
#' @param pooled Logical. If TRUE, compute pooled RSS; if FALSE, within-group RSS.
#'
#' @return RSS value.
#'
#' @keywords internal
.rss_two_sample_gau <- function(y, blocks, i, p, pooled = TRUE) {
    g1 <- which(blocks == 1)
    g2 <- which(blocks == 2)

    if (pooled) {
        # Pooled: common intercept across groups
        .rss_gau(y, i, p, mu = NULL, include_intercept = TRUE)
    } else {
        # Separate: sum of within-group RSS
        rss1 <- .rss_gau(y[g1, , drop = FALSE], i, p, mu = NULL, include_intercept = TRUE)
        rss2 <- .rss_gau(y[g2, , drop = FALSE], i, p, mu = NULL, include_intercept = TRUE)
        rss1 + rss2
    }
}


#' Compute innovation variance estimates under AD(p)
#'
#' The MLE of the innovation variance at time i is RSS_i(p) / N.
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param p Antedependence order.
#' @param mu Optional mean vector for centering.
#'
#' @return Numeric vector of innovation variance estimates.
#'
#' @keywords internal
.innovation_var_gau <- function(y, p, mu = NULL) {
    n <- nrow(y)
    rss <- .rss_vector_gau(y, p, mu = mu)
    rss / n
}


#' Compute MLE of mean vector under AD(p)
#'
#' Under AD(p), the MLE of mu is the sample mean vector when there are no
#' covariates.
#'
#' @param y Numeric matrix n_subjects by n_time.
#'
#' @return Numeric vector of mean estimates.
#'
#' @keywords internal
.mle_mean_gau <- function(y) {
    colMeans(y)
}


#' Compute group-specific means for two-sample case
#'
#' @param y Numeric matrix n_subjects by n_time.
#' @param blocks Integer vector indicating group membership.
#'
#' @return List with mu1 and mu2 vectors.
#'
#' @keywords internal
.group_means_gau <- function(y, blocks) {
    g1 <- which(blocks == 1)
    g2 <- which(blocks == 2)

    list(
        mu1 = colMeans(y[g1, , drop = FALSE]),
        mu2 = colMeans(y[g2, , drop = FALSE])
    )
}


#' Modified psi function for test statistic correction
#'
#' Computes the psi function used in the modified likelihood ratio test
#' (Kenward, 1987) for better small-sample approximation.
#'
#' @param x First argument.
#' @param y Second argument.
#'
#' @return Value of psi(x, y).
#'
#' @details
#' psi(x, y) is defined as:
#' \itemize{
#'   \item 0 if x = 0
#'   \item (2y - 1) / (2y(y - 1)) if x = 1
#'   \item 2 * sum from l=1 to x/2 of (y + 2l - 2)^(-1) if x > 0 and even
#'   \item psi(1, y) + psi(x - 1, y + 1) if x > 1 and odd
#' }
#'
#' @keywords internal
.psi_kenward <- function(x, y) {
    if (x == 0) {
        return(0)
    } else if (x == 1) {
        return((2 * y - 1) / (2 * y * (y - 1)))
    } else if (x > 0 && x %% 2 == 0) {
        # Even
        return(2 * sum(1 / (y + 2 * (1:(x/2)) - 2)))
    } else {
        # Odd and > 1
        return(.psi_kenward(1, y) + .psi_kenward(x - 1, y + 1))
    }
}


#' Count parameters for AD model
#'
#' @param order Model order.
#' @param n_time Number of time points.
#' @param n_blocks Number of blocks (1 if no blocks).
#' @param estimate_mu Whether mu is estimated.
#'
#' @return Number of free parameters.
#'
#' @keywords internal
.count_params_gau <- function(order, n_time, n_blocks = 1L, estimate_mu = FALSE) {
    k <- 0L

    # mu parameters
    if (estimate_mu) {
        k <- k + n_time
    }

    # sigma (innovation variance) parameters: n_time
    k <- k + n_time

    # phi parameters
    # Order p: for i = 2, ..., n, we have min(p, i-1) phi parameters
    # Total = sum_{i=2}^{n} min(p, i-1) = (n-1) + (n-2) + ... for p >= n-1
    # For constant order p: (p)(n - (p+1)/2) approximately
    if (order >= 1L) {
        for (i in 2:n_time) {
            k <- k + min(order, i - 1L)
        }
    }

    # tau (block effect) parameters
    if (n_blocks > 1L) {
        k <- k + (n_blocks - 1L)
    }

    k
}
