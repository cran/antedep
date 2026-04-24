#' Simulate Gaussian antedependence series
#'
#' Generate longitudinal continuous data from a Gaussian antedependence (AD)
#' model of order 0, 1, or 2 using a conditional regression on predecessors.
#'
#' For `order = 0`, each time point is generated independently as
#' `Y[, t] = mu[t] + tau[block] + eps`, with `eps ~ N(0, sigma[t]^2)`.
#'
#' For `order = 1`, for `t >= 2`:
#' `Y[, t] = m_t + phi[t] * (Y[, t - 1] - m_{t-1}) + eps_t`,
#' where `m_t = mu[t] + tau[block]` and `eps_t ~ N(0, sigma[t]^2)`.
#'
#' For `order = 2`, for `t >= 3`:
#' `Y[, t] = m_t + phi[1, t] * (Y[, t - 1] - m_{t-1}) +
#'                 phi[2, t] * (Y[, t - 2] - m_{t-2}) + eps_t`.
#'
#' If `blocks` is provided, each subject s belongs to a block and receives a
#' mean shift `tau[blocks[s]]`. `tau[1]` is forced to 0.
#'
#' @param n_subjects number of subjects
#' @param n_time number of time points
#' @param order antedependence order, 0, 1 or 2
#' @param mu mean parameter; `NULL`, scalar, or length `n_time`
#' @param phi dependence parameter; ignored when `order = 0`.
#'   For `order = 1`, `NULL`, scalar, or length `n_time`.
#'   For `order = 2`, `NULL` or a 2 by `n_time` matrix.
#' @param sigma innovation standard deviation; `NULL`, scalar, or length `n_time`
#' @param blocks integer vector of length `n_subjects` indicating block membership
#'   for each subject; if `NULL`, no block effect is applied
#' @param tau group effect vector indexed by block; `tau[1]` is forced to 0.
#'   If scalar x, it is expanded to `c(0, x, ..., x)` with length equal to the
#'   number of blocks
#' @param seed optional random seed for reproducibility
#'
#' @return numeric matrix with dimension `n_subjects` by `n_time`
#'
#' @examples
#' y <- simulate_gau(
#'   n_subjects = 20,
#'   n_time = 6,
#'   order = 1,
#'   phi = 0.4,
#'   seed = 42
#' )
#' dim(y)
#' @export
simulate_gau <- function(n_subjects,
                        n_time,
                        order = 1L,
                        mu = NULL,
                        phi = NULL,
                        sigma = NULL,
                        blocks = NULL,
                        tau = 0,
                        seed = NULL) {

    if (!is.null(seed)) set.seed(seed)

    if (!order %in% c(0L, 1L, 2L)) {
        stop("order must be 0, 1 or 2")
    }

    if (length(n_subjects) != 1L || n_subjects <= 0) {
        stop("n_subjects must be a positive integer")
    }
    if (length(n_time) != 1L || n_time <= 0) {
        stop("n_time must be a positive integer")
    }

    n_subjects <- as.integer(n_subjects)
    n_time <- as.integer(n_time)
    order <- as.integer(order)

    if (!is.null(blocks)) {
        if (length(blocks) != n_subjects) {
            stop("blocks must have length n_subjects")
        }
        if (any(is.na(blocks))) {
            stop("blocks must have no missing values")
        }
        blocks <- as.integer(factor(blocks))
        B <- max(blocks)
        if (length(tau) == 1L) {
            tau <- c(0, rep(as.numeric(tau), max(B - 1L, 0L)))
        } else {
            tau <- as.numeric(tau)
            if (length(tau) != B) {
                stop("when blocks is provided, tau must be scalar or have length equal to max(blocks)")
            }
            tau[1L] <- 0
        }
    } else {
        B <- 1L
        blocks <- rep(1L, n_subjects)
        tau <- 0
    }

    if (is.null(mu)) {
        mu <- rep(0, n_time)
    } else if (length(mu) == 1L) {
        mu <- rep(as.numeric(mu), n_time)
    } else {
        mu <- as.numeric(mu)
    }
    if (length(mu) != n_time) {
        stop("mu must be NULL, scalar, or length n_time")
    }

    if (is.null(sigma)) {
        sigma <- rep(1, n_time)
    } else if (length(sigma) == 1L) {
        sigma <- rep(as.numeric(sigma), n_time)
    } else {
        sigma <- as.numeric(sigma)
    }
    if (length(sigma) != n_time) {
        stop("sigma must be NULL, scalar, or length n_time")
    }
    if (any(!is.finite(sigma)) || any(sigma <= 0)) {
        stop("sigma must be positive and finite")
    }

    if (order == 0L) {
        phi_use <- 0
    } else if (order == 1L) {
        if (is.null(phi)) {
            phi_use <- rep(0.5, n_time)
            phi_use[1L] <- 0
        } else if (length(phi) == 1L) {
            phi_use <- rep(as.numeric(phi), n_time)
        } else {
            phi_use <- as.numeric(phi)
        }
        if (length(phi_use) != n_time) {
            stop("for order 1, phi must be NULL, scalar, or length n_time")
        }
        if (any(!is.finite(phi_use))) {
            stop("phi must be finite")
        }
    } else {
        if (is.null(phi)) {
            phi_use <- matrix(0, nrow = 2L, ncol = n_time)
            if (n_time >= 2L) {
                phi_use[1L, 2:n_time] <- 0.5
            }
            if (n_time >= 3L) {
                phi_use[2L, 3:n_time] <- 0.25
            }
        } else {
            if (!is.matrix(phi) || nrow(phi) != 2L || ncol(phi) != n_time) {
                stop("for order 2, phi must be NULL or a 2 by n_time matrix")
            }
            phi_use <- phi
        }
        if (any(!is.finite(phi_use))) {
            stop("phi must be finite")
        }
    }

    y <- matrix(NA_real_, nrow = n_subjects, ncol = n_time)

    m_at <- function(t_index) {
        mu[t_index] + tau[blocks]
    }

    if (n_time >= 1L) {
        y[, 1L] <- m_at(1L) + stats::rnorm(n_subjects, mean = 0, sd = sigma[1L])
    }

    if (order == 0L) {
        if (n_time >= 2L) {
            for (t in 2L:n_time) {
                y[, t] <- m_at(t) + stats::rnorm(n_subjects, mean = 0, sd = sigma[t])
            }
        }
        return(y)
    }

    if (order == 1L) {
        if (n_time >= 2L) {
            for (t in 2L:n_time) {
                mt <- m_at(t)
                mtm1 <- m_at(t - 1L)
                y[, t] <- mt +
                    phi_use[t] * (y[, t - 1L] - mtm1) +
                    stats::rnorm(n_subjects, mean = 0, sd = sigma[t])
            }
        }
        return(y)
    }

    if (n_time >= 2L) {
        t <- 2L
        mt <- m_at(t)
        mtm1 <- m_at(t - 1L)
        y[, t] <- mt +
            phi_use[1L, t] * (y[, t - 1L] - mtm1) +
            stats::rnorm(n_subjects, mean = 0, sd = sigma[t])
    }

    if (n_time >= 3L) {
        for (t in 3L:n_time) {
            mt <- m_at(t)
            mtm1 <- m_at(t - 1L)
            mtm2 <- m_at(t - 2L)
            y[, t] <- mt +
                phi_use[1L, t] * (y[, t - 1L] - mtm1) +
                phi_use[2L, t] * (y[, t - 2L] - mtm2) +
                stats::rnorm(n_subjects, mean = 0, sd = sigma[t])
        }
    }

    y
}
