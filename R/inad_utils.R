# Internal utility functions for INAD models
# These functions are shared across loglik_inad.R and other INAD functions

#' Log-sum-exp for numerical stability
#' @keywords internal
.log_sum_exp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(-Inf)
    m + log(sum(exp(x - m)))
}

#' Convert Bell parameter to Bell mean
#' @keywords internal
.bell_mean_from_theta <- function(theta) {
    theta * exp(theta)
}

#' Convert Bell mean to Bell parameter
#' @keywords internal
.bell_theta_from_mean <- function(mu) {
    mu <- as.numeric(mu)
    out <- rep(NA_real_, length(mu))
    ok <- is.finite(mu) & mu > 0
    if (!any(ok)) return(out)

    mu_ok <- mu[ok]
    uniq <- unique(mu_ok)
    th_uniq <- vapply(uniq, solve_theta_bell, numeric(1))
    out[ok] <- th_uniq[match(mu_ok, uniq)]
    out
}

#' Effective innovation distribution parameter under block effects
#'
#' For Poisson and negative binomial innovations, block effects are additive on
#' the mean. For Bell innovations, block effects are additive on the innovation
#' mean `theta * exp(theta)`, then mapped back to the Bell parameter.
#'
#' @keywords internal
.inad_effective_innovation_param <- function(theta, tau, innovation) {
    if (innovation == "bell") {
        mu <- .bell_mean_from_theta(theta) + tau
        if (any(!is.finite(mu) | mu <= 0)) return(rep(NA_real_, length(mu)))
        return(.bell_theta_from_mean(mu))
    }
    theta + tau
}

#' Effective innovation mean under block effects
#' @keywords internal
.inad_effective_innovation_mean <- function(theta, tau, innovation) {
    if (innovation == "bell") {
        return(.bell_mean_from_theta(theta) + tau)
    }
    theta + tau
}

#' Thinning probability vector
#'
#' Computes the probability mass for thinned counts.
#'
#' @param k_vals Integer vector of possible thinned values.
#' @param yprev Previous count (size parameter for thinning).
#' @param a Thinning parameter alpha.
#' @param thinning One of "binom", "pois", "nbinom".
#'
#' @return Numeric vector of probabilities.
#' @keywords internal
.thin_vec <- function(k_vals, yprev, a, thinning) {
    if (!is.finite(a) || a < 0) return(rep(0, length(k_vals)))
    if (!is.finite(yprev) || yprev < 0) return(rep(0, length(k_vals)))

    if (thinning == "binom") {
        if (a > 1) return(rep(0, length(k_vals)))
        return(stats::dbinom(k_vals, size = yprev, prob = a))
    }

    if (thinning == "pois") {
        return(stats::dpois(k_vals, lambda = a * yprev))
    }

    # nbinom thinning: alpha * Y | Y ~ NegBin(Y, 1/(1+alpha))
    # See Supplement A.1 of the paper
    if (yprev == 0) {
        return(as.numeric(k_vals == 0))
    }
    p_nb <- 1 / (1 + a)
    stats::dnbinom(k_vals, size = yprev, prob = p_nb)
}

#' Innovation probability vector
#'
#' Computes the probability mass for innovation counts.
#'
#' @param u_vals Integer vector of possible innovation values.
#' @param lam Innovation parameter used by the PMF. This is the mean for
#'   `"pois"` and `"nbinom"`, and the Bell parameter for `"bell"`.
#' @param innovation One of "pois", "bell", "nbinom".
#' @param sz Size parameter for negative binomial innovation (ignored otherwise).
#'
#' @return Numeric vector of probabilities.
#' @keywords internal
.innov_vec <- function(u_vals, lam, innovation, sz) {
    if (!is.finite(lam) || lam <= 0) return(rep(0, length(u_vals)))

    if (innovation == "pois") return(stats::dpois(u_vals, lambda = lam))
    if (innovation == "bell") return(dbell(u_vals, theta = lam))
    # nbinom: size = sz (dispersion), mu = lam (mean)
    # See Supplement B.3 of the paper
    stats::dnbinom(u_vals, size = sz, mu = lam)
}

#' Unpack alpha parameters
#'
#' Converts alpha from various input formats to a standardized list.
#'
#' @param order Model order (0, 1, or 2).
#' @param alpha Alpha parameter in various formats.
#' @param N Number of time points.
#'
#' @return List with elements a1 and a2.
#' @keywords internal
.unpack_alpha <- function(order, alpha, N) {
    if (order == 0) return(list(a1 = NULL, a2 = NULL))

    if (order == 1) {
        a <- as.numeric(alpha)
        if (length(a) == 1L) a <- rep(a, N)
        if (length(a) != N) stop("alpha must be length 1 or ncol(y) for order=1")
        return(list(a1 = a, a2 = NULL))
    }

    if (is.matrix(alpha) && ncol(alpha) >= 2) {
        if (nrow(alpha) == 1L) alpha <- alpha[rep.int(1L, N), , drop = FALSE]
        if (nrow(alpha) != N) stop("for order=2, alpha matrix must have nrow 1 or ncol(y)")
        return(list(a1 = as.numeric(alpha[, 1]), a2 = as.numeric(alpha[, 2])))
    }

    if (is.list(alpha) && length(alpha) >= 2) {
        a1 <- as.numeric(alpha[[1]])
        a2 <- as.numeric(alpha[[2]])
        if (length(a1) == 1L) a1 <- rep(a1, N)
        if (length(a2) == 1L) a2 <- rep(a2, N)
        if (length(a1) != N || length(a2) != N) stop("alpha[[1]] and alpha[[2]] must be length 1 or ncol(y)")
        return(list(a1 = a1, a2 = a2))
    }

    stop("for order=2, alpha must be a matrix with 2 columns or a list(alpha1, alpha2)")
}
