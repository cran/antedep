#' Bell numbers
#'
#' Compute Bell numbers up to order z using a truncated series representation.
#'
#' This function is for internal use. The results are stored in \code{.Bell_num}
#' and used by the Bell distribution functions.
#'
#' @param z nonnegative integer giving the maximum order.
#'
#' @return numeric vector of length \code{z + 1} containing B_0, ..., B_z.
#' @keywords internal
Bz <- function(z) {
    v <- seq(0, z, by = 1L)
    k <- seq(0, 20, by = 1L)
    sapply(
        v,
        function(n) {
            sum(k ^ n / factorial(k)) / exp(1)
        }
    )
}

.BELL_MAX_Z <- 236L
.Bell_num <- Bz(.BELL_MAX_Z)

# Precompute log of Bell numbers for numerical stability
.log_Bell_num <- log(.Bell_num)

#' The Bell distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Bell distribution with parameter \code{theta}.
#'
#' Let \eqn{B_x} denote the xth Bell number. The Bell distribution has
#' probability mass function
#' \deqn{
#' P(X = x) = \theta^x \exp(-\exp(\theta) + 1) \frac{B_x}{x!},
#' }
#' for nonnegative integers \eqn{x} and \eqn{\theta \ge 0}.
#'
#' For \eqn{\theta > 0}, the Bell mean is \eqn{E[X] = \theta e^\theta}.
#' At \eqn{\theta = 0}, the distribution is degenerate at 0.
#'
#' The functions follow the standard naming used in base R:
#' \code{dbell} for the density, \code{pbell} for the distribution function,
#' \code{qbell} for the quantile function and \code{rbell} for random
#' generation.
#'
#' @name Bell
#' @importFrom stats rmultinom
#' @param x vector of nonnegative integers (for \code{dbell} and \code{pbell}).
#' @param p numeric vector of probabilities between 0 and 1 inclusive (for \code{qbell}).
#' @param n number of observations to generate (for \code{rbell}).
#' @param theta scalar nonnegative Bell parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param max_z maximum support value used for approximation in
#'   \code{rbell} and \code{qbell}.
#'
#' @return
#' For \code{dbell}, a numeric vector of probabilities.
#' For \code{pbell}, a numeric vector of cumulative probabilities.
#' For \code{qbell}, an integer vector of quantiles.
#' For \code{rbell}, an integer vector of random values.
#'
#' @examples
#' dbell(0:5, theta = 1)
#' pbell(0:5, theta = 1)
#' qbell(c(0.25, 0.5, 0.9), theta = 1)
#' set.seed(1)
#' rbell(10, theta = 1)
#'
#' @export
dbell <- function(x, theta, log = FALSE) {
    x_int <- as.integer(x)
    if (any(x_int < 0L)) {
        stop("x must be nonnegative")
    }
    if (any(x_int > .BELL_MAX_Z)) {
        stop("x exceeds precomputed Bell numbers; increase .BELL_MAX_Z if needed")
    }
    theta <- .check_bell_theta(theta)

    if (theta == 0) {
        log_p0 <- ifelse(x_int == 0L, 0, -Inf)
        if (log) return(log_p0)
        return(exp(log_p0))
    }

    # Use log-scale computation to avoid overflow for large x
    # log P(X=x) = x*log(theta) - exp(theta) + 1 + log(B_x) - log(x!)
    log_p <- x_int * log(theta) - exp(theta) + 1 + .log_Bell_num[x_int + 1L] - lfactorial(x_int)

    if (log) {
        return(log_p)
    }
    exp(log_p)
}

#' @rdname Bell
#' @export
pbell <- function(x, theta) {
    theta <- .check_bell_theta(theta)
    x_int <- as.integer(x)
    if (any(x_int < 0L)) {
        stop("x must be nonnegative")
    }
    sapply(
        x_int,
        function(z) {
            sum(dbell(0:z, theta))
        }
    )
}

#' @rdname Bell
#' @export
rbell <- function(n, theta, max_z = 100L) {
    theta <- .check_bell_theta(theta)
    if (length(n) != 1L || n < 0L) {
        stop("n must be a nonnegative integer")
    }
    if (max_z > .BELL_MAX_Z) {
        stop("max_z exceeds precomputed Bell numbers; increase .BELL_MAX_Z if needed")
    }
    z_vals <- 0:max_z
    p <- dbell(z_vals, theta)
    p <- p / sum(p)
    if (n == 0L) {
        return(integer(0L))
    }
    sampling <- rmultinom(n, size = 1L, prob = p)
    out <- integer(n)
    for (i in seq_len(n)) {
        out[i] <- which(sampling[, i] == 1L) - 1L
    }
    out
}

#' @rdname Bell
#' @export
qbell <- function(p, theta, max_z = 100L) {
    theta <- .check_bell_theta(theta)
    if (any(p < 0 | p > 1)) {
        stop("p must be between 0 and 1 inclusive")
    }
    if (max_z > .BELL_MAX_Z) {
        stop("max_z exceeds precomputed Bell numbers; increase .BELL_MAX_Z if needed")
    }
    z_vals <- 0:max_z
    pmf <- dbell(z_vals, theta)
    pmf <- pmf / sum(pmf)
    cdf <- cumsum(pmf)
    sapply(
        p,
        function(prob) {
            which(cdf >= prob)[1L] - 1L
        }
    )
}

#' @keywords internal
.check_bell_theta <- function(theta) {
    if (length(theta) != 1L || !is.finite(theta) || theta < 0) {
        stop("theta must be a finite scalar >= 0")
    }
    as.numeric(theta)
}
