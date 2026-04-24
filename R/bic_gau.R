# File: R/bic_gau.R

#' Bayesian information criterion for fitted Gaussian AD models
#'
#' Computes BIC using the fitted log likelihood and a parameter count that
#' respects identifiability constraints for the Gaussian antedependence
#' parameters.
#'
#' @param fit A fitted model object returned by \code{\link{fit_gau}}.
#' @param n_subjects Number of subjects, typically \code{nrow(y)}. If
#'   \code{NULL}, inferred from \code{fit$settings$n_subjects}.
#'
#' @return A numeric scalar BIC value.
#'
#' @details
#' The BIC is computed as:
#' \deqn{BIC = -2 \times \ell + k \times \log(N)}
#' where \eqn{\ell} is the log-likelihood, \eqn{k} is the number of free parameters,
#' and \eqn{N} is the number of subjects.
#'
#' This function applies to Gaussian AD fits from \code{\link{fit_gau}}.
#' For categorical and INAD models, use \code{\link{bic_cat}} and
#' \code{\link{bic_inad}}.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 30, n_time = 5, order = 1, phi = 0.3)
#' fit <- fit_gau(y, order = 1)
#' bic_gau(fit, n_subjects = nrow(y))
#' @export
bic_gau <- function(fit, n_subjects = NULL) {
    if (is.null(n_subjects)) n_subjects <- fit$settings$n_subjects
    if (is.null(n_subjects) || length(n_subjects) != 1 || !is.finite(n_subjects) || n_subjects <= 0) {
        stop("n_subjects must be a positive finite scalar")
    }

    k <- .count_params_gau_fit(fit)
    -2 * fit$log_l + k * log(n_subjects)
}

#' Akaike information criterion for fitted Gaussian AD models
#'
#' Computes AIC using the fitted log likelihood and a parameter count that
#' respects identifiability constraints for the Gaussian antedependence
#' parameters.
#'
#' @param fit A fitted model object returned by \code{\link{fit_gau}}.
#'
#' @return A numeric scalar AIC value.
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2 \times \ell + 2k}
#' where \eqn{\ell} is the log-likelihood and \eqn{k} is the number of free
#' parameters.
#'
#' This function applies to Gaussian AD fits from \code{\link{fit_gau}}.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 30, n_time = 5, order = 1, phi = 0.3)
#' fit <- fit_gau(y, order = 1)
#' aic_gau(fit)
#' @export
aic_gau <- function(fit) {
    k <- .count_params_gau_fit(fit)
    -2 * fit$log_l + 2 * k
}

#' @keywords internal
.count_params_gau_fit <- function(fit) {
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")

    ord <- fit$settings$order
    if (!(ord %in% c(0, 1, 2))) stop("fit$settings$order must be 0, 1, or 2")

    N <- NULL
    if (!is.null(fit$sigma)) N <- length(as.numeric(fit$sigma))
    if (is.null(N) && !is.null(fit$mu)) N <- length(as.numeric(fit$mu))
    if (is.null(N) || N < 1) stop("Cannot infer n_time from fit; need fit$sigma or fit$mu")

    estimate_mu <- fit$settings$estimate_mu
    if (is.null(estimate_mu)) estimate_mu <- !is.null(fit$mu)
    estimate_mu <- isTRUE(estimate_mu)

    k <- 0

    if (estimate_mu) k <- k + N

    k <- k + N

    if (ord == 1) k <- k + (N - 1)
    if (ord == 2) k <- k + (2 * N - 3)

    if (!is.null(fit$tau)) {
        B <- length(as.numeric(fit$tau))
        if (B >= 1) k <- k + (B - 1)
    }

    k
}
