#' Bayesian information criterion for fitted INAD models
#'
#' Computes BIC using the fitted log likelihood and a parameter count that
#' respects structural zeros and identifiability constraints.
#'
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#' @param n_subjects Number of subjects, typically \code{nrow(y)}. If
#'   \code{NULL}, inferred from \code{fit$settings$n_subjects} or legacy
#'   \code{length(fit$settings$blocks)} when available (with a warning).
#'
#' @return A numeric scalar BIC value.
#'
#' @details
#' The BIC is computed as:
#' \deqn{BIC = -2 \times \ell + k \times \log(N)}
#' where \eqn{\ell} is the log-likelihood, \eqn{k} is the number of free parameters,
#' and \eqn{N} is the number of subjects.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 40,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 20)
#' bic_inad(fit, n_subjects = nrow(y))
#' @export
bic_inad <- function(fit, n_subjects = NULL) {
    if (is.null(n_subjects)) {
        n_subjects <- fit$settings$n_subjects
        if (is.null(n_subjects) && !is.null(fit$settings$blocks)) {
            warning(
                "n_subjects inferred from legacy fit$settings$blocks; please refit or pass n_subjects explicitly.",
                call. = FALSE
            )
            n_subjects <- length(fit$settings$blocks)
        }
    }
    if (is.null(n_subjects) || length(n_subjects) != 1 || !is.finite(n_subjects) || n_subjects <= 0) {
        stop("n_subjects must be a positive finite scalar")
    }

    k <- .count_params_inad_fit(fit)
    -2 * fit$log_l + k * log(n_subjects)
}

#' Akaike information criterion for fitted INAD models
#'
#' Computes AIC using the fitted log likelihood and a parameter count that
#' respects structural zeros and identifiability constraints.
#'
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#'
#' @return A numeric scalar AIC value.
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2 \times \ell + 2k}
#' where \eqn{\ell} is the log-likelihood and \eqn{k} is the number of free
#' parameters.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 40,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 20)
#' aic_inad(fit)
#' @export
aic_inad <- function(fit) {
    k <- .count_params_inad_fit(fit)
    -2 * fit$log_l + 2 * k
}

#' @keywords internal
.count_params_inad_fit <- function(fit) {
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")
    if (is.null(fit$settings$innovation)) stop("fit$settings$innovation is missing")

    ord <- fit$settings$order
    innovation <- fit$settings$innovation

    N <- length(fit$theta)
    if (N < 1) stop("fit$theta must be nonempty")

    k <- 0
    k <- k + N

    if (ord == 1) k <- k + (N - 1)
    if (ord == 2) k <- k + (2 * N - 3)

    if (!is.null(fit$tau)) {
        B <- length(fit$tau)
        if (B >= 1) k <- k + (B - 1)
    }

    if (!is.null(fit$nb_inno_size) && innovation == "nbinom") {
        if (length(fit$nb_inno_size) == 1L) k <- k + 1 else k <- k + N
    }

    k
}
