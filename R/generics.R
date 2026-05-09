# R/generics.R
# S3 generic methods for antedependence model fit objects:
#   logLik, coef, nobs, fitted, residuals
#
# Once logLik methods are registered, stats::AIC() and stats::BIC() work
# automatically for all three fit families without any extra code.

# ===========================================================================
# logLik
# ===========================================================================

#' Log-likelihood for antedependence model fits
#'
#' Extracts the maximised log-likelihood from a fitted antedependence model
#' and returns a \code{"logLik"} object compatible with \code{\link[stats]{AIC}}
#' and \code{\link[stats]{BIC}}.
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, or \code{inad_fit}.
#' @param ... Unused; present for S3 consistency.
#'
#' @return A scalar of class \code{"logLik"} with attributes
#'   \code{df} (number of free parameters) and \code{nobs}
#'   (number of subjects).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' logLik(fit)
#' AIC(fit)
#' BIC(fit)
#' }
#'
#' @name logLik.antedep
NULL

#' @rdname logLik.antedep
#' @export
logLik.gau_fit <- function(object, ...) {
    ll <- object$log_l
    attr(ll, "df")   <- object$n_params
    attr(ll, "nobs") <- object$settings$n_subjects
    class(ll) <- "logLik"
    ll
}

#' @rdname logLik.antedep
#' @export
logLik.cat_fit <- function(object, ...) {
    ll <- object$log_l
    attr(ll, "df")   <- object$n_params
    attr(ll, "nobs") <- object$settings$n_subjects
    class(ll) <- "logLik"
    ll
}

#' @rdname logLik.antedep
#' @export
logLik.inad_fit <- function(object, ...) {
    ll <- object$log_l
    attr(ll, "df")   <- object$n_params
    attr(ll, "nobs") <- object$settings$n_subjects
    class(ll) <- "logLik"
    ll
}

# ===========================================================================
# nobs
# ===========================================================================

#' Number of observations for antedependence model fits
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, or \code{inad_fit}.
#' @param ... Unused.
#'
#' @return Integer scalar: number of subjects.
#'
#' @name nobs.antedep
NULL

#' @rdname nobs.antedep
#' @importFrom stats nobs
#' @export
nobs.gau_fit <- function(object, ...) object$settings$n_subjects

#' @rdname nobs.antedep
#' @export
nobs.cat_fit <- function(object, ...) object$settings$n_subjects

#' @rdname nobs.antedep
#' @export
nobs.inad_fit <- function(object, ...) object$settings$n_subjects

# ===========================================================================
# coef
# ===========================================================================

#' Extract model coefficients from antedependence fits
#'
#' Returns the estimated parameters of a fitted antedependence model as a
#' named list.
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, or \code{inad_fit}.
#' @param ... Unused.
#'
#' @return For \code{gau_fit}: a named list with elements \code{mu},
#'   \code{phi}, \code{sigma}, and (if a block effect was estimated)
#'   \code{tau}.
#'   For \code{cat_fit}: a list with \code{marginal} and \code{transition}
#'   probability arrays.
#'   For \code{inad_fit}: a named list with elements \code{alpha},
#'   \code{theta}, and (if applicable) \code{tau} and \code{nb_inno_size}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' coef(fit)
#' }
#'
#' @name coef.antedep
NULL

#' @rdname coef.antedep
#' @export
coef.gau_fit <- function(object, ...) {
    out <- list()
    if (!is.null(object$mu))    out$mu    <- object$mu
    if (!is.null(object$phi))   out$phi   <- object$phi
    if (!is.null(object$sigma)) out$sigma <- object$sigma
    if (!is.null(object$tau) && any(as.numeric(object$tau) != 0))
        out$tau <- object$tau
    out
}

#' @rdname coef.antedep
#' @export
coef.cat_fit <- function(object, ...) {
    list(marginal = object$marginal, transition = object$transition)
}

#' @rdname coef.antedep
#' @export
coef.inad_fit <- function(object, ...) {
    out <- list()
    if (!is.null(object$alpha))         out$alpha         <- object$alpha
    if (!is.null(object$theta))         out$theta         <- object$theta
    if (!is.null(object$tau) && any(as.numeric(object$tau) != 0))
        out$tau <- object$tau
    if (!is.null(object$nb_inno_size))  out$nb_inno_size  <- object$nb_inno_size
    out
}

# ===========================================================================
# fitted / residuals  (Gaussian only -- meaningful for complete-data fits)
# ===========================================================================

#' Fitted values and residuals for Gaussian AD fits
#'
#' Computes fitted (conditional mean) values and residuals for a
#' complete-data Gaussian antedependence fit.  The fitted value for subject
#' \eqn{s} at time \eqn{t} is the conditional mean
#' \eqn{E[Y_{st} | Y_{s,t-p}, \ldots, Y_{s,t-1}]}.
#'
#' @param object A \code{gau_fit} object (complete-data fit only).
#' @param y The original data matrix used to produce \code{object}
#'   (n_subjects \eqn{\times} n_time).  Required because the fit object
#'   does not store the raw data.
#' @param ... Unused.
#'
#' @return A numeric matrix of the same dimensions as \code{y}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' yhat <- fitted(fit, y)
#' resid <- residuals(fit, y)
#' }
#'
#' @name fitted.gau_fit
NULL

#' @rdname fitted.gau_fit
#' @export
fitted.gau_fit <- function(object, y, ...) {
    if (missing(y)) stop("'y' must be supplied; it is not stored in the fit object.")
    if (!is.matrix(y)) y <- as.matrix(y)
    n <- nrow(y)
    T <- ncol(y)
    ord <- object$settings$order
    mu  <- as.numeric(object$mu)
    yhat <- matrix(NA_real_, nrow = n, ncol = T)

    # block mean shifts
    tau_vec <- rep(0, n)
    if (!is.null(object$settings$blocks) && !is.null(object$tau)) {
        tau_vec <- as.numeric(object$tau)[object$settings$blocks]
    }
    m <- function(t) mu[t] + tau_vec   # length-n mean vector at time t

    yhat[, 1] <- m(1)

    if (ord == 0) {
        for (t in seq_len(T)) yhat[, t] <- m(t)
        return(yhat)
    }

    if (ord == 1) {
        phi <- as.numeric(object$phi)
        for (t in 2:T)
            yhat[, t] <- m(t) + phi[t] * (y[, t - 1] - m(t - 1))
        return(yhat)
    }

    if (ord == 2) {
        phi <- as.matrix(object$phi)   # 2 x T
        yhat[, 2] <- m(2) + phi[1, 2] * (y[, 1] - m(1))
        for (t in 3:T)
            yhat[, t] <- m(t) +
                phi[1, t] * (y[, t - 1] - m(t - 1)) +
                phi[2, t] * (y[, t - 2] - m(t - 2))
        return(yhat)
    }

    yhat
}

#' @rdname fitted.gau_fit
#' @importFrom stats fitted
#' @export
residuals.gau_fit <- function(object, y, ...) {
    if (missing(y)) stop("'y' must be supplied; it is not stored in the fit object.")
    if (!is.matrix(y)) y <- as.matrix(y)
    y - fitted(object, y)
}

# ===========================================================================
# as.ts  (simulate_* output -> ts / mts)
# ===========================================================================

#' Convert antedependence simulation output to a time-series object
#'
#' Converts the panel matrix returned by \code{\link{simulate_gau}},
#' \code{\link{simulate_cat}}, or \code{\link{simulate_inad}} into an
#' \code{mts} / \code{\link[stats]{ts}} object suitable for time-series
#' functions.  Rows of the original matrix (subjects) become the columns of
#' the returned object; columns (time points) become the rows.
#'
#' @param x A matrix of class \code{gau_sim}, \code{cat_sim}, or
#'   \code{inad_sim} as returned by the corresponding \code{simulate_*}
#'   function.
#' @param ... Additional arguments passed to \code{\link[stats]{ts}}
#'   (e.g., \code{start}, \code{frequency}).
#'
#' @return A \code{ts} (\code{mts}) object with one column per subject.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 5, n_time = 10, order = 1, phi = 0.4)
#' y_ts <- as.ts(y)
#' class(y_ts)    # "mts" "ts" "matrix" "array"
#' dim(y_ts)      # 10 x 5
#' }
#'
#' @name as.ts.antedep_sim
NULL

#' @rdname as.ts.antedep_sim
#' @importFrom stats as.ts
#' @export
as.ts.gau_sim <- function(x, ...) stats::ts(t(unclass(x)), ...)

#' @rdname as.ts.antedep_sim
#' @export
as.ts.cat_sim <- function(x, ...) stats::ts(t(unclass(x)), ...)

#' @rdname as.ts.antedep_sim
#' @export
as.ts.inad_sim <- function(x, ...) stats::ts(t(unclass(x)), ...)

# ===========================================================================
# deviance
# ===========================================================================

#' Deviance for antedependence model fits
#'
#' Returns the deviance \eqn{-2 \ell} for a fitted antedependence model.
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, or \code{inad_fit}.
#' @param ... Unused.
#'
#' @return A scalar deviance value.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' deviance(fit)
#' }
#'
#' @name deviance.antedep
NULL

#' @rdname deviance.antedep
#' @importFrom stats deviance
#' @export
deviance.gau_fit <- function(object, ...) -2 * as.numeric(logLik(object))

#' @rdname deviance.antedep
#' @export
deviance.cat_fit <- function(object, ...) -2 * as.numeric(logLik(object))

#' @rdname deviance.antedep
#' @export
deviance.inad_fit <- function(object, ...) -2 * as.numeric(logLik(object))

# ===========================================================================
# confint
# ===========================================================================

#' Confidence intervals for antedependence model fits
#'
#' Computes Wald confidence intervals for the parameters of a fitted
#' antedependence model by delegating to the corresponding
#' \code{ci_gau}, \code{ci_cat}, or \code{ci_inad} function.
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, or \code{inad_fit}.
#' @param parm Unused (all parameters are returned).
#' @param level Confidence level; default \code{0.95}.
#' @param y (\code{inad_fit} only) The original data matrix used for
#'   fitting, required by \code{\link{ci_inad}}.
#' @param ... Passed to the underlying \code{ci_*} function.
#'
#' @return For \code{gau_fit} and \code{cat_fit}: a matrix with columns
#'   \code{lower} and \code{upper} and one row per parameter.
#'   For \code{inad_fit}: the \code{inad_ci} list returned by
#'   \code{\link{ci_inad}}.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' confint(fit)
#' }
#'
#' @name confint.antedep
NULL

#' @rdname confint.antedep
#' @importFrom stats confint
#' @export
confint.gau_fit <- function(object, parm, level = 0.95, ...) {
    ci <- ci_gau(fit = object, level = level)
    rows <- do.call(rbind, Filter(Negate(is.null), list(ci$mu, ci$phi, ci$sigma)))
    if (is.null(rows)) return(matrix(NA_real_, nrow = 0, ncol = 2,
                                     dimnames = list(NULL, c("lower", "upper"))))
    mat <- matrix(c(rows$lower, rows$upper), ncol = 2,
                  dimnames = list(rows$param, c("lower", "upper")))
    mat
}

#' @rdname confint.antedep
#' @export
confint.cat_fit <- function(object, parm, level = 0.95, y = NULL, ...) {
    ci_cat(fit = object, y = y, level = level)
}

#' @rdname confint.antedep
#' @export
confint.inad_fit <- function(object, parm, level = 0.95, y, ...) {
    if (missing(y)) {
        stop(
            "confint() for inad_fit requires the original data matrix 'y'; ",
            "call ci_inad(y, fit, ...) directly.",
            call. = FALSE
        )
    }
    ci_inad(y = y, fit = object, level = level, ...)
}

# ===========================================================================
# vcov
# ===========================================================================

#' Variance-covariance matrix for Gaussian AD fits
#'
#' Returns an approximate variance-covariance matrix for the estimated
#' parameters of a \code{gau_fit} object.  The matrix is block-diagonal
#' (asymptotically exact for the Gaussian AD model given the factored
#' likelihood) with entries derived from the closed-form standard errors
#' used by \code{\link{ci_gau}}.
#'
#' @param object A \code{gau_fit} object.
#' @param ... Unused.
#'
#' @return A named square matrix of parameter variances (diagonal entries)
#'   and covariances (off-diagonal entries, all zero in the current
#'   implementation).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' vcov(fit)
#' }
#'
#' @importFrom stats vcov
#' @export
vcov.gau_fit <- function(object, ...) {
    ci <- ci_gau(fit = object, level = 0.95)
    rows <- do.call(rbind, Filter(Negate(is.null), list(ci$mu, ci$phi, ci$sigma)))
    if (is.null(rows)) return(matrix(NA_real_, 0, 0))
    vars <- rows$se^2
    names(vars) <- rows$param
    diag(vars)
}
