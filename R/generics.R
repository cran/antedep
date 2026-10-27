# R/generics.R
# S3 methods for antedependence model fit objects.

.ad_logLik <- function(object) {
    ll <- object$log_l
    attr(ll, "df") <- object$n_params
    attr(ll, "nobs") <- object$settings$n_subjects
    class(ll) <- "logLik"
    ll
}

.ad_nobs <- function(object) object$settings$n_subjects

.ad_empty_numeric <- function() {
    stats::setNames(numeric(0), character(0))
}

.ad_named_vector <- function(x, prefix, drop_first = FALSE) {
    if (is.null(x)) return(.ad_empty_numeric())

    if (is.matrix(x) || length(dim(x)) > 1L) {
        dims <- dim(x)
        idx <- arrayInd(seq_along(x), .dim = dims)
        vals <- as.numeric(x)
        keep <- rep(TRUE, length(vals))
        if (drop_first) keep <- rowSums(idx > 1L) > 0L
        vals <- vals[keep]
        idx <- idx[keep, , drop = FALSE]
        nms <- apply(idx, 1L, function(z) paste0(prefix, "[", paste(z, collapse = ","), "]"))
        return(stats::setNames(vals, nms))
    }

    vals <- as.numeric(x)
    idx <- seq_along(vals)
    if (drop_first && length(idx) > 0L) {
        vals <- vals[-1L]
        idx <- idx[-1L]
    }
    if (length(vals) == 0L) return(.ad_empty_numeric())
    stats::setNames(vals, paste0(prefix, "[", idx, "]"))
}

.ad_flatten_list <- function(x, prefix = NULL) {
    if (is.null(x)) return(.ad_empty_numeric())
    if (!is.list(x) || is.data.frame(x)) {
        nm <- if (is.null(prefix) || !nzchar(prefix)) "param" else prefix
        return(.ad_named_vector(x, nm))
    }

    nms <- names(x)
    if (is.null(nms)) nms <- paste0("item", seq_along(x))
    parts <- vector("list", length(x))
    for (i in seq_along(x)) {
        child_prefix <- if (is.null(prefix) || !nzchar(prefix)) nms[[i]] else paste(prefix, nms[[i]], sep = ".")
        parts[[i]] <- .ad_flatten_list(x[[i]], child_prefix)
    }
    out <- unlist(parts, use.names = TRUE)
    if (is.null(out)) return(.ad_empty_numeric())
    out
}

.ad_ci_summary <- function(ci) {
    out <- summary(ci)
    if (is.null(out) || !is.data.frame(out) || nrow(out) == 0L) {
        return(data.frame(param = character(0), est = numeric(0),
                          se = numeric(0), lower = numeric(0),
                          upper = numeric(0)))
    }
    if (!"param" %in% names(out) && "parameter" %in% names(out)) out$param <- out$parameter
    if (!"est" %in% names(out) && "estimate" %in% names(out)) out$est <- out$estimate
    out
}

.ad_confint_matrix <- function(ci, parm = NULL) {
    tab <- .ad_ci_summary(ci)
    if (!all(c("param", "lower", "upper") %in% names(tab))) {
        return(matrix(NA_real_, nrow = 0L, ncol = 2L,
                      dimnames = list(NULL, c("lower", "upper"))))
    }

    mat <- as.matrix(tab[, c("lower", "upper"), drop = FALSE])
    rownames(mat) <- tab$param
    colnames(mat) <- c("lower", "upper")

    if (!missing(parm) && !is.null(parm)) {
        mat <- mat[intersect(as.character(parm), rownames(mat)), , drop = FALSE]
    }
    mat
}

.ad_vcov_from_ci <- function(coefs, ci) {
    vars <- stats::setNames(rep(NA_real_, length(coefs)), names(coefs))
    tab <- .ad_ci_summary(ci)
    if (all(c("param", "se") %in% names(tab))) {
        valid <- !is.na(tab$se)
        hit <- intersect(names(vars), tab$param[valid])
        vars[hit] <- tab$se[match(hit, tab$param)]^2
    }
    mat <- diag(vars, nrow = length(vars), ncol = length(vars))
    dimnames(mat) <- list(names(vars), names(vars))
    mat
}

.ad_anova <- function(object, ..., test = "Chisq") {
    dots <- list(...)
    objects <- c(list(object), dots)
    if (length(objects) == 0L) stop("at least one model object is required")
    if (!all(vapply(objects, inherits, logical(1), what = "ad_fit"))) {
        stop("all objects must inherit from 'ad_fit'", call. = FALSE)
    }

    calls <- as.list(match.call(expand.dots = FALSE)$...)
    labels <- c(deparse(match.call()$object), vapply(calls, deparse, character(1)))
    if (length(labels) != length(objects) || any(!nzchar(labels))) {
        labels <- paste0("model", seq_along(objects))
    }

    ll <- vapply(objects, function(x) as.numeric(logLik(x)), numeric(1))
    df <- vapply(objects, function(x) attr(logLik(x), "df"), numeric(1))
    aic <- vapply(objects, stats::AIC, numeric(1))
    bic <- vapply(objects, stats::BIC, numeric(1))

    chi_df <- c(NA_real_, diff(df))
    chisq <- c(NA_real_, 2 * diff(ll))
    pval <- ifelse(is.finite(chisq) & is.finite(chi_df) & chi_df > 0 & chisq >= 0,
                   stats::pchisq(chisq, df = chi_df, lower.tail = FALSE),
                   NA_real_)

    out <- data.frame(
        Df = df,
        logLik = ll,
        AIC = aic,
        BIC = bic,
        Chisq = chisq,
        `Chi Df` = chi_df,
        `Pr(>Chisq)` = pval,
        check.names = FALSE
    )
    rownames(out) <- labels
    structure(out, heading = "Likelihood-ratio comparison of antedependence fits",
              class = c("anova", "data.frame"))
}

#' Log-Likelihood for Antedependence Model Fits
#'
#' Extracts the maximized log-likelihood from a fitted antedependence model
#' and returns a \code{"logLik"} object compatible with \code{\link[stats]{AIC}}
#' and \code{\link[stats]{BIC}}.
#'
#' @param object A fitted model object of class \code{gau_fit},
#'   \code{cat_fit}, \code{inad_fit}, or common parent class \code{ad_fit}.
#' @param ... Unused; present for S3 consistency.
#'
#' @return A scalar of class \code{"logLik"} with attributes \code{df}
#'   and \code{nobs}.
#'
#' @name logLik.antedep
NULL

#' @rdname logLik.antedep
#' @importFrom stats logLik
#' @export
logLik.ad_fit <- function(object, ...) .ad_logLik(object)

#' @rdname logLik.antedep
#' @export
logLik.gau_fit <- function(object, ...) .ad_logLik(object)

#' @rdname logLik.antedep
#' @export
logLik.cat_fit <- function(object, ...) .ad_logLik(object)

#' @rdname logLik.antedep
#' @export
logLik.inad_fit <- function(object, ...) .ad_logLik(object)

#' Number of Subjects for Antedependence Model Fits
#'
#' @param object A fitted model object.
#' @param ... Unused.
#'
#' @return Integer scalar: number of subjects.
#'
#' @name nobs.antedep
NULL

#' @rdname nobs.antedep
#' @importFrom stats nobs
#' @export
nobs.ad_fit <- function(object, ...) .ad_nobs(object)

#' @rdname nobs.antedep
#' @export
nobs.gau_fit <- function(object, ...) .ad_nobs(object)

#' @rdname nobs.antedep
#' @export
nobs.cat_fit <- function(object, ...) .ad_nobs(object)

#' @rdname nobs.antedep
#' @export
nobs.inad_fit <- function(object, ...) .ad_nobs(object)

#' Extract Model Coefficients from Antedependence Fits
#'
#' By default, returns a named numeric vector so standard tools can consume
#' fitted objects. Use \code{type = "list"} to inspect the family-specific
#' structured representation.
#'
#' @param object A fitted model object.
#' @param type Return type: \code{"vector"} (default) or \code{"list"}.
#' @param ... Unused.
#'
#' @return A named numeric vector or a family-specific named list.
#'
#' @name coef.antedep
NULL

#' @rdname coef.antedep
#' @importFrom stats coef
#' @export
coef.gau_fit <- function(object, type = c("vector", "list"), ...) {
    type <- match.arg(type)
    out <- list()
    if (!is.null(object$mu)) out$mu <- object$mu
    if (!is.null(object$sigma)) out$sigma <- object$sigma
    if (!is.null(object$phi)) out$phi <- object$phi
    if (!is.null(object$tau) && any(as.numeric(object$tau) != 0)) {
        tau <- as.numeric(object$tau)
        if (length(tau) > 1L) out$tau <- tau[-1L]
    }
    if (type == "list") return(out)

    ord <- object$settings$order
    parts <- list(.ad_named_vector(object$mu, "mu"),
                  .ad_named_vector(object$sigma, "sigma"))
    if (ord == 1L && !is.null(object$phi)) {
        parts[[length(parts) + 1L]] <- .ad_named_vector(object$phi, "phi1", drop_first = TRUE)
    } else if (ord == 2L && !is.null(object$phi)) {
        phi <- as.matrix(object$phi)
        if (nrow(phi) >= 1L) parts[[length(parts) + 1L]] <- .ad_named_vector(phi[1L, ], "phi1", drop_first = TRUE)
        if (nrow(phi) >= 2L) parts[[length(parts) + 1L]] <- .ad_named_vector(phi[2L, ], "phi2", drop_first = TRUE)
    }
    if (!is.null(out$tau)) {
        parts[[length(parts) + 1L]] <- stats::setNames(out$tau, paste0("tau[", seq_along(out$tau) + 1L, "]"))
    }
    unlist(parts, use.names = TRUE)
}

#' @rdname coef.antedep
#' @export
coef.cat_fit <- function(object, type = c("vector", "list"), ...) {
    type <- match.arg(type)
    out <- list(marginal = object$marginal, transition = object$transition)
    if (type == "list") return(out)

    ci_vec <- tryCatch({
        tab <- .ad_ci_summary(ci_cat(fit = object))
        if (all(c("param", "est") %in% names(tab))) stats::setNames(tab$est, tab$param) else NULL
    }, error = function(e) NULL)
    if (!is.null(ci_vec)) return(ci_vec)

    .ad_flatten_list(out)
}

#' @rdname coef.antedep
#' @export
coef.inad_fit <- function(object, type = c("vector", "list"), ...) {
    type <- match.arg(type)
    out <- list()
    if (!is.null(object$alpha)) out$alpha <- object$alpha
    if (!is.null(object$theta)) out$theta <- object$theta
    if (!is.null(object$tau) && any(as.numeric(object$tau) != 0)) {
        tau <- as.numeric(object$tau)
        if (length(tau) > 1L) out$tau <- tau[-1L]
    }
    if (!is.null(object$nb_inno_size)) out$nb_inno_size <- object$nb_inno_size
    if (type == "list") return(out)

    ord <- object$settings$order
    parts <- list()
    if (ord == 1L && !is.null(object$alpha)) {
        parts[[length(parts) + 1L]] <- .ad_named_vector(object$alpha, "alpha", drop_first = TRUE)
    } else if (ord == 2L && !is.null(object$alpha)) {
        alpha <- as.matrix(object$alpha)
        if (ncol(alpha) >= 1L) parts[[length(parts) + 1L]] <- .ad_named_vector(alpha[, 1L], "alpha1", drop_first = TRUE)
        if (ncol(alpha) >= 2L) {
            a2 <- .ad_named_vector(alpha[, 2L], "alpha2", drop_first = TRUE)
            a2 <- a2[!grepl("^alpha2\\[2\\]$", names(a2))]
            parts[[length(parts) + 1L]] <- a2
        }
    }
    parts[[length(parts) + 1L]] <- .ad_named_vector(object$theta, "theta")
    if (!is.null(out$tau)) {
        parts[[length(parts) + 1L]] <- stats::setNames(out$tau, paste0("tau[", seq_along(out$tau) + 1L, "]"))
    }
    if (!is.null(object$nb_inno_size)) {
        parts[[length(parts) + 1L]] <- .ad_named_vector(object$nb_inno_size, "nb_inno_size")
    }
    unlist(parts, use.names = TRUE)
}

#' Fitted Values and Residuals for Gaussian AD Fits
#'
#' @param object A \code{gau_fit} object.
#' @param y Original data matrix.
#' @param ... Unused.
#'
#' @return A numeric matrix of fitted values or residuals.
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
    mu <- as.numeric(object$mu)
    yhat <- matrix(NA_real_, nrow = n, ncol = T)

    tau_vec <- rep(0, n)
    if (!is.null(object$settings$blocks) && !is.null(object$tau)) {
        tau_vec <- as.numeric(object$tau)[object$settings$blocks]
    }
    m <- function(t) mu[t] + tau_vec

    yhat[, 1] <- m(1)
    if (ord == 0) {
        for (t in seq_len(T)) yhat[, t] <- m(t)
        return(yhat)
    }
    if (ord == 1) {
        phi <- as.numeric(object$phi)
        for (t in 2:T) yhat[, t] <- m(t) + phi[t] * (y[, t - 1] - m(t - 1))
        return(yhat)
    }
    if (ord == 2) {
        phi <- as.matrix(object$phi)
        yhat[, 2] <- m(2) + phi[1, 2] * (y[, 1] - m(1))
        for (t in 3:T) {
            yhat[, t] <- m(t) +
                phi[1, t] * (y[, t - 1] - m(t - 1)) +
                phi[2, t] * (y[, t - 2] - m(t - 2))
        }
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

#' Convert Antedependence Simulation Output to a Time-Series Object
#'
#' @param x A matrix returned by \code{\link{simulate_gau}},
#'   \code{\link{simulate_cat}}, or \code{\link{simulate_inad}}.
#' @param ... Additional arguments passed to \code{\link[stats]{ts}}.
#'
#' @return A \code{ts} or \code{mts} object.
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

#' Deviance for Antedependence Model Fits
#'
#' @param object A fitted model object.
#' @param ... Unused.
#'
#' @return A scalar deviance value.
#'
#' @name deviance.antedep
NULL

#' @rdname deviance.antedep
#' @importFrom stats deviance
#' @export
deviance.ad_fit <- function(object, ...) -2 * as.numeric(logLik(object))

#' @rdname deviance.antedep
#' @export
deviance.gau_fit <- function(object, ...) deviance.ad_fit(object, ...)

#' @rdname deviance.antedep
#' @export
deviance.cat_fit <- function(object, ...) deviance.ad_fit(object, ...)

#' @rdname deviance.antedep
#' @export
deviance.inad_fit <- function(object, ...) deviance.ad_fit(object, ...)

#' Confidence Intervals for Antedependence Model Fits
#'
#' Returns a standard two-column matrix with \code{lower} and \code{upper}
#' columns. The structured helpers \code{ci_gau()}, \code{ci_cat()}, and
#' \code{ci_inad()} remain available for detailed inspection.
#'
#' @param object A fitted model object.
#' @param parm Optional character vector of parameter names.
#' @param level Confidence level.
#' @param y Original data matrix, required for \code{inad_fit}.
#' @param ... Passed to the underlying \code{ci_*} helper.
#'
#' @return A matrix with columns \code{lower} and \code{upper}.
#'
#' @name confint.antedep
NULL

#' @rdname confint.antedep
#' @importFrom stats confint
#' @export
confint.gau_fit <- function(object, parm = NULL, level = 0.95, ...) {
    .ad_confint_matrix(ci_gau(fit = object, level = level, ...), parm = parm)
}

#' @rdname confint.antedep
#' @export
confint.cat_fit <- function(object, parm = NULL, level = 0.95, y = NULL, ...) {
    .ad_confint_matrix(ci_cat(fit = object, y = y, level = level, ...), parm = parm)
}

#' @rdname confint.antedep
#' @export
confint.inad_fit <- function(object, parm = NULL, level = 0.95, y, ...) {
    if (missing(y)) {
        stop("confint() for inad_fit requires the original data matrix 'y'.", call. = FALSE)
    }
    .ad_confint_matrix(ci_inad(y = y, fit = object, level = level, ...), parm = parm)
}

#' Variance-Covariance Matrices for Antedependence Model Fits
#'
#' Returns a named diagonal matrix using the standard errors from the
#' corresponding confidence-interval helper. Entries are \code{NA} when a
#' standard error is not available.
#'
#' @param object A fitted model object.
#' @param y Original data matrix, required for \code{inad_fit}.
#' @param ... Passed to the underlying \code{ci_*} helper.
#'
#' @return A named square matrix.
#'
#' @name vcov.antedep
NULL

#' @rdname vcov.antedep
#' @importFrom stats vcov
#' @export
vcov.gau_fit <- function(object, ...) {
    .ad_vcov_from_ci(coef(object), ci_gau(fit = object, ...))
}

#' @rdname vcov.antedep
#' @export
vcov.cat_fit <- function(object, y = NULL, ...) {
    .ad_vcov_from_ci(coef(object), ci_cat(fit = object, y = y, ...))
}

#' @rdname vcov.antedep
#' @export
vcov.inad_fit <- function(object, y, ...) {
    if (missing(y)) {
        stop("vcov() for inad_fit requires the original data matrix 'y'.", call. = FALSE)
    }
    .ad_vcov_from_ci(coef(object), ci_inad(y = y, fit = object, ...))
}

#' Likelihood-Ratio Comparisons for Nested AD Fits
#'
#' Compares one or more nested antedependence model fits in the order supplied.
#'
#' @param object First fitted model object.
#' @param ... Additional fitted model objects.
#' @param test Test label; only \code{"Chisq"} is currently supported.
#'
#' @details For INAD fits, a comparison between models of different orders is
#' a boundary problem (the added thinning coefficients are zero under the
#' smaller model), so no chi-square p-value is reported for it; see
#' \code{\link{test_order_inad}}.
#'
#' @return An \code{anova} data frame.
#'
#' @name anova.antedep
NULL

#' @rdname anova.antedep
#' @importFrom stats anova
#' @export
anova.gau_fit <- function(object, ..., test = "Chisq") .ad_anova(object, ..., test = test)

#' @rdname anova.antedep
#' @export
anova.cat_fit <- function(object, ..., test = "Chisq") .ad_anova(object, ..., test = test)

#' @rdname anova.antedep
#' @export
anova.inad_fit <- function(object, ..., test = "Chisq") {
    out <- .ad_anova(object, ..., test = test)
    ord <- vapply(c(list(object), list(...)), function(f) as.numeric(f$settings$order), numeric(1))
    boundary <- c(FALSE, diff(ord) != 0)
    if (any(boundary)) {
        out[["Pr(>Chisq)"]][boundary] <- NA_real_
        attr(out, "heading") <- c(
            attr(out, "heading"),
            paste("Order comparisons are boundary problems; no chi-square p-value is given",
                  "for them (see test_order_inad).")
        )
    }
    out
}
