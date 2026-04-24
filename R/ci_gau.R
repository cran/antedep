#' Confidence intervals for fitted Gaussian AD models
#'
#' Computes approximate Wald confidence intervals for selected parameters from a
#' fitted Gaussian AD model.
#'
#' @param fit A fitted model object returned by \code{\link{fit_gau}}.
#' @param level Confidence level between 0 and 1.
#' @param parameters Which parameters to include: \code{"all"} (default),
#'   \code{"mu"}, \code{"phi"}, or \code{"sigma"}.
#'
#' @return An object of class \code{gau_ci}, a list with elements
#'   \code{settings}, \code{level}, \code{mu}, \code{phi}, and \code{sigma}.
#'   Each non-NULL element is a data frame with columns \code{param},
#'   \code{est}, \code{se}, \code{lower}, \code{upper}, and \code{level}.
#'
#' @details
#' This helper currently supports complete-data Gaussian AD fits.
#'
#' Standard errors are based on large-sample approximations:
#' \itemize{
#'   \item \eqn{SE(\hat{\mu}_t) \approx \hat{\sigma}_t / \sqrt{n}}
#'   \item \eqn{SE(\hat{\sigma}_t) \approx \hat{\sigma}_t / \sqrt{2n}}
#'   \item \eqn{SE(\hat{\phi}) \approx \sqrt{(1-\hat{\phi}^2)/n}} for free \eqn{\phi} entries
#' }
#'
#' @examples
#' \donttest{
#' y <- simulate_gau(n_subjects = 80, n_time = 6, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' ci <- ci_gau(fit)
#' ci$mu
#' ci$phi
#' ci$sigma
#' }
#'
#' @seealso \code{\link{fit_gau}}, \code{\link{ci_cat}}, \code{\link{ci_inad}}
#' @export
ci_gau <- function(fit, level = 0.95, parameters = "all") {
    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")
    if (is.null(fit$settings$n_subjects)) stop("fit$settings$n_subjects is missing")
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (level <= 0 || level >= 1) stop("level must be between 0 and 1")

    if (!is.null(fit$n_missing) && fit$n_missing > 0) {
        stop(
            "ci_gau currently supports complete-data fits only. Missing-data Gaussian AD confidence intervals are not implemented yet.",
            call. = FALSE
        )
    }
    if (!is.null(fit$settings$na_action) && identical(fit$settings$na_action, "em")) {
        stop(
            "ci_gau currently supports complete-data fits only. Missing-data Gaussian AD confidence intervals are not implemented yet.",
            call. = FALSE
        )
    }

    parameters <- match.arg(parameters, c("all", "mu", "phi", "sigma"))
    z <- stats::qnorm(1 - (1 - level) / 2)

    n <- fit$settings$n_subjects
    ord <- fit$settings$order
    N <- if (!is.null(fit$sigma)) length(as.numeric(fit$sigma)) else length(as.numeric(fit$mu))

    out <- list(settings = fit$settings, level = level, mu = NULL, phi = NULL, sigma = NULL)

    if (parameters %in% c("all", "mu") && !is.null(fit$mu)) {
        mu_hat <- as.numeric(fit$mu)
        sigma_ref <- if (!is.null(fit$sigma)) as.numeric(fit$sigma) else rep(stats::sd(mu_hat), length(mu_hat))
        se_mu <- sigma_ref / sqrt(n)
        out$mu <- data.frame(
            param = paste0("mu[", seq_along(mu_hat), "]"),
            est = mu_hat,
            se = se_mu,
            lower = mu_hat - z * se_mu,
            upper = mu_hat + z * se_mu,
            level = level,
            row.names = NULL
        )
    }

    if (parameters %in% c("all", "sigma") && !is.null(fit$sigma)) {
        sigma_hat <- as.numeric(fit$sigma)
        se_sigma <- sigma_hat / sqrt(2 * n)
        out$sigma <- data.frame(
            param = paste0("sigma[", seq_along(sigma_hat), "]"),
            est = sigma_hat,
            se = se_sigma,
            lower = pmax(0, sigma_hat - z * se_sigma),
            upper = sigma_hat + z * se_sigma,
            level = level,
            row.names = NULL
        )
    }

    if (parameters %in% c("all", "phi")) {
        phi_rows <- list()

        if (ord == 1 && !is.null(fit$phi)) {
            phi_hat <- as.numeric(fit$phi)
            if (length(phi_hat) != N) stop("fit$phi has unexpected length for order 1")

            if (N >= 2) for (t in 2:N) {
                est <- phi_hat[t]
                se <- sqrt(max(1e-12, 1 - est^2) / n)
                phi_rows[[length(phi_rows) + 1L]] <- data.frame(
                    param = paste0("phi1[", t, "]"),
                    est = est,
                    se = se,
                    lower = est - z * se,
                    upper = est + z * se,
                    level = level,
                    row.names = NULL
                )
            }
        }

        if (ord == 2 && !is.null(fit$phi)) {
            phi_hat <- as.matrix(fit$phi)
            if (nrow(phi_hat) != 2 || ncol(phi_hat) != N) stop("fit$phi has unexpected dimensions for order 2")

            if (N >= 2) for (t in 2:N) {
                est <- phi_hat[1, t]
                se <- sqrt(max(1e-12, 1 - est^2) / n)
                phi_rows[[length(phi_rows) + 1L]] <- data.frame(
                    param = paste0("phi1[", t, "]"),
                    est = est,
                    se = se,
                    lower = est - z * se,
                    upper = est + z * se,
                    level = level,
                    row.names = NULL
                )
            }

            if (N >= 3) for (t in 3:N) {
                est <- phi_hat[2, t]
                se <- sqrt(max(1e-12, 1 - est^2) / n)
                phi_rows[[length(phi_rows) + 1L]] <- data.frame(
                    param = paste0("phi2[", t, "]"),
                    est = est,
                    se = se,
                    lower = est - z * se,
                    upper = est + z * se,
                    level = level,
                    row.names = NULL
                )
            }
        }

        if (length(phi_rows) > 0) out$phi <- do.call(rbind, phi_rows)
    }

    class(out) <- "gau_ci"
    out
}

#' Print method for AD confidence intervals
#'
#' @param x An object of class \code{gau_ci}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.gau_ci <- function(x, ...) {
    cat("CI level:", x$level, "\n")
    if (!is.null(x$mu)) cat("mu rows:", nrow(x$mu), "\n")
    if (!is.null(x$phi)) cat("phi rows:", nrow(x$phi), "\n")
    if (!is.null(x$sigma)) cat("sigma rows:", nrow(x$sigma), "\n")
    invisible(x)
}

#' Summary method for gau_ci objects
#'
#' @param object A \code{gau_ci} object.
#' @param ... Unused.
#'
#' @return A data frame stacking available confidence intervals with columns
#'   \code{param}, \code{component}, \code{est}, \code{se}, \code{lower},
#'   \code{upper}, and \code{level}. Returns \code{NULL} if no intervals are available.
#' @export
summary.gau_ci <- function(object, ...) {
    parts <- list()

    add_part <- function(df, component) {
        if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(NULL)
        df$component <- component
        keep <- intersect(
            c("param", "component", "est", "se", "lower", "upper", "level"),
            names(df)
        )
        df[, keep, drop = FALSE]
    }

    parts[[length(parts) + 1L]] <- add_part(object$mu, "mu")
    parts[[length(parts) + 1L]] <- add_part(object$phi, "phi")
    parts[[length(parts) + 1L]] <- add_part(object$sigma, "sigma")
    parts <- Filter(Negate(is.null), parts)

    if (length(parts) == 0L) return(NULL)
    do.call(rbind, parts)
}
