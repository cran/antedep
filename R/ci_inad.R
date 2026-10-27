#' Confidence Intervals for Fitted INAD Models
#'
#' Computes Wald confidence intervals for the free parameters of a fitted INAD
#' model from the observed information.
#'
#' With \code{method = "louis"} (the default), the observed information is
#' computed exactly by Louis' identity: for every time point the latent thinned
#' counts have a finite conditional distribution given the data, so the
#' conditional moments of the complete-data score and Hessian are finite sums.
#' The time-point contributions are combined over all free parameters jointly,
#' including block effects and parameters shared across time by
#' \code{fit_inad(constrain = ...)}. With \code{method = "hessian"}, the
#' observed information is instead the negative Hessian of the observed-data
#' log-likelihood computed numerically (Richardson extrapolation with fixed
#' absolute steps; positive parameters on the log scale); the two methods agree
#' up to numerical error, except that numerical differentiation is unreliable
#' for estimates very close to zero, for which the Louis calculation should be
#' used. Estimates on the boundary of the parameter space (a thinning
#' coefficient at zero or, for binomial thinning, at one, or a negative binomial
#' innovation size at its cap \code{nb_inno_size_ub}) are held fixed and
#' receive no Wald interval. The information is computed from exact finite
#' sums (up to floating-point arithmetic); the resulting Wald intervals are
#' large-sample approximations.
#'
#' By default, intervals for the positive parameters (\code{alpha},
#' \code{theta}, \code{nb_inno_size}) are formed on the log scale and
#' transformed back; intervals for block effects use their natural scale. For
#' unconstrained fits with block effects, profile-likelihood intervals for
#' \code{tau} are also computed (element \code{tau_profile}).
#'
#' @param y Integer matrix with \code{n_subjects} rows and \code{n_time} columns.
#' @param fit A fitted model object returned by \code{\link{fit_inad}} or
#'   \code{\link{em_inad}}.
#' @param blocks Optional integer vector of length \code{n_subjects}. Defaults to
#'   the block assignment stored in \code{fit}.
#' @param level Confidence level between 0 and 1.
#' @param idx_time Optional integer vector of time indices whose time-specific
#'   parameters are reported. Default is all time points. Time-invariant
#'   parameters and block effects are always reported.
#' @param ridge Nonnegative ridge value added to the observed information
#'   matrix.
#' @param profile_maxeval Maximum number of function evaluations used in the
#'   profile likelihood refits for \code{tau}.
#' @param profile_xtol_rel Relative tolerance used in the profile likelihood refits.
#' @param method \code{"louis"} (default) or \code{"hessian"}.
#' @param log_scale Logical; if \code{TRUE} (default), intervals for positive
#'   parameters are computed on the log scale.
#' @param profile_tau Logical; compute profile-likelihood intervals for block
#'   effects of unconstrained fits. Default \code{TRUE}.
#' @param hessian_step Absolute step for \code{method = "hessian"}.
#'
#' @return An object of class \code{inad_ci}, a list with elements
#'   \code{settings}, \code{level}, \code{method}, \code{alpha}, \code{theta},
#'   \code{nb_inno_size}, \code{tau} (Wald), \code{tau_profile} (when computed),
#'   and \code{vcov}, the estimated covariance matrix of the free parameters.
#'   Each interval element is a data frame with columns \code{param},
#'   \code{est}, \code{se}, \code{lower}, \code{upper}, \code{width},
#'   \code{level}, and \code{scale}.
#'
#' @references
#' Louis, T.A. (1982). Finding the observed information matrix when using the
#' EM algorithm. \emph{Journal of the Royal Statistical Society, Series B}, 44,
#' 226-233.
#'
#' @examples
#' \donttest{
#' data("bolus_inad", package = "antedep")
#' y <- bolus_inad$y[, 1:4]  # first four occasions, to keep the example short
#' blocks <- bolus_inad$blocks
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", blocks = blocks)
#' ci <- ci_inad(y, fit, blocks = blocks, profile_tau = FALSE)
#' ci$alpha
#' ci$theta
#' ci$tau
#' }
#'
#' @export
ci_inad <- function(
        y,
        fit,
        blocks = NULL,
        level = 0.95,
        idx_time = NULL,
        ridge = 0,
        profile_maxeval = 2500,
        profile_xtol_rel = 1e-6,
        method = c("louis", "hessian"),
        log_scale = TRUE,
        profile_tau = TRUE,
        hessian_step = 1e-3
) {
    method <- match.arg(method)
    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "ci_inad currently supports complete data only. Missing-data INAD confidence intervals are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) stop("y must be a nonnegative integer matrix")

    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")
    if (is.null(fit$settings$thinning)) stop("fit$settings$thinning is missing")
    if (is.null(fit$settings$innovation)) stop("fit$settings$innovation is missing")
    if (is.null(fit$theta)) stop("fit$theta is missing")
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (!is.null(fit$settings$na_action) && identical(fit$settings$na_action, "marginalize")) {
        stop(
            "ci_inad currently supports complete-data fits only. Missing-data INAD confidence intervals are not implemented yet.",
            call. = FALSE
        )
    }

    n <- nrow(y)
    N <- ncol(y)
    if (is.null(idx_time)) idx_time <- seq_len(N)
    idx_time <- as.integer(idx_time)
    if (any(idx_time < 1) || any(idx_time > N)) stop("idx_time out of range")

    if (is.null(blocks) && !is.null(fit$settings$blocks)) blocks <- fit$settings$blocks
    has_blocks <- !is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L
    if (has_blocks) {
        blocks <- .normalize_blocks(blocks, n)$blocks_id
    } else {
        blocks <- NULL
    }
    is_constrained <- !is.null(fit$settings$constraint_df)

    if (method == "louis") {
        li <- .louis_full_info(y, fit, blocks)
        info <- li$info; tab <- li$fp$tab
    } else {
        ni <- .inad_numeric_info(y, fit, blocks, step = hessian_step)
        info <- ni$info; tab <- ni$fp$tab
    }
    res <- .inad_ci_from_info(tab, info, level = level, log_scale = log_scale, ridge = ridge)

    keep_time <- function(df) {
        if (is.null(df)) return(NULL)
        k <- tab$index[match(df$param, tab$param)]
        df <- df[is.na(k) | k %in% idx_time, , drop = FALSE]
        rownames(df) <- NULL
        df
    }
    out <- list(settings = fit$settings, level = level, method = method,
                alpha = keep_time(res$alpha), theta = keep_time(res$theta),
                nb_inno_size = keep_time(res$nb_inno_size), tau = res$tau,
                vcov = res$vcov, info_eigen_range = res$info_eigen_range)

    if (has_blocks && profile_tau && !is_constrained && fit$settings$order %in% c(1, 2)) {
        out$tau_profile <- .ci_tau_profile_inad(
            y = y,
            fit = fit,
            blocks = blocks,
            level = level,
            maxeval = profile_maxeval,
            xtol_rel = profile_xtol_rel
        )
    }

    class(out) <- "inad_ci"
    out
}

#' Print method for INAD confidence intervals
#'
#' @param x An object of class \code{inad_ci}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.inad_ci <- function(x, ...) {
    cat("CI level:", x$level, "\n")
    if (!is.null(x$alpha)) cat("alpha rows:", nrow(x$alpha), "\n")
    if (!is.null(x$theta)) cat("theta rows:", nrow(x$theta), "\n")
    if (!is.null(x$nb_inno_size)) cat("nb_inno_size rows:", nrow(x$nb_inno_size), "\n")
    if (!is.null(x$tau)) cat("tau rows:", nrow(x$tau), "\n")
    invisible(x)
}

#' Summary method for inad_ci objects
#'
#' @param object An \code{inad_ci} object.
#' @param ... Unused.
#'
#' @return A data frame stacking available confidence intervals with columns
#'   \code{param}, \code{component}, \code{est}, \code{se}, \code{lower},
#'   \code{upper}, and \code{level}. Returns \code{NULL} if no intervals are available.
#' @export
summary.inad_ci <- function(object, ...) {
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

    parts[[length(parts) + 1L]] <- add_part(object$alpha, "alpha")
    parts[[length(parts) + 1L]] <- add_part(object$theta, "theta")
    parts[[length(parts) + 1L]] <- add_part(object$nb_inno_size, "nb_inno_size")
    parts[[length(parts) + 1L]] <- add_part(object$tau, "tau")
    parts <- Filter(Negate(is.null), parts)

    if (length(parts) == 0L) return(NULL)
    do.call(rbind, parts)
}

#' @keywords internal
.ci_tau_profile_inad <- function(
        y,
        fit,
        blocks,
        level,
        maxeval,
        xtol_rel,
        expand = 2,
        max_bracket_iter = 50
) {
    if (!requireNamespace("nloptr", quietly = TRUE)) stop("nloptr is required for tau profiling")

    if (is.null(fit$tau) || length(fit$tau) < 2) stop("fit$tau must exist with at least two blocks")
    if (is.null(fit$alpha)) stop("fit$alpha is missing for profiling refits")
    if (is.null(fit$theta)) stop("fit$theta is missing for profiling refits")

    n <- nrow(y)
    N <- ncol(y)
    B <- max(blocks)

    if (length(fit$tau) != B) stop("length(fit$tau) must equal max(blocks)")
    if (length(fit$theta) != N) stop("length(fit$theta) must equal ncol(y)")

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    eps_pos <- 1e-10
    eps_alpha <- 1e-10
    alpha_ub <- 1 - eps_alpha
    alpha_share_tol <- 1e-6

    # If order-1 alpha is supplied as time-constant, profile tau under the same
    # shared-alpha nuisance structure (rather than a fully time-varying alpha).
    alpha_profile_shared <- FALSE
    if (ord == 1 && N >= 2 && length(fit$alpha) == N) {
        a <- as.numeric(fit$alpha)
        if (all(is.finite(a[2:N])) && max(abs(a[2:N] - a[2])) <= alpha_share_tol) {
            alpha_profile_shared <- TRUE
        }
    }

    lambda_ok <- function(theta_vec, tau_vec) {
        mu <- if (innovation == "bell") {
            outer(rep(1, n), .bell_mean_from_theta(theta_vec)) + matrix(tau_vec[blocks], n, N)
        } else {
            outer(rep(1, n), theta_vec) + matrix(tau_vec[blocks], n, N)
        }
        all(is.finite(mu)) && all(mu > eps_pos)
    }

    pack_par <- function(alpha, theta, tau, nb_inno_size, ord, B, N, innovation,
                         tau_fix_idx, alpha_profile_shared = FALSE) {
        out <- numeric(0)
        tau_free_idx <- 2:B
        if (length(tau_fix_idx) > 0) tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
        if (length(tau_free_idx) > 0) out <- c(out, tau[tau_free_idx])

        if (ord == 1) {
            if (N >= 2) {
                if (alpha_profile_shared) {
                    out <- c(out, alpha[2])
                } else {
                    out <- c(out, alpha[2:N])
                }
            }
        }
        if (ord == 2) {
            if (N >= 2) out <- c(out, alpha[2:N, 1])
            if (N >= 3) out <- c(out, alpha[3:N, 2])
        }

        out <- c(out, theta)

        # nb_inno_size is NOT included: held fixed at its full-model MLE during
        # profile refits (Variant 1 / constrained-fit paradigm).
        out
    }

    unpack_par <- function(par, alpha0, theta0, tau0, nb0, ord, B, N, innovation,
                           tau_fix_idx, tau_fix_val, alpha_profile_shared = FALSE) {
        k <- 0

        tau <- tau0
        tau[1] <- 0
        tau_free_idx <- 2:B
        if (length(tau_fix_idx) > 0) tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
        if (length(tau_free_idx) > 0) {
            tau[tau_free_idx] <- par[(k + 1):(k + length(tau_free_idx))]
            k <- k + length(tau_free_idx)
        }
        if (length(tau_fix_idx) > 0) {
            tau[tau_fix_idx] <- tau_fix_val
            tau[1] <- 0
        }

        alpha <- alpha0
        if (ord == 1) {
            if (N >= 2) {
                if (alpha_profile_shared) {
                    alpha[2:N] <- par[k + 1]
                    k <- k + 1
                } else {
                    alpha[2:N] <- par[(k + 1):(k + (N - 1))]
                    k <- k + (N - 1)
                }
            }
            alpha[1] <- 0
            alpha <- pmin(pmax(alpha, 0), alpha_ub)
            alpha[1] <- 0
        }
        if (ord == 2) {
            if (N >= 2) {
                alpha[2:N, 1] <- par[(k + 1):(k + (N - 1))]
                k <- k + (N - 1)
            }
            if (N >= 3) {
                alpha[3:N, 2] <- par[(k + 1):(k + (N - 2))]
                k <- k + (N - 2)
            }
            alpha[1, 1] <- 0
            alpha[1, 2] <- 0
            if (N >= 2) alpha[2, 2] <- 0
            alpha <- pmin(pmax(alpha, 0), alpha_ub)
            alpha[1, 1] <- 0
            alpha[1, 2] <- 0
            if (N >= 2) alpha[2, 2] <- 0
        }

        theta <- par[(k + 1):(k + N)]
        k <- k + N
        theta <- pmax(theta, eps_pos)

        # nb_inno_size is held fixed at its full-model MLE (nb0) throughout the
        # profile refit.  Re-optimising it as a nuisance parameter (Variant 2)
        # would widen the CI inconsistently with the constrained-fit paradigm
        # used everywhere else in the package and can yield intervals that
        # contradict the LRT p-value.
        nb_inno_size <- nb0

        list(alpha = alpha, theta = theta, tau = tau, nb_inno_size = nb_inno_size)
    }

    refit_tau_fixed <- function(tau_fix_idx, tau_fix_val) {
        alpha0 <- fit$alpha
        theta0 <- pmax(as.numeric(fit$theta), eps_pos)

        tau0 <- as.numeric(fit$tau)
        tau0[1] <- 0

        nb0 <- fit$nb_inno_size
        if (innovation == "nbinom") {
            if (is.null(nb0)) nb0 <- rep(1, N)
            nb0 <- as.numeric(nb0)
            if (length(nb0) == 1L) nb0 <- rep(nb0, N)
            if (length(nb0) != N) stop("fit$nb_inno_size must be length 1 or ncol(y)")
            nb0 <- pmax(nb0, eps_pos)
        } else {
            nb0 <- NULL
        }

        par0 <- pack_par(alpha0, theta0, tau0, nb0, ord, B, N, innovation,
                         tau_fix_idx, alpha_profile_shared = alpha_profile_shared)

        obj <- function(par) {
            cur <- unpack_par(
                par = par,
                alpha0 = alpha0,
                theta0 = theta0,
                tau0 = tau0,
                nb0 = nb0,
                ord = ord,
                B = B,
                N = N,
                innovation = innovation,
                tau_fix_idx = tau_fix_idx,
                tau_fix_val = tau_fix_val,
                alpha_profile_shared = alpha_profile_shared
            )
            if (!lambda_ok(cur$theta, cur$tau)) return(1e8)

            ll <- tryCatch(
                logL_inad(
                    y = y,
                    order = ord,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = cur$alpha,
                    theta = cur$theta,
                    nb_inno_size = cur$nb_inno_size,
                    blocks = blocks,
                    tau = cur$tau
                ),
                error = function(e) NA_real_
            )
            if (!is.finite(ll)) return(1e8)
            -ll
        }

        opt <- nloptr::nloptr(
            x0 = par0,
            eval_f = obj,
            opts = list(
                algorithm = "NLOPT_LN_BOBYQA",
                xtol_rel = xtol_rel,
                maxeval = maxeval
            )
        )

        -opt$objective
    }

    q_chi_1 <- qchisq(level, df = 1)
    prof_target <- function(b, tval) {
        ll_prof <- refit_tau_fixed(tau_fix_idx = b, tau_fix_val = tval)
        2 * (fit$log_l - ll_prof) - q_chi_1
    }

    bracket_one_side <- function(b, direction, step0) {
        tau_mle <- fit$tau[b]
        lower_bound <- if (innovation == "bell") {
            -min(.bell_mean_from_theta(fit$theta))
        } else {
            -min(fit$theta)
        }

        x1 <- tau_mle
        step <- step0 * direction

        for (k in 1:max_bracket_iter) {
            x2 <- x1 + step
            # Enforce only the physical lower bound (innovation mean must be
            # positive); no artificial upper cap so the search is free to expand
            # as far as the profile likelihood requires.
            x2 <- max(lower_bound, x2)
            f2 <- prof_target(b, x2)
            if (is.finite(f2) && f2 >= 0) return(sort(c(x1, x2)))
            if (x2 == lower_bound) break   # hit hard physical constraint
            x1 <- x2
            step <- step * expand
        }

        NULL
    }

    tau_ci <- vector("list", length = B - 1)
    for (b in 2:B) {
        tau_mle <- fit$tau[b]
        lower_bound <- if (innovation == "bell") {
            -min(.bell_mean_from_theta(fit$theta))
        } else {
            -min(fit$theta)
        }
        # Initial step: 20 % of |tau_mle| (min 0.1).  With expand = 2 and
        # max_bracket_iter = 50 the bracket can reach tau_mle ± step0*(2^50-1),
        # far beyond any realistic CI bound.
        step0 <- max(0.1, abs(tau_mle) * 0.2)
        if (!is.finite(step0) || step0 <= 0) step0 <- 0.1

        f_mle <- prof_target(b, tau_mle)
        if (!is.finite(f_mle) || f_mle > 1e-6) stop("Profile target at MLE not near zero for tau[", b, "]")

        left_int <- bracket_one_side(b, direction = -1, step0 = step0)
        right_int <- bracket_one_side(b, direction = 1, step0 = step0)

        left <- NA_real_
        right <- NA_real_

        if (!is.null(left_int)) left <- uniroot(function(x) prof_target(b, x), interval = left_int)$root
        if (!is.null(right_int)) right <- uniroot(function(x) prof_target(b, x), interval = right_int)$root

        # Report tau SE on the same profile-likelihood scale as its CI. This
        # keeps tabular SE and CI summaries aligned for the block effect.
        se <- NA_real_
        if (is.finite(left) && is.finite(right)) {
            z <- qnorm(1 - (1 - level) / 2)
            if (is.finite(z) && z > 0) se <- (right - left) / (2 * z)
        }

        tau_ci[[b - 1]] <- data.frame(
            param = paste0("tau[", b, "]"),
            est = tau_mle,
            se = se,
            lower = left,
            upper = right,
            width = if (is.finite(left) && is.finite(right)) right - left else NA_real_,
            level = level,
            row.names = NULL
        )
    }

    do.call(rbind, tau_ci)
}

