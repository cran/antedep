#' Confidence intervals for fitted INAD models
#'
#' Computes confidence intervals for selected parameters from a fitted INAD model.
#' For the fixed effect case, Wald intervals for time varying alpha and theta are
#' computed via Louis identity for supported thinning-innovation combinations.
#' For block effects tau, profile likelihood intervals are computed by fixing one
#' component of tau and re maximizing the log likelihood over nuisance parameters.
#' For negative binomial innovations, Wald intervals for the innovation size
#' parameter are computed using a one dimensional observed information approximation
#' per time point, holding other parameters fixed at their fitted values.
#'
#' @param y Integer matrix with \code{n_subjects} rows and \code{n_time} columns.
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#' @param blocks Optional integer vector of length \code{n_subjects}. Required for
#'   block effect intervals. If provided, should match \code{fit$settings$blocks}.
#' @param level Confidence level between 0 and 1.
#' @param idx_time Optional integer vector of time indices for which to compute
#'   intervals. Default is all time points.
#' @param ridge Nonnegative ridge value added to the observed information matrix
#'   used for Louis based Wald intervals.
#' @param profile_maxeval Maximum number of function evaluations used in the
#'   profile likelihood refits.
#' @param profile_xtol_rel Relative tolerance used in the profile likelihood refits.
#'
#' @return An object of class \code{inad_ci}, a list with elements
#'   \code{settings}, \code{level}, \code{alpha}, \code{theta}, \code{nb_inno_size},
#'   and \code{tau}. Each non NULL interval element is a data frame with columns
#'   \code{param}, \code{est}, \code{lower}, \code{upper}, and possibly \code{se}
#'   and \code{width}.
#'
#' @examples
#' \donttest{
#' data("bolus_inad", package = "antedep")
#' y <- bolus_inad$y
#' blocks <- bolus_inad$blocks
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "bell", blocks = blocks)
#' ci <- ci_inad(y, fit, blocks = blocks)
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
        profile_xtol_rel = 1e-6
) {
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

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    n <- nrow(y)
    N <- ncol(y)

    if (is.null(idx_time)) idx_time <- seq_len(N)
    idx_time <- as.integer(idx_time)
    if (any(idx_time < 1) || any(idx_time > N)) stop("idx_time out of range")

    has_blocks <- !is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 0L
    if (has_blocks) {
        blocks <- as.integer(blocks)
        if (length(blocks) != n) stop("blocks must have length nrow(y)")
        if (any(blocks < 1)) stop("blocks must be positive integers starting at 1")
    }

    out <- list(settings = fit$settings, level = level)

    # Compute alpha and theta CIs using Louis identity for supported combinations
    if (has_blocks && ord %in% c(1, 2)) {
        tmp <- .ci_alpha_theta_louis_fe(
            y = y,
            fit = fit,
            blocks = blocks,
            level = level,
            idx_time = idx_time,
            ridge = ridge,
            thinning = thinning,
            innovation = innovation
        )
        out$alpha <- tmp$alpha
        out$theta <- tmp$theta
    } else {
        out$alpha <- NULL
        out$theta <- NULL
    }

    if (innovation == "nbinom") {
        out$nb_inno_size <- .ci_nb_inno_size_obsinfo_1d(
            y = y,
            fit = fit,
            blocks = if (has_blocks) blocks else NULL,
            level = level,
            idx_time = idx_time
        )
    } else {
        out$nb_inno_size <- NULL
    }

    if (has_blocks) {
        out$tau <- .ci_tau_profile_inad(
            y = y,
            fit = fit,
            blocks = blocks,
            level = level,
            maxeval = profile_maxeval,
            xtol_rel = profile_xtol_rel
        )
    } else {
        out$tau <- NULL
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
.lse <- function(v) {
    m <- max(v)
    m + log(sum(exp(v - m)))
}

#' Posterior computations for all thinning-innovation combinations
#'
#' Computes conditional expectations and variances needed for Louis' identity.
#' See Appendix C of the paper for the derivation.
#'
#' @keywords internal
.poster_general <- function(m, n, alpha, lambda, thinning, innovation, nb_inno_size = NULL) {
    alpha <- max(alpha, 1e-8)
    lambda <- max(lambda, 1e-8)

    # Compute posterior weights
    k <- 0:n
    u <- n - k  # innovation values

    # Thinning log-probabilities
    if (thinning == "binom") {
        alpha <- min(alpha, 1 - 1e-8)
        lthin <- stats::dbinom(k, size = m, prob = alpha, log = TRUE)
    } else if (thinning == "pois") {
        lthin <- stats::dpois(k, lambda = alpha * m, log = TRUE)
    } else {  # nbinom
        p_nb <- 1 / (1 + alpha)
        lthin <- stats::dnbinom(k, size = m, prob = p_nb, log = TRUE)
    }

    # Innovation log-probabilities
    if (innovation == "pois") {
        linno <- stats::dpois(u, lambda = lambda, log = TRUE)
    } else if (innovation == "bell") {
        linno <- sapply(u, function(z) log(dbell(z, theta = lambda)))
    } else {  # nbinom
        if (is.null(nb_inno_size)) nb_inno_size <- 1
        linno <- stats::dnbinom(u, size = nb_inno_size, mu = lambda, log = TRUE)
    }

    lw <- lthin + linno
    lZ <- .lse(lw)
    w <- exp(lw - lZ)

    # Score functions for thinning parameter alpha
    # See Appendix A.1 for the conditional distributions
    if (thinning == "binom") {
        # Binomial thinning: k | Y ~ Bin(Y, alpha)
        # Score: k/alpha - (m-k)/(1-alpha)
        score_alpha <- k / alpha - (m - k) / (1 - alpha)
        # Second derivative (negative Hessian): k/alpha^2 + (m-k)/(1-alpha)^2
        h_alpha <- k / (alpha^2) + (m - k) / ((1 - alpha)^2)
    } else if (thinning == "pois") {
        # Poisson thinning: k | Y ~ Pois(alpha * Y)
        # Score: k/alpha - m
        score_alpha <- k / alpha - m
        # Second derivative: k/alpha^2
        h_alpha <- k / (alpha^2)
    } else {  # nbinom
        # NB thinning: k | Y ~ NB(Y, 1/(1+alpha))
        # p = 1/(1+alpha), so dp/dalpha = -1/(1+alpha)^2 = -p^2
        # Score and Hessian for NB(size=m, prob=p) w.r.t. alpha
        p <- 1 / (1 + alpha)
        q <- 1 - p
        score_alpha <- -m * p + k * (p^2 / q)
        h_alpha <- -m * p^2 + k * (p^3 * (2 - p) / (q^2))
    }

    # Score functions for innovation parameter lambda (theta)
    if (innovation == "pois") {
        # Poisson: score = u/lambda - 1
        score_lambda <- u / lambda - 1
        h_lambda <- u / (lambda^2)
    } else if (innovation == "bell") {
        # Bell: score = u/lambda - exp(lambda)
        score_lambda <- u / lambda - exp(lambda)
        h_lambda <- u / (lambda^2) + exp(lambda)
    } else {  # nbinom
        # NB(size=r, mu=lambda): score = u/lambda - (r+u)/(r+lambda)
        r <- if (is.null(nb_inno_size)) 1 else nb_inno_size
        score_lambda <- u / lambda - (r + u) / (r + lambda)
        h_lambda <- u / (lambda^2) - (r + u) / ((r + lambda)^2)
    }

    # Compute expectations
    E_sA <- sum(w * score_alpha)
    E_hA <- sum(w * h_alpha)
    E_sL <- sum(w * score_lambda)
    E_hL <- sum(w * h_lambda)

    E_sA2 <- sum(w * (score_alpha^2))
    E_sL2 <- sum(w * (score_lambda^2))
    E_sAsL <- sum(w * (score_alpha * score_lambda))

    Var_sA <- max(0, E_sA2 - E_sA^2)
    Var_sL <- max(0, E_sL2 - E_sL^2)
    Cov_sA_sL <- E_sAsL - E_sA * E_sL

    c(
        E_hA = E_hA,
        Var_sA = Var_sA,
        E_hL = E_hL,
        Var_sL = Var_sL,
        Cov_sA_sL = Cov_sA_sL
    )
}

#' Louis' identity observed information at time i
#'
#' @keywords internal
.louis_info_i_fe <- function(y, i, alpha_hat, theta_hat, tau_hat, blocks,
                              thinning, innovation, nb_inno_size = NULL) {
    n <- nrow(y)

    if (i == 1) {
        # At i=1, only innovation contributes
        Itt <- 0
        for (s in 1:n) {
            b <- blocks[s]
            z <- y[s, 1]
            if (innovation == "bell") {
                mu <- .bell_mean_from_theta(theta_hat) + tau_hat[b]
                la <- max(.bell_theta_from_mean(mu), 1e-8)
            } else {
                la <- max(theta_hat + tau_hat[b], 1e-8)
            }

            if (innovation == "pois") {
                Itt <- Itt + z / (la^2)
            } else if (innovation == "bell") {
                Itt <- Itt + z / (la^2) + exp(la)
            } else {  # nbinom
                r <- if (is.null(nb_inno_size)) 1 else nb_inno_size[1]
                Itt <- Itt + z / (la^2) - (r + z) / ((r + la)^2)
            }
        }
        mat <- matrix(Itt, 1, 1)
        dimnames(mat) <- list("theta[1]", "theta[1]")
        return(mat)
    }

    m <- y[, i - 1]
    nvec <- y[, i]

    a <- max(alpha_hat, 1e-8)
    th <- max(theta_hat, 1e-8)

    E_hA_sum <- 0
    Var_sA_sum <- 0
    E_hL_theta <- 0
    Var_sL_theta <- 0
    Cov_theta <- 0

    nb_sz <- if (innovation == "nbinom" && !is.null(nb_inno_size)) nb_inno_size[i] else NULL

    for (s in 1:n) {
        b <- blocks[s]
        if (innovation == "bell") {
            mu <- .bell_mean_from_theta(th) + tau_hat[b]
            la <- max(.bell_theta_from_mean(mu), 1e-8)
        } else {
            la <- max(th + tau_hat[b], 1e-8)
        }
        acc <- .poster_general(m = m[s], n = nvec[s], alpha = a, lambda = la,
                               thinning = thinning, innovation = innovation,
                               nb_inno_size = nb_sz)

        E_hA_sum <- E_hA_sum + acc["E_hA"]
        Var_sA_sum <- Var_sA_sum + acc["Var_sA"]
        E_hL_theta <- E_hL_theta + acc["E_hL"]
        Var_sL_theta <- Var_sL_theta + acc["Var_sL"]
        Cov_theta <- Cov_theta + acc["Cov_sA_sL"]
    }

    # Louis' identity: I_obs = E[B|Y] - Var[S|Y]
    Iaa <- E_hA_sum - Var_sA_sum
    Itt <- E_hL_theta - Var_sL_theta
    Iat <- -Cov_theta

    mat <- matrix(c(Iaa, Iat, Iat, Itt), 2, 2, byrow = TRUE)
    par_names <- c(paste0("alpha[", i, "]"), paste0("theta[", i, "]"))
    dimnames(mat) <- list(par_names, par_names)
    mat
}

#' Wald CI at time i using Louis' identity
#'
#' @keywords internal
.ci_wald_i_louis_fe <- function(y, i, alpha_hat, theta_hat, tau_hat, blocks,
                                 thinning, innovation, nb_inno_size = NULL,
                                 level = 0.95, ridge = 0) {
    Iobs <- .louis_info_i_fe(y, i, alpha_hat, theta_hat, tau_hat, blocks,
                              thinning, innovation, nb_inno_size)
    if (ridge > 0) Iobs <- Iobs + diag(ridge, nrow(Iobs))

    V <- tryCatch(solve(Iobs), error = function(e) NULL)
    if (is.null(V)) stop("Singular information at i = ", i, ".")

    z <- qnorm(1 - (1 - level) / 2)

    if (i == 1) {
        est <- c(theta_hat)
        nms <- c("theta[1]")
    } else {
        est <- c(alpha_hat, theta_hat)
        nms <- c(paste0("alpha[", i, "]"), paste0("theta[", i, "]"))
    }

    se <- sqrt(pmax(diag(V), 0))
    L <- est - z * se
    U <- est + z * se

    if (i != 1) {
        L[1] <- max(L[1], 0)
        U[1] <- max(U[1], 0)
    }

    data.frame(
        param = nms,
        est = est,
        se = se,
        lower = L,
        upper = U,
        width = U - L,
        level = level,
        row.names = NULL
    )
}

#' Compute alpha and theta CIs using Louis' identity (general version)
#'
#' @keywords internal
.ci_alpha_theta_louis_fe <- function(y, fit, blocks, level, idx_time, ridge,
                                      thinning, innovation) {
    alpha_rows <- list()
    theta_rows <- list()

    nb_inno_size <- fit$nb_inno_size

    for (i in idx_time) {
        if (i == 1) {
            tmp <- .ci_wald_i_louis_fe(
                y = y,
                i = 1,
                alpha_hat = NA,
                theta_hat = fit$theta[1],
                tau_hat = fit$tau,
                blocks = blocks,
                thinning = thinning,
                innovation = innovation,
                nb_inno_size = nb_inno_size,
                level = level,
                ridge = ridge
            )
            theta_rows[[length(theta_rows) + 1L]] <- tmp
        } else {
            alpha_i <- if (is.matrix(fit$alpha)) fit$alpha[i, 1] else fit$alpha[i]
            tmp <- .ci_wald_i_louis_fe(
                y = y,
                i = i,
                alpha_hat = alpha_i,
                theta_hat = fit$theta[i],
                tau_hat = fit$tau,
                blocks = blocks,
                thinning = thinning,
                innovation = innovation,
                nb_inno_size = nb_inno_size,
                level = level,
                ridge = ridge
            )
            alpha_rows[[length(alpha_rows) + 1L]] <- tmp[tmp$param == paste0("alpha[", i, "]"), , drop = FALSE]
            theta_rows[[length(theta_rows) + 1L]] <- tmp[tmp$param == paste0("theta[", i, "]"), , drop = FALSE]
        }
    }

    list(
        alpha = if (length(alpha_rows) > 0) do.call(rbind, alpha_rows) else NULL,
        theta = if (length(theta_rows) > 0) do.call(rbind, theta_rows) else NULL
    )
}

#' @keywords internal
.ci_nb_inno_size_obsinfo_1d <- function(y, fit, blocks, level, idx_time) {
    eps_pos <- 1e-10
    z <- qnorm(1 - (1 - level) / 2)

    if (is.null(fit$nb_inno_size)) return(NULL)

    nb <- as.numeric(fit$nb_inno_size)
    if (length(nb) == 1L) nb <- rep(nb, ncol(y))
    if (length(nb) != ncol(y)) stop("fit$nb_inno_size must be length 1 or ncol(y)")

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    if (innovation != "nbinom") return(NULL)

    res <- vector("list", length(idx_time))

    for (jj in seq_along(idx_time)) {
        i <- idx_time[jj]
        s0 <- max(nb[i], eps_pos)

        ll_at <- function(sz) {
            sz <- max(sz, eps_pos)
            nb_try <- nb
            nb_try[i] <- sz

            ll <- tryCatch(
                logL_inad(
                    y = y,
                    order = ord,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = fit$alpha,
                    theta = fit$theta,
                    nb_inno_size = nb_try,
                    blocks = blocks,
                    tau = if (!is.null(blocks)) fit$tau else NULL
                ),
                error = function(e) NA_real_
            )
            ll
        }

        h0 <- max(1e-4, 1e-3 * s0)
        f0 <- ll_at(s0)
        d2 <- NA_real_

        for (hfac in c(1, 0.1, 0.01)) {
            h <- h0 * hfac
            fp <- ll_at(s0 + h)
            fm <- ll_at(max(eps_pos, s0 - h))
            if (is.finite(f0) && is.finite(fp) && is.finite(fm)) {
                d2 <- (fp - 2 * f0 + fm) / (h^2)
                break
            }
        }

        if (is.na(d2)) {
            res[[jj]] <- data.frame(
                param = paste0("nb_inno_size[", i, "]"),
                est = s0,
                se = NA_real_,
                lower = NA_real_,
                upper = NA_real_,
                width = NA_real_,
                level = level,
                row.names = NULL
            )
            next
        }
        I <- max(0, -d2)

        if (I <= 0 || !is.finite(I)) {
            se <- NA_real_
            L <- NA_real_
            U <- NA_real_
        } else {
            se <- sqrt(1 / I)
            L <- max(0, s0 - z * se)
            U <- max(0, s0 + z * se)
        }

        res[[jj]] <- data.frame(
            param = paste0("nb_inno_size[", i, "]"),
            est = s0,
            se = se,
            lower = L,
            upper = U,
            width = if (is.finite(L) && is.finite(U)) (U - L) else NA_real_,
            level = level,
            row.names = NULL
        )
    }

    do.call(rbind, res)
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

        # Approximate SE from observed profile curvature at tau_mle.
        # Fall back to CI-width-based approximation if curvature is unavailable.
        se <- NA_real_
        h <- max(1e-4, min(0.1, abs(tau_mle) * 0.05))
        if (tau_mle - h > lower_bound + 1e-10) {
            ll_m <- refit_tau_fixed(tau_fix_idx = b, tau_fix_val = tau_mle - h)
            ll_0 <- fit$log_l - (f_mle + q_chi_1) / 2
            ll_p <- refit_tau_fixed(tau_fix_idx = b, tau_fix_val = tau_mle + h)
            if (is.finite(ll_m) && is.finite(ll_0) && is.finite(ll_p)) {
                d2 <- (ll_p - 2 * ll_0 + ll_m) / (h^2)
                info <- -d2
                if (is.finite(info) && info > 0) se <- sqrt(1 / info)
            }
        }
        if (!is.finite(se) && is.finite(left) && is.finite(right)) {
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
