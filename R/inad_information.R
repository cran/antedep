# Observed information for INAD models.
#
# Louis' identity gives the observed information exactly: for every time point
# the latent thinned counts have a finite conditional distribution, so the
# conditional moments of the complete-data score and Hessian are finite sums.
# The time-point contributions are mapped onto the free parameters of the fit,
# which may share parameters across time (time-invariance constraints) and
# include block effects. A numerical Hessian with fixed absolute steps is kept
# as an independent check.

#' Free parameters of a fitted INAD model
#'
#' @return A list with \code{tab}, a data frame with one row per free
#'   parameter (columns \code{component}, \code{group}, \code{param},
#'   \code{est}, \code{index}), and the unpacked parameter vectors.
#' @keywords internal
.inad_free_params <- function(fit, N, B) {
    ord <- fit$settings$order
    innovation <- fit$settings$innovation
    cinfo <- fit$settings$constraint_info
    if (is.null(cinfo)) {
        cinfo <- list(alpha1_const = FALSE, alpha2_const = FALSE,
                      theta_const = FALSE, kappa_const = FALSE)
    }
    au <- if (ord >= 1) .unpack_alpha(ord, fit$alpha, N) else list(a1 = NULL, a2 = NULL)
    theta <- rep_len(as.numeric(fit$theta), N)
    kappa <- if (innovation == "nbinom") rep_len(as.numeric(fit$nb_inno_size), N) else NULL
    tau <- if (B > 1L) as.numeric(fit$tau)[seq_len(B)] else 0

    rows <- list()
    add <- function(component, group, param, est, index) {
        rows[[length(rows) + 1L]] <<- data.frame(
            component = component, group = group, param = param, est = est,
            index = index, stringsAsFactors = FALSE
        )
    }
    lab1 <- if (ord == 2) "alpha1" else "alpha"
    if (ord >= 1) {
        if (cinfo$alpha1_const) add("alpha", "a1", lab1, au$a1[2], NA_integer_)
        else add("alpha", "a1", paste0(lab1, "[", 2:N, "]"), au$a1[2:N], 2:N)
    }
    if (ord == 2) {
        if (cinfo$alpha2_const) add("alpha", "a2", "alpha2", au$a2[3], NA_integer_)
        else add("alpha", "a2", paste0("alpha2[", 3:N, "]"), au$a2[3:N], 3:N)
    }
    if (cinfo$theta_const) add("theta", "theta", "theta", theta[1], NA_integer_)
    else add("theta", "theta", paste0("theta[", 1:N, "]"), theta, 1:N)
    if (!is.null(kappa)) {
        if (cinfo$kappa_const) add("nb_inno_size", "kappa", "kappa", kappa[1], NA_integer_)
        else add("nb_inno_size", "kappa", paste0("kappa[", 1:N, "]"), kappa, 1:N)
    }
    if (B > 1L) add("tau", "tau", paste0("tau[", 2:B, "]"), tau[2:B], 2:B)
    tab <- do.call(rbind, rows)
    # Estimates on the boundary of the parameter space (a thinning coefficient
    # at zero or, for binomial thinning, at one, an innovation parameter at
    # zero, or a negative binomial size at its upper bound) have no Wald
    # interval; they are held fixed in information calculations.
    kappa_ub <- if (!is.null(fit$settings$nb_inno_size_ub)) fit$settings$nb_inno_size_ub else 50
    is_alpha <- tab$group %in% c("a1", "a2")
    tab$boundary <- (is_alpha & tab$est <= 1e-8) |
        (is_alpha & identical(fit$settings$thinning, "binom") & tab$est >= 1 - 1e-8) |
        (tab$group == "theta" & tab$est <= 1e-8) |
        (tab$group == "kappa" & tab$est >= kappa_ub * (1 - 1e-6))
    list(tab = tab, a1 = au$a1, a2 = au$a2, theta = theta,
         kappa = kappa, tau = tau, cinfo = cinfo)
}

#' Row of the free-parameter table for a local (time or block) parameter
#' @keywords internal
.inad_free_index <- function(tab, group, index) {
    k <- which(tab$group == group & (is.na(tab$index) | tab$index == index))
    if (length(k) != 1L) stop("internal error: cannot map ", group, "[", index, "]")
    k
}

#' Write free-parameter values back into a fit-like list
#' @keywords internal
.inad_set_free <- function(fit, fp, values, N) {
    tab <- fp$tab
    a1 <- fp$a1; a2 <- fp$a2; theta <- fp$theta; kappa <- fp$kappa; tau <- fp$tau
    for (k in seq_len(nrow(tab))) {
        idx <- tab$index[k]
        v <- values[k]
        switch(tab$group[k],
               a1 = if (is.na(idx)) a1[2:N] <- v else a1[idx] <- v,
               a2 = if (is.na(idx)) a2[3:N] <- v else a2[idx] <- v,
               theta = if (is.na(idx)) theta[] <- v else theta[idx] <- v,
               kappa = if (is.na(idx)) kappa[] <- v else kappa[idx] <- v,
               tau = tau[idx] <- v)
    }
    ord <- fit$settings$order
    fit$alpha <- if (ord == 2) cbind(a1, a2) else a1
    fit$theta <- theta
    if (!is.null(kappa)) fit$nb_inno_size <- kappa
    if (length(tau) > 1L) fit$tau <- tau
    fit
}

#' Bell parameter for a given mean, to machine precision
#' @keywords internal
.bell_theta_newton <- function(mu, start = NULL) {
    th <- if (is.null(start) || !is.finite(start)) .bell_theta_from_mean(mu) else start
    for (k in seq_len(50)) {
        step <- (th * exp(th) - mu) / (exp(th) * (1 + th))
        th <- th - step
        if (abs(step) <= 1e-15 * max(1, abs(th))) break
    }
    th
}

#' Complete-data score and negative Hessian of one thinned count
#'
#' Valid for \code{0 < alpha} (and \code{alpha < 1} for binomial thinning). At
#' \code{alpha = 0}, or binomial \code{alpha = 1}, the thinned count is
#' degenerate and the conditional mean of this score is 0/0; see
#' \code{.louis_thin_boundary}.
#'
#' @return A list with the score \code{s} and the negative second derivative
#'   \code{b} of \code{log P(alpha odot yprev = z)} with respect to alpha.
#' @keywords internal
.louis_thin_terms <- function(z, yprev, alpha, thinning) {
    a <- alpha
    if (thinning == "binom") {
        list(s = z / a - (yprev - z) / (1 - a), b = z / a^2 + (yprev - z) / (1 - a)^2)
    } else if (thinning == "pois") {
        list(s = z / a - yprev, b = z / a^2)
    } else {
        list(s = z / a - (yprev + z) / (1 + a), b = z / a^2 - (yprev + z) / (1 + a)^2)
    }
}

#' Is a thinning coefficient on a bound where the thinned count is degenerate?
#' @keywords internal
.thin_at_bound <- function(alpha, thinning) {
    alpha == 0 || (thinning == "binom" && alpha == 1)
}

#' One-sided observed-data derivatives for a thinning coefficient on a bound
#'
#' At \code{alpha = 0} the thinned count \code{Z} is identically 0, and at
#' binomial \code{alpha = 1} it equals \code{yprev}; the support of the
#' complete-data distribution changes there, so the conditional mean of the
#' complete-data score is 0/0 and Fisher's identity cannot be evaluated at the
#' bound. The observed-data likelihood is \eqn{f(y) = \sum_z p_\alpha(z)
#' h(y - z)}, where \eqn{h} is the pmf of the other components of \eqn{y}.
#' With \eqn{n} = \code{yprev} and \eqn{H_k = h(y - k)}, the derivatives of
#' \eqn{p_\alpha} at \eqn{\alpha = 0} are \eqn{-n, n, 0, \ldots} for \eqn{z =
#' 0, 1, 2, \ldots} (all three thinning operators), so the right derivative is
#' \eqn{n (H_1 / H_0 - 1)}. The second derivatives are \eqn{c, -2c, c} for
#' \eqn{z = 0, 1, 2}, with \eqn{c = n(n - 1)}, \eqn{n^2}, and \eqn{n(n + 1)} for
#' binomial, Poisson, and negative binomial thinning. Binomial thinning at
#' \eqn{\alpha = 1} is the mirror image in \eqn{1 - \alpha}, with \eqn{z = n, n
#' - 1, n - 2}, giving the left derivative \eqn{-n (H_{n-1} / H_n - 1)}.
#'
#' @param z Values of the thinned count on the latent grid.
#' @param lp_other Log of the joint pmf of the other components at each grid
#'   point.
#' @return A named vector: the one-sided derivative \code{s} of
#'   \code{log f(y)} with respect to alpha, and \code{b}, minus the one-sided
#'   second derivative.
#' @keywords internal
.louis_thin_boundary <- function(z, lp_other, yprev, alpha, thinning) {
    n <- as.numeric(yprev)  # drop any name, which c(s = ...) would extend
    if (n == 0) return(c(s = 0, b = 0))
    lse_at <- function(k) {
        j <- z == k
        if (any(j)) .log_sum_exp(lp_other[j]) else -Inf
    }
    if (alpha == 0) {
        k0 <- 0; dir <- 1
    } else {
        k0 <- n; dir <- -1
    }
    l0 <- lse_at(k0)
    r1 <- exp(lse_at(k0 + dir) - l0)
    r2 <- exp(lse_at(k0 + 2 * dir) - l0)
    cc <- switch(thinning, binom = n * (n - 1), pois = n^2, n * (n + 1))
    d1 <- n * (r1 - 1)
    c(s = dir * d1, b = -(cc * (1 - 2 * r1 + r2) - d1^2))
}

#' Louis score and information for the parameters entering time i
#'
#' Local parameters are \code{a1} (i >= 2), \code{a2} (order 2, i >= 3),
#' \code{theta}, \code{kappa} (negative binomial innovations), and the block
#' effects \code{tau2}, ..., \code{tauB}. Returns the observed score (the
#' conditional mean of the complete-data score) and the observed information
#' \eqn{E(B \mid y) - \mathrm{Cov}(S \mid y)}, both summed over subjects.
#'
#' @keywords internal
.louis_time_terms <- function(y, i, ord, thinning, innovation, a1, a2, theta_i,
                              kappa_i, tau, blocks) {
    n <- nrow(y)
    B <- if (is.null(blocks)) 1L else length(tau)
    has_a1 <- ord >= 1 && i >= 2
    has_a2 <- ord == 2 && i >= 3
    is_nb <- innovation == "nbinom"
    tau_nm <- if (B > 1L) paste0("tau", 2:B) else character(0)
    nm <- c(if (has_a1) "a1", if (has_a2) "a2", "theta", if (is_nb) "kappa", tau_nm)
    P <- length(nm)
    info <- matrix(0, P, P, dimnames = list(nm, nm))
    score <- stats::setNames(numeric(P), nm)
    theta_bell <- if (innovation == "bell") theta_i else NA_real_

    for (s in seq_len(n)) {
        yi <- y[s, i]
        b <- if (is.null(blocks)) 1L else blocks[s]
        tb <- if (b > 1L) tau[b] else 0
        tnm <- if (b > 1L) paste0("tau", b) else NULL

        if (has_a2) {
            g <- expand.grid(z1 = 0:yi, z2 = 0:yi)
            g <- g[g$z1 + g$z2 <= yi, , drop = FALSE]
            z1 <- g$z1; z2 <- g$z2
        } else if (has_a1) {
            z1 <- 0:yi; z2 <- 0 * z1
        } else {
            z1 <- 0; z2 <- 0
        }
        w <- yi - z1 - z2

        if (innovation == "bell") {
            v <- if (tb == 0) theta_bell else
                .bell_theta_newton(.bell_mean_from_theta(theta_bell) + tb)
            par_inno <- v
        } else {
            m <- theta_i + tb
            par_inno <- m
        }

        lp_w <- log(.innov_vec(w, par_inno, innovation, if (is_nb) kappa_i else NA_real_))
        lp_1 <- if (has_a1) log(.thin_vec(z1, y[s, i - 1], a1[i], thinning)) else 0
        lp_2 <- if (has_a2) log(.thin_vec(z2, y[s, i - 2], a2[i], thinning)) else 0
        lp <- lp_w + lp_1 + lp_2
        lZ <- .log_sum_exp(lp)
        if (!is.finite(lZ)) stop("Zero likelihood contribution at time ", i, ", subject ", s)
        post <- exp(lp - lZ)

        S <- matrix(0, length(w), P, dimnames = list(NULL, nm))
        Bm <- matrix(0, P, P, dimnames = list(nm, nm))
        # A coefficient on a bound where its thinned count is degenerate gets
        # its exact one-sided score from .louis_thin_boundary and contributes
        # nothing to the complete-data score matrix S; its row and column of
        # the information hold only the one-sided curvature. Such coefficients
        # are held fixed or moved off the bound by a line search, never by a
        # Newton step, so their cross terms are not needed.
        bnd_score <- stats::setNames(numeric(P), nm)
        if (has_a1) {
            if (.thin_at_bound(a1[i], thinning)) {
                tb1 <- .louis_thin_boundary(z1, lp_w + lp_2, y[s, i - 1], a1[i], thinning)
                bnd_score["a1"] <- tb1[["s"]]; Bm["a1", "a1"] <- tb1[["b"]]
            } else {
                tt <- .louis_thin_terms(z1, y[s, i - 1], a1[i], thinning)
                S[, "a1"] <- tt$s; Bm["a1", "a1"] <- sum(post * tt$b)
            }
        }
        if (has_a2) {
            if (.thin_at_bound(a2[i], thinning)) {
                tb2 <- .louis_thin_boundary(z2, lp_w + lp_1, y[s, i - 2], a2[i], thinning)
                bnd_score["a2"] <- tb2[["s"]]; Bm["a2", "a2"] <- tb2[["b"]]
            } else {
                tt <- .louis_thin_terms(z2, y[s, i - 2], a2[i], thinning)
                S[, "a2"] <- tt$s; Bm["a2", "a2"] <- sum(post * tt$b)
            }
        }

        # Innovation terms. For Poisson and negative binomial innovations the
        # mean m = theta_i + tau_b, so d m / d theta_i = d m / d tau_b = 1.
        # The innovation pmf is zero for m <= 0, so m > 0 here, and the
        # conditional means below are exact for any m > 0.
        mean_nm <- c("theta", tnm)
        if (innovation == "pois") {
            mm <- m
            s_m <- w / mm - 1
            b_mm <- sum(post * w / mm^2)
            for (a in mean_nm) S[, a] <- s_m
            Bm[mean_nm, mean_nm] <- b_mm
        } else if (innovation == "nbinom") {
            mm <- m; k <- max(kappa_i, 1e-12)
            s_m <- w / mm - (w + k) / (k + mm)
            for (a in mean_nm) S[, a] <- s_m
            S[, "kappa"] <- digamma(w + k) - digamma(k) + log(k / (k + mm)) + (mm - w) / (k + mm)
            Bm[mean_nm, mean_nm] <- sum(post * (w / mm^2 - (w + k) / (k + mm)^2))
            Bm["kappa", "kappa"] <- -sum(post * (trigamma(w + k) - trigamma(k) + 1 / k -
                                                    1 / (k + mm) - (mm - w) / (k + mm)^2))
            b_mk <- sum(post * (mm - w) / (k + mm)^2)
            Bm["kappa", mean_nm] <- b_mk
            Bm[mean_nm, "kappa"] <- b_mk
        } else {
            # Bell: the Bell parameter v solves v exp(v) = theta exp(theta) + tau_b.
            th <- theta_bell
            q <- exp(v) * (1 + v); q1 <- exp(v) * (2 + v)
            c_th <- exp(th) * (1 + th) / q
            c_thth <- (exp(th) * (2 + th) - c_th^2 * q1) / q
            dl <- w / v - exp(v)
            d2l <- -w / v^2 - exp(v)
            S[, "theta"] <- dl * c_th
            Bm["theta", "theta"] <- -sum(post * (d2l * c_th^2 + dl * c_thth))
            if (!is.null(tnm)) {
                c_t <- 1 / q
                c_tt <- -q1 * c_t^3
                c_tht <- -q1 * c_th * c_t^2
                S[, tnm] <- dl * c_t
                Bm[tnm, tnm] <- -sum(post * (d2l * c_t^2 + dl * c_tt))
                Bm["theta", tnm] <- Bm[tnm, "theta"] <- -sum(post * (d2l * c_th * c_t + dl * c_tht))
            }
        }
        ES <- colSums(post * S)
        score <- score + ES + bnd_score
        info <- info + Bm - (crossprod(S, post * S) - tcrossprod(ES))
    }
    list(score = score, info = info)
}

#' Observed score and information for all free parameters, by Louis' identity
#'
#' @keywords internal
.louis_full_info <- function(y, fit, blocks = NULL) {
    N <- ncol(y)
    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    fp <- .inad_free_params(fit, N, B)
    tab <- fp$tab
    P <- nrow(tab)
    info <- matrix(0, P, P, dimnames = list(tab$param, tab$param))
    score <- stats::setNames(numeric(P), tab$param)
    for (i in seq_len(N)) {
        tt <- .louis_time_terms(y, i, ord, thinning, innovation, fp$a1, fp$a2, fp$theta[i],
                                if (!is.null(fp$kappa)) fp$kappa[i] else NULL, fp$tau, blocks)
        loc <- names(tt$score)
        map <- vapply(loc, function(p) {
            if (p %in% c("a1", "a2", "theta", "kappa")) .inad_free_index(tab, p, i)
            else .inad_free_index(tab, "tau", as.integer(sub("tau", "", p, fixed = TRUE)))
        }, integer(1))
        score[map] <- score[map] + tt$score
        info[map, map] <- info[map, map] + tt$info
    }
    list(score = score, info = info, fp = fp, B = B, blocks = blocks)
}

#' Wald intervals from an observed information matrix
#'
#' @keywords internal
.inad_ci_from_info <- function(tab, info, level = 0.95, log_scale = TRUE, ridge = 0) {
    inner <- if (is.null(tab$boundary)) rep(TRUE, nrow(tab)) else !tab$boundary
    if (any(!inner)) {
        warning("Estimates on the boundary of the parameter space (",
                paste(tab$param[!inner], collapse = ", "),
                ") are held fixed and have no Wald interval.", call. = FALSE)
    }
    I_in <- info[inner, inner, drop = FALSE]
    if (ridge > 0) I_in <- I_in + diag(ridge, nrow(I_in))
    I_in <- (I_in + t(I_in)) / 2
    V_in <- tryCatch(chol2inv(chol(I_in)), error = function(e) NULL)
    if (is.null(V_in)) stop("Observed information is not positive definite at the fitted values.")
    V <- matrix(NA_real_, nrow(tab), nrow(tab), dimnames = list(tab$param, tab$param))
    V[inner, inner] <- V_in
    info <- I_in
    z <- stats::qnorm(1 - (1 - level) / 2)
    se <- sqrt(pmax(diag(V), 0))
    positive <- tab$group != "tau"
    use_log <- log_scale & positive & tab$est > 0
    lower <- ifelse(use_log, tab$est * exp(-z * se / tab$est), tab$est - z * se)
    upper <- ifelse(use_log, tab$est * exp(z * se / tab$est), tab$est + z * se)
    lower[positive & !use_log] <- pmax(lower[positive & !use_log], 0)
    part <- function(component) {
        k <- tab$component == component
        if (!any(k)) return(NULL)
        data.frame(param = tab$param[k], est = tab$est[k], se = se[k], lower = lower[k],
                   upper = upper[k], width = upper[k] - lower[k], level = level,
                   scale = ifelse(use_log[k], "log", "natural"),
                   row.names = NULL, stringsAsFactors = FALSE)
    }
    list(alpha = part("alpha"), theta = part("theta"), nb_inno_size = part("nb_inno_size"),
         tau = part("tau"), vcov = V,
         info_eigen_range = range(eigen(info, symmetric = TRUE, only.values = TRUE)$values))
}

#' Numerical observed-data Hessian with fixed absolute steps
#'
#' Positive parameters are differentiated on the log scale and block effects on
#' their natural scale. Steps are absolute (\code{step}, refined by Richardson
#' extrapolation), because steps proportional to the parameter value become too
#' small for log-scale parameters near zero.
#'
#' @return The observed information for the free parameters on their natural
#'   scale.
#' @keywords internal
.inad_numeric_info <- function(y, fit, blocks = NULL, step = 1e-3) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("numDeriv is required for method = 'hessian'")
    N <- ncol(y)
    ord <- fit$settings$order
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    fp <- .inad_free_params(fit, N, B)
    tab <- fp$tab
    use_log <- tab$group != "tau" & tab$est > 0
    eta_hat <- ifelse(use_log, log(pmax(tab$est, 1e-300)), tab$est)
    inner <- which(!tab$boundary)  # boundary estimates are held fixed
    negll <- function(delta_inner) {
        delta <- numeric(nrow(tab)); delta[inner] <- delta_inner
        eta <- eta_hat + delta
        f <- .inad_set_free(fit, fp, ifelse(use_log, exp(eta), eta), N)
        ll <- tryCatch(
            logL_inad(y = y, order = ord, thinning = fit$settings$thinning,
                      innovation = fit$settings$innovation, alpha = f$alpha, theta = f$theta,
                      nb_inno_size = f$nb_inno_size, blocks = blocks,
                      tau = if (B > 1L) f$tau else 0),
            error = function(e) -Inf
        )
        if (!is.finite(ll)) return(1e20)
        -ll
    }
    H <- numDeriv::hessian(negll, rep(0, length(inner)), method = "Richardson",
                           method.args = list(eps = step, d = step, r = 4, v = 2))
    H <- (H + t(H)) / 2
    # Back to the natural scale (the score vanishes at the MLE, so only the
    # Jacobian enters). Rows and columns of boundary estimates are NA.
    jac <- ifelse(use_log, tab$est, 1)[inner]
    info <- matrix(NA_real_, nrow(tab), nrow(tab), dimnames = list(tab$param, tab$param))
    info[inner, inner] <- H / outer(jac, jac)
    list(info = info, fp = fp)
}

#' Observed-data log-likelihood of a fit-like list, as a function of the fit
#' @keywords internal
.inad_loglik_fun <- function(y, fit, blocks = NULL) {
    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    function(f) {
        tryCatch(
            logL_inad(y = y, order = ord, thinning = thinning, innovation = innovation,
                      alpha = f$alpha, theta = f$theta, nb_inno_size = f$nb_inno_size,
                      blocks = blocks, tau = if (B > 1L) f$tau else 0),
            error = function(e) -Inf
        )
    }
}

#' Line search along one free parameter, starting from a bound
#'
#' Maximizes the observed-data log-likelihood over the parameter in row
#' \code{k} of the free-parameter table, all other parameters fixed, between
#' its current value and \code{lim}. The search starts with a step of size
#' \code{t0} towards \code{lim}, shortens it until the log-likelihood
#' increases, lengthens it while the log-likelihood keeps increasing, and
#' finishes with \code{stats::optimize} on the bracket found.
#'
#' @return A list with the best value \code{par}, the log-likelihood increase
#'   \code{gain} (0 if no increase was found), and the corresponding \code{fit}.
#' @keywords internal
.inad_coord_search <- function(loglik, fit, fp, k, lim, t0, N, ll0 = loglik(fit)) {
    v0 <- fp$tab$est
    a <- v0[k]
    fit_at <- function(x) {
        v <- v0
        v[k] <- x
        .inad_set_free(fit, fp, v, N)
    }
    phi <- function(x) {
        ll <- loglik(fit_at(x))
        if (is.finite(ll)) ll - ll0 else -Inf
    }
    d <- lim - a
    none <- list(par = a, gain = 0, fit = fit)
    if (!is.finite(d) || d == 0) return(none)
    t <- sign(d) * min(abs(t0), abs(d))
    p1 <- phi(a + t)
    shrink <- 0
    while (!(p1 > 0) && shrink < 40) {
        t <- t / 4
        p1 <- phi(a + t)
        shrink <- shrink + 1
    }
    if (!(p1 > 0)) return(none)
    while (abs(2 * t) <= abs(d)) {
        p2 <- phi(a + 2 * t)
        if (!(p2 > p1)) break
        t <- 2 * t
        p1 <- p2
    }
    hi <- if (abs(2 * t) <= abs(d)) a + 2 * t else lim
    opt <- stats::optimize(phi, interval = sort(c(a, hi)), maximum = TRUE,
                           tol = 1e-10 * max(abs(t), 1e-12))
    if (is.finite(opt$objective) && opt$objective > p1) {
        best <- opt$maximum; gain <- opt$objective
    } else {
        best <- a + t; gain <- p1
    }
    list(par = best, gain = gain, fit = fit_at(best))
}

#' Newton polish of an INAD fit using the exact score and information
#'
#' Takes Newton-Raphson steps on the free parameters, using the Louis score and
#' observed information, with step halving so that the log-likelihood never
#' decreases and parameters stay in their parameter space. Used to finish
#' derivative-free fits to full precision.
#'
#' @return A list with the polished \code{fit} (parameters only) and its
#'   \code{log_l}.
#' @keywords internal
.inad_newton_polish <- function(y, fit, blocks = NULL, nb_inno_size_ub = 50,
                                max_iter = 25, tol = 1e-10) {
    N <- ncol(y)
    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    loglik <- .inad_loglik_fun(y, fit, blocks)
    in_space <- function(tab, v) {
        a <- v[tab$group %in% c("a1", "a2")]
        th <- v[tab$group == "theta"]
        kp <- v[tab$group == "kappa"]
        all(a >= 0) && (thinning != "binom" || all(a < 1)) &&
            all(th > 0) && all(kp > 0) && all(kp <= nb_inno_size_ub)
    }
    ll <- loglik(fit)
    for (it in seq_len(max_iter)) {
        li <- tryCatch(.louis_full_info(y, fit, blocks), error = function(e) NULL)
        if (is.null(li)) break
        tab <- li$fp$tab
        est <- tab$est
        # Coefficients exactly on a bound where the thinned count is degenerate
        # (alpha = 0, or binomial alpha = 1) have a one-sided score. If it
        # points into the parameter space, the coefficient is moved off the
        # bound by a line search on the observed log-likelihood, and Newton
        # steps continue from the interior on the next iteration. Coefficients
        # that the search cannot move are held on the bound.
        is_a <- tab$group %in% c("a1", "a2")
        on_bnd <- is_a & vapply(est, .thin_at_bound, logical(1), thinning = thinning)
        inward <- on_bnd & ((est == 0 & li$score > 0) | (est == 1 & li$score < 0))
        stuck <- rep(FALSE, length(est))
        if (any(inward)) {
            moved <- FALSE
            for (k in which(inward)) {
                lim <- if (est[k] == 1) 0 else if (thinning == "binom") 1 - 1e-10 else est[k] + 1e3
                ikk <- li$info[k, k]
                t0 <- if (is.finite(ikk) && ikk > 0) abs(li$score[k]) / ikk else 1e-4
                cs <- .inad_coord_search(loglik, fit, li$fp, k, lim, t0, N, ll0 = ll)
                if (cs$gain > 0) {
                    fit <- cs$fit
                    ll <- ll + cs$gain
                    # Later searches rebuild the fit from these values.
                    li$fp$tab$est[k] <- cs$par
                    moved <- TRUE
                } else {
                    stuck[k] <- TRUE
                }
            }
            if (moved) next
        }
        # Active set: parameters on the boundary whose score points out of the
        # parameter space stay fixed (Karush-Kuhn-Tucker conditions). A zero
        # score at zero also fixes a coefficient that is not identifiable.
        at_lower <- (tab$group %in% c("a1", "a2", "theta", "kappa") & est <= 1e-8 & li$score <= 0) |
            stuck | (on_bnd & est == 0)
        at_upper <- (tab$group == "kappa" & est >= nb_inno_size_ub * (1 - 1e-8) & li$score > 0) |
            (thinning == "binom" & is_a & est >= 1 - 1e-8 & li$score > 0) | (on_bnd & est == 1)
        free <- !(at_lower | at_upper)
        if (!any(free)) break
        Iff <- li$info[free, free, drop = FALSE]
        sf <- li$score[free]
        try_step <- function(step, lam_min) {
            # Step halving; returns the accepted candidate or NULL.
            lam <- 1
            while (lam >= lam_min) {
                v <- est + lam * step
                if (in_space(li$fp$tab, v)) {
                    cand <- .inad_set_free(fit, li$fp, v, N)
                    ll_c <- loglik(cand)
                    if (is.finite(ll_c) && ll_c >= ll - 1e-9) {
                        return(list(fit = cand, ll = ll_c, change = max(abs(lam * step))))
                    }
                }
                lam <- lam / 2
            }
            NULL
        }
        lower_grp <- tab$group %in% c("a1", "a2", "theta", "kappa")
        step <- numeric(length(est))
        sol <- tryCatch(solve(Iff, sf), error = function(e) NULL)
        # A parameter already on its lower bound whose Newton step points out of
        # the parameter space is held there and the step recomputed.
        for (k in seq_len(5)) {
            if (is.null(sol) || any(!is.finite(sol))) break
            step[] <- 0
            step[free] <- sol
            hold <- free & lower_grp & est <= 1e-8 & step < 0
            if (!any(hold)) break
            free <- free & !hold
            if (!any(free)) {
                sol <- NULL
                break
            }
            Iff <- li$info[free, free, drop = FALSE]
            sf <- li$score[free]
            sol <- tryCatch(solve(Iff, sf), error = function(e) NULL)
        }
        if (!any(free)) break
        acc <- NULL
        if (!is.null(sol) && all(is.finite(sol))) {
            step[] <- 0
            step[free] <- sol
            # Projection: parameters that the full step would take below their
            # lower bound are moved onto the bound, and the move is kept if the
            # log-likelihood does not decrease; the active-set rule then holds
            # them there. Without this, a coefficient slightly above zero whose
            # score points outward makes the step-halving search stall.
            cross <- free & lower_grp & est > 1e-8 & est + step < 0
            if (any(cross)) {
                v <- est
                v[cross] <- ifelse(tab$group[cross] %in% c("a1", "a2"), 0, 1e-10)
                cand <- .inad_set_free(fit, li$fp, v, N)
                ll_c <- loglik(cand)
                if (is.finite(ll_c) && ll_c >= ll - 1e-9) {
                    fit <- cand
                    ll <- max(ll, ll_c)
                    next
                }
            }
            acc <- try_step(step, 1 / 1024)
            if (!is.null(acc)) {
                fit <- acc$fit
                ll <- max(ll, acc$ll)
                if (acc$change < tol) break
                next
            }
        }
        # Damped (Levenberg-Marquardt) steps when the Newton step fails, as far
        # from the maximum, where the observed information need not be positive
        # definite: a large damping gives a short step along the score, which
        # increases the log-likelihood.
        d <- pmax(abs(diag(Iff)), 1e-8)
        mu <- 1e-3
        while (is.null(acc) && mu <= 1e10) {
            sol <- tryCatch(solve(Iff + mu * diag(d, nrow = length(d)), sf), error = function(e) NULL)
            if (!is.null(sol) && all(is.finite(sol))) {
                step <- numeric(length(est))
                step[free] <- sol
                acc <- try_step(step, 1 / 8)
                if (!is.null(acc) && acc$ll <= ll + 1e-12) acc <- NULL
            }
            mu <- mu * 10
        }
        if (is.null(acc)) break
        fit <- acc$fit
        ll <- max(ll, acc$ll)
    }
    list(fit = fit, log_l = ll)
}
