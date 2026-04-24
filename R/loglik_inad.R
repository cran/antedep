#' Log-likelihood for INAD models (with missing data support)
#'
#' If blocks is NULL, this computes the log likelihood as the sum of per time
#' contributions from logL_inad_i for computational convenience.
#'
#' @param y Integer matrix n_sub by n_time.
#' @param order Integer in \{0, 1, 2\}.
#' @param thinning One of "binom", "pois", "nbinom".
#' @param innovation One of "pois", "bell", "nbinom".
#' @param alpha Thinning parameters. For order 1, numeric length 1 or n_time.
#'   For order 2, either a matrix n_time by 2 or a list(alpha1, alpha2).
#' @param theta Innovation parameter(s). Numeric length 1 or n_time. For
#'   Poisson and negative binomial innovations, this is the innovation mean.
#'   For Bell innovations, this is the Bell rate parameter
#'   (mean \eqn{\theta e^\theta}).
#' @param nb_inno_size Size parameter for innovation "nbinom". Numeric length 1 or n_time.
#' @param blocks Optional integer vector of length n_sub. If NULL, no fixed effect.
#' @param tau Optional numeric vector. Only used if blocks is not NULL.
#' @param na_action How to handle missing values:
#'   \itemize{
#'     \item \code{"fail"}: error if any NA is present.
#'     \item \code{"complete"}: use only complete-case subjects.
#'     \item \code{"marginalize"}: observed-data likelihood under MAR via truncated-state recursion.
#'   }
#'
#' @return A scalar log likelihood.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 60,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 20)
#' logL_inad(
#'   y,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = fit$alpha,
#'   theta = fit$theta
#' )
#' @export
logL_inad <- function(
        y,
        order = 1,
        thinning = c("binom", "pois", "nbinom"),
        innovation = c("pois", "bell", "nbinom"),
        alpha,
        theta,
        nb_inno_size = NULL,
        blocks = NULL,
        tau = 0,
        na_action = c("fail", "complete", "marginalize")
) {
    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)
    na_action <- match.arg(na_action)

    if (!is.matrix(y)) y <- as.matrix(y)
    y_obs <- y[!is.na(y)]
    if (any(!is.finite(y_obs)) || any(y_obs < 0) || any(y_obs != floor(y_obs))) return(-Inf)

    n <- nrow(y)
    N <- ncol(y)
    if (!(order %in% c(0, 1, 2))) stop("order must be 0, 1, or 2")
    has_missing <- any(is.na(y))

    if (has_missing) {
        if (na_action == "fail") {
            stop("y contains NA values. Use na_action = 'complete' or 'marginalize'.", call. = FALSE)
        }

        if (na_action == "complete") {
            cc <- .extract_complete_cases(y, blocks, warn = FALSE)
            y <- cc$y
            blocks <- cc$blocks
            n <- nrow(y)
            has_missing <- FALSE
        }

        if (na_action == "marginalize" && has_missing) {
            return(.logL_inad_missing(
                y = y,
                order = order,
                thinning = thinning,
                innovation = innovation,
                alpha = alpha,
                theta = theta,
                nb_inno_size = nb_inno_size,
                blocks = blocks,
                tau = tau
            ))
        }
    }

    if (is.null(blocks)) {
        return(sum(vapply(
            seq_len(N),
            function(i) logL_inad_i(
                y = y,
                i = i,
                order = order,
                thinning = thinning,
                innovation = innovation,
                alpha = alpha,
                theta = theta,
                nb_inno_size = nb_inno_size
            ),
            numeric(1)
        )))
    }

    if (length(blocks) != n) stop("length(blocks) must equal nrow(y)")
    blocks <- as.integer(blocks)

    B <- max(blocks)
    tau <- as.numeric(tau)

    if (B <= 1L) {
        tau <- 0
    } else if (length(tau) == 1L) {
        tau <- c(0, rep(tau, B - 1L))
    } else {
        if (length(tau) < B) stop("tau must be length 1 or at least max(blocks)")
        tau <- tau[seq_len(B)]
        tau[1] <- 0
    }

    theta <- as.numeric(theta)
    if (length(theta) == 1L) theta <- rep(theta, N)
    if (length(theta) != N) stop("theta must be length 1 or ncol(y)")

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) stop("nb_inno_size is required for innovation='nbinom'")
        nb_inno_size <- as.numeric(nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("nb_inno_size must be length 1 or ncol(y)")
        if (any(nb_inno_size <= 0)) return(-Inf)
    }

    ap <- .unpack_alpha(order, alpha, N)
    alpha1 <- ap$a1
    alpha2 <- ap$a2

    loglik <- 0

    lam_eff <- matrix(NA_real_, nrow = B, ncol = N)
    for (b in seq_len(B)) {
        lam_b <- .inad_effective_innovation_param(theta, tau[b], innovation)
        if (any(!is.finite(lam_b)) || any(lam_b <= 0)) return(-Inf)
        lam_eff[b, ] <- lam_b
    }

    for (s in 1:n) {
        b <- blocks[s]
        for (i in 1:N) {
            lam <- lam_eff[b, i]
            if (!is.finite(lam) || lam <= 0) return(-Inf)

            y_si <- y[s, i]

            if (i == 1 || order == 0) {
                p0 <- .innov_vec(y_si, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)
                if (!is.finite(p0) || p0 <= 0) return(-Inf)
                loglik <- loglik + log(p0)
                next
            }

            if (order == 1 || i == 2) {
                y_prev <- y[s, i - 1]
                k_vals <- 0:y_si

                thin_vals <- .thin_vec(k_vals, y_prev, alpha1[i], thinning)
                innov_vals <- .innov_vec(y_si - k_vals, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

                conv <- sum(thin_vals * innov_vals)
                if (!is.finite(conv) || conv <= 0) return(-Inf)

                loglik <- loglik + log(conv)
                next
            }

            # Order 2, i >= 3
            y1 <- y[s, i - 1]
            y2 <- y[s, i - 2]

            k1_vals <- 0:y_si
            log_terms <- rep(-Inf, length(k1_vals))

            for (idx in seq_along(k1_vals)) {
                k1 <- k1_vals[idx]
                thin1 <- .thin_vec(k1, y1, alpha1[i], thinning)
                if (!is.finite(thin1) || thin1 <= 0) next

                remain <- y_si - k1
                k2_vals <- 0:remain

                thin2_vals <- .thin_vec(k2_vals, y2, alpha2[i], thinning)
                innov_vals <- .innov_vec(remain - k2_vals, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

                inner_sum <- sum(thin2_vals * innov_vals)
                if (!is.finite(inner_sum) || inner_sum <= 0) next

                log_terms[idx] <- log(thin1) + log(inner_sum)
            }

            log_conv <- .log_sum_exp(log_terms)
            if (!is.finite(log_conv)) return(-Inf)
            loglik <- loglik + log_conv
        }
    }

    loglik
}

#' INAD log likelihood contribution at time i (no fixed effect)
#'
#' Returns the time i contribution, summed over subjects, for the no fixed effect model.
#'
#' @param y Integer matrix n_sub by n_time.
#' @param i Time index in 1..ncol(y).
#' @param order Integer in \{0, 1, 2\}.
#' @param thinning One of "binom", "pois", "nbinom".
#' @param innovation One of "pois", "bell", "nbinom".
#' @param alpha Thinning parameters. For order 1, numeric length 1 or n_time.
#'   For order 2, either a matrix n_time by 2 or a list(alpha1, alpha2).
#' @param theta Innovation parameter at time i, or a vector length 1 or n_time.
#'   For Poisson and negative binomial innovations, this is the innovation
#'   mean. For Bell innovations, this is the Bell rate parameter
#'   (mean \eqn{\theta e^\theta}).
#' @param nb_inno_size Size parameter for innovation "nbinom". Numeric length 1 or n_time.
#'
#' @return A scalar log likelihood contribution for time i.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 50,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 20)
#' logL_inad_i(
#'   y = y,
#'   i = 3,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = fit$alpha,
#'   theta = fit$theta
#' )
#' @export
logL_inad_i <- function(
        y,
        i,
        order = 1,
        thinning = c("binom", "pois", "nbinom"),
        innovation = c("pois", "bell", "nbinom"),
        alpha,
        theta,
        nb_inno_size = NULL
) {
    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)

    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) return(-Inf)

    n <- nrow(y)
    N <- ncol(y)

    if (!(order %in% c(0, 1, 2))) stop("order must be 0, 1, or 2")
    if (length(i) != 1L || i < 1L || i > N) stop("i must be in 1..ncol(y)")

    theta <- as.numeric(theta)
    if (length(theta) == 1L) theta <- rep(theta, N)
    if (length(theta) != N) stop("theta must be length 1 or ncol(y)")

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) stop("nb_inno_size is required for innovation='nbinom'")
        nb_inno_size <- as.numeric(nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("nb_inno_size must be length 1 or ncol(y)")
        if (any(nb_inno_size <= 0)) return(-Inf)
    }

    ap <- .unpack_alpha(order, alpha, N)
    alpha1 <- ap$a1
    alpha2 <- ap$a2

    loglik_i <- 0
    th_i <- theta[i]
    if (!is.finite(th_i) || th_i <= 0) return(-Inf)

    for (s in 1:n) {
        y_si <- y[s, i]

        if (i == 1 || order == 0) {
            p0 <- .innov_vec(y_si, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)
            if (!is.finite(p0) || p0 <= 0) return(-Inf)
            loglik_i <- loglik_i + log(p0)
            next
        }

        if (order == 1 || i == 2) {
            y_prev <- y[s, i - 1]
            k_vals <- 0:y_si

            thin_vals <- .thin_vec(k_vals, y_prev, alpha1[i], thinning)
            innov_vals <- .innov_vec(y_si - k_vals, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

            conv <- sum(thin_vals * innov_vals)
            if (!is.finite(conv) || conv <= 0) return(-Inf)

            loglik_i <- loglik_i + log(conv)
            next
        }

        # Order 2, i >= 3
        y1 <- y[s, i - 1]
        y2 <- y[s, i - 2]

        k1_vals <- 0:y_si
        log_terms <- rep(-Inf, length(k1_vals))

        for (idx in seq_along(k1_vals)) {
            k1 <- k1_vals[idx]
            thin1 <- .thin_vec(k1, y1, alpha1[i], thinning)
            if (!is.finite(thin1) || thin1 <= 0) next

            remain <- y_si - k1
            k2_vals <- 0:remain

            thin2_vals <- .thin_vec(k2_vals, y2, alpha2[i], thinning)
            innov_vals <- .innov_vec(remain - k2_vals, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

            inner_sum <- sum(thin2_vals * innov_vals)
            if (!is.finite(inner_sum) || inner_sum <= 0) next

            log_terms[idx] <- log(thin1) + log(inner_sum)
        }

        log_conv <- .log_sum_exp(log_terms)
        if (!is.finite(log_conv)) return(-Inf)
        loglik_i <- loglik_i + log_conv
    }

    loglik_i
}

#' Observed-data INAD likelihood with missing values under MAR
#'
#' Uses forward recursion over a truncated count state space.
#'
#' @keywords internal
.logL_inad_missing <- function(
        y,
        order,
        thinning,
        innovation,
        alpha,
        theta,
        nb_inno_size = NULL,
        blocks = NULL,
        tau = 0
) {
    n <- nrow(y)
    N <- ncol(y)

    theta <- as.numeric(theta)
    if (length(theta) == 1L) theta <- rep(theta, N)
    if (length(theta) != N) stop("theta must be length 1 or ncol(y)")

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) stop("nb_inno_size is required for innovation='nbinom'")
        nb_inno_size <- as.numeric(nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("nb_inno_size must be length 1 or ncol(y)")
        if (any(nb_inno_size <= 0) || any(!is.finite(nb_inno_size))) return(-Inf)
    }

    ap <- .unpack_alpha(order, alpha, N)
    alpha1 <- ap$a1
    alpha2 <- ap$a2

    if (is.null(blocks)) {
        blocks <- rep(1L, n)
        tau_use <- 0
        n_blocks <- 1L
    } else {
        if (length(blocks) != n) stop("length(blocks) must equal nrow(y)")
        blocks <- as.integer(blocks)
        if (any(!is.finite(blocks)) || any(blocks < 1L)) stop("blocks must be positive integers")

        n_blocks <- max(blocks)
        tau_use <- as.numeric(tau)
        if (n_blocks <= 1L) {
            tau_use <- 0
        } else if (length(tau_use) == 1L) {
            tau_use <- c(0, rep(tau_use, n_blocks - 1L))
        } else {
            if (length(tau_use) < n_blocks) stop("tau must be length 1 or at least max(blocks)")
            tau_use <- tau_use[seq_len(n_blocks)]
            tau_use[1] <- 0
        }
    }

    K <- .inad_state_max(
        y = y,
        order = order,
        alpha1 = alpha1,
        alpha2 = alpha2,
        theta = theta,
        tau = tau_use,
        innovation = innovation
    )

    trans_by_block <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
        tau_b <- if (length(tau_use) == 1L) tau_use else tau_use[b]
        lam <- .inad_effective_innovation_param(theta, tau_b, innovation)
        if (any(!is.finite(lam)) || any(lam <= 0)) return(-Inf)
        trans_by_block[[b]] <- .inad_build_transitions(
            order = order,
            N = N,
            K = K,
            alpha1 = alpha1,
            alpha2 = alpha2,
            thinning = thinning,
            innovation = innovation,
            lambda = lam,
            nb_inno_size = nb_inno_size
        )
    }

    ll <- 0
    for (s in seq_len(n)) {
        y_s <- y[s, ]
        if (all(is.na(y_s))) next

        tr <- trans_by_block[[blocks[s]]]
        ll_s <- if (order == 0) {
            .inad_subject_ll_order0(y_s, tr$innov_list, K)
        } else if (order == 1) {
            .inad_subject_ll_order1(y_s, tr$innov_list[[1]], tr$trans1, K)
        } else {
            .inad_subject_ll_order2(y_s, tr$innov_list[[1]], tr$trans1[[2]], tr$order2_info, K)
        }

        if (!is.finite(ll_s)) return(-Inf)
        ll <- ll + ll_s
    }

    ll
}

#' Build transition objects for truncated-state missing-data likelihood
#' @keywords internal
.inad_build_transitions <- function(
        order, N, K, alpha1, alpha2, thinning, innovation, lambda, nb_inno_size
) {
    innov_list <- vector("list", N)
    for (t in seq_len(N)) {
        sz_t <- if (innovation == "nbinom") nb_inno_size[t] else NA_real_
        innov_list[[t]] <- .innov_vec(0:K, lambda[t], innovation, sz_t)
    }

    if (order == 0) {
        return(list(innov_list = innov_list, trans1 = NULL, order2_info = NULL))
    }

    trans1 <- vector("list", N)
    order2_info <- vector("list", N)

    for (t in 2:N) {
        thin1 <- lapply(0:K, function(prev) .inad_make_thin_pmf(prev, alpha1[t], thinning, K))

        if (order == 1 || t == 2) {
            M <- matrix(0, nrow = K + 1L, ncol = K + 1L)
            for (prev in 0:K) {
                M[prev + 1L, ] <- .inad_conv_trunc(thin1[[prev + 1L]], innov_list[[t]], K)
            }
            trans1[[t]] <- M
        } else {
            thin2 <- lapply(0:K, function(prev) .inad_make_thin_pmf(prev, alpha2[t], thinning, K))
            order2_info[[t]] <- list(
                thin1 = thin1,
                thin2 = thin2,
                innov = innov_list[[t]],
                cache = new.env(parent = emptyenv())
            )
        }
    }

    list(innov_list = innov_list, trans1 = trans1, order2_info = order2_info)
}

#' Subject-level observed-data likelihood for INAD(0)
#' @keywords internal
.inad_subject_ll_order0 <- function(y_s, innov_list, K) {
    obs_idx <- which(!is.na(y_s))
    if (length(obs_idx) == 0) return(0)

    ll <- 0
    for (t in obs_idx) {
        y_t <- y_s[t]
        if (y_t > K) return(-Inf)
        p <- innov_list[[t]][y_t + 1L]
        if (!is.finite(p) || p <= 0) return(-Inf)
        ll <- ll + log(p)
    }
    ll
}

#' Subject-level observed-data likelihood for INAD(1)
#' @keywords internal
.inad_subject_ll_order1 <- function(y_s, init_pmf, trans1, K) {
    N <- length(y_s)
    if (N == 0) return(0)

    y1 <- y_s[1]
    if (is.na(y1)) {
        p_prev <- init_pmf
    } else {
        if (y1 > K) return(-Inf)
        p_prev <- rep(0, K + 1L)
        p_prev[y1 + 1L] <- init_pmf[y1 + 1L]
    }

    z0 <- sum(p_prev)
    if (!is.finite(z0) || z0 <= 0) return(-Inf)
    ll <- log(z0)
    p_prev <- p_prev / z0

    for (t in 2:N) {
        M <- trans1[[t]]
        if (is.null(M)) return(-Inf)

        y_t <- y_s[t]
        p_next <- rep(0, K + 1L)

        if (is.na(y_t)) {
            p_next <- as.numeric(p_prev %*% M)
        } else {
            if (y_t > K) return(-Inf)
            v <- sum(p_prev * M[, y_t + 1L])
            p_next[y_t + 1L] <- v
        }

        z <- sum(p_next)
        if (!is.finite(z) || z <= 0) return(-Inf)
        ll <- ll + log(z)
        p_prev <- p_next / z
    }

    ll
}

#' Subject-level observed-data likelihood for INAD(2)
#' @keywords internal
.inad_subject_ll_order2 <- function(y_s, init_pmf, trans_t2, order2_info, K) {
    N <- length(y_s)
    if (N == 0) return(0)
    if (N == 1) {
        y1 <- y_s[1]
        if (is.na(y1)) return(log(sum(init_pmf)))
        if (y1 > K) return(-Inf)
        p <- init_pmf[y1 + 1L]
        if (!is.finite(p) || p <= 0) return(-Inf)
        return(log(p))
    }

    y1 <- y_s[1]
    if (is.na(y1)) {
        p_y1 <- init_pmf
    } else {
        if (y1 > K) return(-Inf)
        p_y1 <- rep(0, K + 1L)
        p_y1[y1 + 1L] <- init_pmf[y1 + 1L]
    }
    z1 <- sum(p_y1)
    if (!is.finite(z1) || z1 <= 0) return(-Inf)
    ll <- log(z1)
    p_y1 <- p_y1 / z1

    y2 <- y_s[2]
    pair <- matrix(0, nrow = K + 1L, ncol = K + 1L)
    for (prev2 in 0:K) {
        w <- p_y1[prev2 + 1L]
        if (w <= 0) next
        if (is.na(y2)) {
            pair[prev2 + 1L, ] <- pair[prev2 + 1L, ] + w * trans_t2[prev2 + 1L, ]
        } else {
            if (y2 > K) return(-Inf)
            pair[prev2 + 1L, y2 + 1L] <- pair[prev2 + 1L, y2 + 1L] + w * trans_t2[prev2 + 1L, y2 + 1L]
        }
    }

    z <- sum(pair)
    if (!is.finite(z) || z <= 0) return(-Inf)
    ll <- ll + log(z)
    pair <- pair / z

    if (N == 2) return(ll)

    for (t in 3:N) {
        info <- order2_info[[t]]
        if (is.null(info)) return(-Inf)
        y_t <- y_s[t]

        next_pair <- matrix(0, nrow = K + 1L, ncol = K + 1L)
        nz <- which(pair > 0, arr.ind = TRUE)

        if (is.na(y_t)) {
            for (k in seq_len(nrow(nz))) {
                prev2 <- nz[k, 1] - 1L
                prev1 <- nz[k, 2] - 1L
                w <- pair[nz[k, 1], nz[k, 2]]
                dist <- .inad_order2_dist(prev2, prev1, info, K)
                next_pair[prev1 + 1L, ] <- next_pair[prev1 + 1L, ] + w * dist
            }
        } else {
            if (y_t > K) return(-Inf)
            col_idx <- y_t + 1L
            for (k in seq_len(nrow(nz))) {
                prev2 <- nz[k, 1] - 1L
                prev1 <- nz[k, 2] - 1L
                w <- pair[nz[k, 1], nz[k, 2]]
                dist <- .inad_order2_dist(prev2, prev1, info, K)
                next_pair[prev1 + 1L, col_idx] <- next_pair[prev1 + 1L, col_idx] + w * dist[col_idx]
            }
        }

        z <- sum(next_pair)
        if (!is.finite(z) || z <= 0) return(-Inf)
        ll <- ll + log(z)
        pair <- next_pair / z
    }

    ll
}

#' Retrieve/calculate order-2 transition distribution for one previous state pair
#' @keywords internal
.inad_order2_dist <- function(prev2, prev1, info, K) {
    key <- paste0(prev2, "_", prev1)
    if (exists(key, envir = info$cache, inherits = FALSE)) {
        return(get(key, envir = info$cache, inherits = FALSE))
    }

    p12 <- .inad_conv_trunc(info$thin1[[prev1 + 1L]], info$thin2[[prev2 + 1L]], K)
    dist <- .inad_conv_trunc(p12, info$innov, K)
    assign(key, dist, envir = info$cache)
    dist
}

#' Truncated discrete convolution
#' @keywords internal
.inad_conv_trunc <- function(p, q, K) {
    out <- numeric(K + 1L)
    idx_p <- which(p > 0) - 1L
    for (i in idx_p) {
        pi <- p[i + 1L]
        max_j <- min(K - i, K)
        if (max_j < 0L) next
        out[(i + 0:max_j) + 1L] <- out[(i + 0:max_j) + 1L] + pi * q[1:(max_j + 1L)]
    }
    out
}

#' Build truncated thinning pmf for one previous count
#' @keywords internal
.inad_make_thin_pmf <- function(prev, alpha, thinning, K) {
    out <- numeric(K + 1L)
    k_max <- min(prev, K)
    if (k_max >= 0L) {
        vals <- .thin_vec(0:k_max, prev, alpha, thinning)
        out[1:(k_max + 1L)] <- vals
    }
    out
}

#' Heuristic state-space bound for INAD missing-data recursion
#' @keywords internal
.inad_state_max <- function(y, order, alpha1, alpha2, theta, tau, innovation = "pois") {
    y_obs <- y[!is.na(y)]
    obs_max <- if (length(y_obs) == 0L) 0 else max(y_obs)

    if (length(tau) == 1L) {
        mu_eff <- .inad_effective_innovation_mean(theta, tau, innovation = innovation)
    } else {
        mu_eff <- as.numeric(outer(
            tau,
            theta,
            FUN = function(tb, th) .inad_effective_innovation_mean(th, tb, innovation = innovation)
        ))
    }
    lam_max <- max(mu_eff, na.rm = TRUE)

    a_sum <- 0
    if (order == 1) {
        if (!is.null(alpha1) && length(alpha1) >= 2L) a_sum <- max(alpha1[2:length(alpha1)], na.rm = TRUE)
    }
    if (order == 2) {
        a2 <- if (is.null(alpha2)) rep(0, length(alpha1)) else alpha2
        a_sum <- max(alpha1 + a2, na.rm = TRUE)
    }

    if (!is.finite(a_sum)) a_sum <- 0
    a_sum <- min(max(a_sum, 0), 0.98)
    gain <- 1 / (1 - a_sum)

    mu_upper <- max(obs_max, lam_max * gain)
    k <- ceiling(mu_upper + 8 * sqrt(mu_upper + 1) + 10)
    k <- max(k, obs_max + 10, 20)
    k <- min(k, 300L)
    if (obs_max > k) k <- obs_max
    as.integer(k)
}
