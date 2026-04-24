#' Fit INAD antedependence model by maximum likelihood
#'
#' Fits INAD models by maximum likelihood.
#'
#' No fixed effect: time separable optimization using logL_inad_i with theta eliminated
#' by moment equations for order 1 and 2.
#'
#' Fixed effect: block coordinate descent using nloptr BOBYQA, updating tau, alpha,
#' theta, and nb_inno_size if needed.
#'
#' @param y Integer matrix n_sub by n_time.
#' @param order Integer in \{0, 1, 2\}.
#' @param thinning One of "binom", "pois", "nbinom".
#' @param innovation One of "pois", "bell", "nbinom".
#' @param blocks Optional integer vector length n_sub. Default NULL.
#' @param max_iter Max iterations for FE coordinate descent.
#' @param tol Tolerance for FE log likelihood stopping.
#' @param verbose Logical.
#' @param init_alpha Optional initial alpha. For order 1 numeric length 1 or n_time.
#'   For order 2 matrix n_time by 2 or list(alpha1, alpha2).
#' @param init_theta Optional initial theta numeric length 1 or n_time.
#' @param init_tau Optional initial tau. Scalar expands to c(0, x, ..., x). Vector forces first to 0.
#' @param init_nb_inno_size Optional initial size for innovation nbinom, length 1 or n_time.
#' @param nb_inno_size_ub Upper bound for innovation negative binomial size
#'   parameter when \code{innovation = "nbinom"}. Default is 50.
#' @param na_action How to handle missing values:
#'   \itemize{
#'     \item \code{"fail"}: stop if any NA is present.
#'     \item \code{"complete"}: fit using complete-case subjects only.
#'     \item \code{"marginalize"}: maximize observed-data likelihood under MAR.
#'   }
#'
#' @return A list of class \code{"inad_fit"} containing:
#' \describe{
#'   \item{alpha}{Estimated antedependence parameter(s)}
#'   \item{theta}{Estimated innovation parameter(s)}
#'   \item{tau}{Estimated block effects (if applicable)}
#'   \item{nb_inno_size}{Estimated innovation NB size parameter(s), when
#'     \code{innovation = "nbinom"}}
#'   \item{log_l}{Maximized log-likelihood}
#'   \item{aic}{Akaike information criterion}
#'   \item{bic}{Bayesian information criterion}
#'   \item{n_params}{Number of free parameters}
#'   \item{convergence}{Convergence code}
#'   \item{settings}{Model and fitting settings}
#' }
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
#' fit$log_l
#' @export
fit_inad <- function(
        y,
        order = 1,
        thinning = c("binom", "pois", "nbinom"),
        innovation = c("pois", "bell", "nbinom"),
        blocks = NULL,
        max_iter = 50,
        tol = 1e-6,
        verbose = FALSE,
        init_alpha = NULL,
        init_theta = NULL,
        init_tau = 0.4,
        init_nb_inno_size = 1,
        nb_inno_size_ub = 50,
        na_action = c("fail", "complete", "marginalize")
) {
    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)
    na_action <- match.arg(na_action)

    if (!is.matrix(y)) y <- as.matrix(y)
    y_obs <- y[!is.na(y)]
    if (any(!is.finite(y_obs)) || any(y_obs < 0) || any(y_obs != floor(y_obs))) {
        stop("y must be a nonnegative integer matrix")
    }
    if (!(order %in% c(0, 1, 2))) stop("order must be 0, 1, or 2")
    if (!is.numeric(nb_inno_size_ub) || length(nb_inno_size_ub) != 1L ||
        !is.finite(nb_inno_size_ub) || nb_inno_size_ub <= 0) {
        stop("nb_inno_size_ub must be a positive finite scalar")
    }
    has_missing <- any(is.na(y))

    if (has_missing) {
        if (na_action == "fail") {
            stop("y contains NA values. Use na_action = 'complete' or 'marginalize'.", call. = FALSE)
        }

        if (na_action == "complete") {
            cc <- .extract_complete_cases(y, blocks, warn = TRUE)
            y <- cc$y
            blocks <- cc$blocks
            has_missing <- FALSE
        }

        if (na_action == "marginalize" && has_missing) {
            block_info <- .normalize_blocks(blocks, nrow(y))
            blocks_id <- if (is.null(blocks)) NULL else block_info$blocks_id
            fit <- .fit_inad_missing(
                y = y,
                order = order,
                thinning = thinning,
                innovation = innovation,
                blocks = blocks_id,
                max_iter = max_iter,
                tol = tol,
                verbose = verbose,
                init_alpha = init_alpha,
                init_theta = init_theta,
                init_tau = init_tau,
                init_nb_inno_size = init_nb_inno_size,
                nb_inno_size_ub = nb_inno_size_ub
            )
            return(.finalize_inad_fit(
                fit, y = y, na_action = na_action,
                blocks_id = blocks_id,
                block_levels = if (is.null(blocks)) NULL else block_info$block_levels
            ))
        }
    }

    if (is.null(blocks)) {
        fit <- fit_inad_no_fe(
            y = y,
            order = order,
            thinning = thinning,
            innovation = innovation,
            init_alpha = init_alpha,
            init_theta = init_theta,
            init_nb_inno_size = init_nb_inno_size,
            nb_inno_size_ub = nb_inno_size_ub
        )
        return(.finalize_inad_fit(fit, y = y, na_action = na_action, blocks_id = NULL, block_levels = NULL))
    }

    block_info <- .normalize_blocks(blocks, nrow(y))
    blocks_id <- block_info$blocks_id
    fit <- fit_inad_fe(
        y = y,
        order = order,
        thinning = thinning,
        innovation = innovation,
        blocks = blocks_id,
        max_iter = max_iter,
        tol = tol,
        verbose = verbose,
        init_alpha = init_alpha,
        init_theta = init_theta,
        init_tau = init_tau,
        init_nb_inno_size = init_nb_inno_size,
        nb_inno_size_ub = nb_inno_size_ub
    )
    .finalize_inad_fit(
        fit, y = y, na_action = na_action,
        blocks_id = blocks_id,
        block_levels = block_info$block_levels
    )
}

#' @keywords internal
.finalize_inad_fit <- function(fit, y, na_action = NULL, blocks_id = NULL, block_levels = NULL) {
    if (is.null(fit$settings)) fit$settings <- list()
    if (is.null(fit$settings$n_subjects)) fit$settings$n_subjects <- nrow(y)
    if (is.null(fit$settings$n_time)) fit$settings$n_time <- ncol(y)
    if (!is.null(na_action)) fit$settings$na_action <- na_action
    if (!is.null(blocks_id)) {
        n_blocks <- length(unique(blocks_id))
        fit$settings$blocks <- if (n_blocks > 1L) blocks_id else NULL
        fit$settings$n_blocks <- n_blocks
        fit$settings$block_levels <- block_levels
    }

    fit$n_obs <- sum(!is.na(y))
    fit$n_missing <- sum(is.na(y))
    fit$pct_missing <- mean(is.na(y)) * 100

    # Keep detailed convergence diagnostics while exposing a comparable code.
    if (is.list(fit$convergence)) {
        fit$convergence_info <- fit$convergence
        if (!is.null(fit$convergence$code)) {
            fit$convergence <- as.integer(fit$convergence$code)
        } else if (!is.null(fit$convergence$per_time)) {
            per_time <- as.integer(fit$convergence$per_time)
            if (length(per_time) == 0L || all(is.na(per_time))) {
                fit$convergence <- NA_integer_
            } else {
                fit$convergence <- as.integer(max(per_time, na.rm = TRUE))
            }
        } else {
            fit$convergence <- NA_integer_
        }
    } else if (is.null(fit$convergence)) {
        fit$convergence <- NA_integer_
    } else {
        fit$convergence <- as.integer(fit$convergence)
    }

    fit$n_params <- tryCatch(
        as.integer(.count_params_inad_fit(fit)),
        error = function(e) NA_integer_
    )

    if (is.finite(fit$log_l) && is.finite(fit$n_params)) {
        fit$aic <- -2 * fit$log_l + 2 * fit$n_params
        fit$bic <- -2 * fit$log_l + fit$n_params * log(fit$settings$n_subjects)
    } else {
        fit$aic <- NA_real_
        fit$bic <- NA_real_
    }

    class(fit) <- "inad_fit"
    fit
}

.fit_inad_missing <- function(
        y,
        order,
        thinning,
        innovation,
        blocks,
        max_iter,
        tol,
        verbose,
        init_alpha,
        init_theta,
        init_tau,
        init_nb_inno_size,
        nb_inno_size_ub
) {
    n <- nrow(y)
    N <- ncol(y)
    eps_pos <- 1e-8
    eps_alpha <- 1e-8
    alpha_ub <- 1 - eps_alpha

    # Normalize blocks and tau setup.
    if (!is.null(blocks)) {
        if (length(blocks) != n) stop("length(blocks) must equal nrow(y)")
        blocks <- as.integer(factor(blocks))
        B <- max(blocks)
    } else {
        B <- 1L
    }

    # Build a fallback complete-data initializer from complete cases when possible.
    fit_init <- NULL
    cc <- tryCatch(.extract_complete_cases(y, blocks, warn = FALSE), error = function(e) NULL)
    if (!is.null(cc) && nrow(cc$y) >= max(5L, order + 2L)) {
        fit_init <- tryCatch(
            fit_inad(
                y = cc$y,
                order = order,
                thinning = thinning,
                innovation = innovation,
                blocks = cc$blocks,
                max_iter = max_iter,
                tol = tol,
                verbose = FALSE,
                init_alpha = init_alpha,
                init_theta = init_theta,
                init_tau = init_tau,
                init_nb_inno_size = init_nb_inno_size,
                nb_inno_size_ub = nb_inno_size_ub,
                na_action = "fail"
            ),
            error = function(e) NULL
        )
    }

    theta_init <- if (!is.null(fit_init)) {
        as.numeric(fit_init$theta)
    } else if (!is.null(init_theta)) {
        th <- as.numeric(init_theta)
        if (length(th) == 1L) th <- rep(th, N)
        if (length(th) != N) stop("init_theta must be length 1 or ncol(y)")
        th
    } else {
        col_mean <- colMeans(y, na.rm = TRUE)
        col_mean[!is.finite(col_mean)] <- 1
        pmax(col_mean, 0.5)
    }
    theta_init <- pmax(theta_init, eps_pos)

    nb_init <- NULL
    if (innovation == "nbinom") {
        nb_init <- if (!is.null(fit_init) && !is.null(fit_init$nb_inno_size)) {
            as.numeric(fit_init$nb_inno_size)
        } else {
            nb <- as.numeric(init_nb_inno_size)
            if (length(nb) == 1L) nb <- rep(nb, N)
            if (length(nb) != N) stop("init_nb_inno_size must be length 1 or ncol(y)")
            nb
        }
        nb_init <- pmin(pmax(nb_init, eps_pos), nb_inno_size_ub)
    }

    alpha_init <- NULL
    if (order == 1) {
        if (!is.null(fit_init) && !is.null(fit_init$alpha)) {
            alpha_init <- as.numeric(fit_init$alpha)
        } else if (!is.null(init_alpha)) {
            a <- as.numeric(init_alpha)
            if (length(a) == 1L) a <- rep(a, N)
            if (length(a) != N) stop("init_alpha must be length 1 or ncol(y) for order=1")
            alpha_init <- a
        } else {
            alpha_init <- rep(0, N)
            if (N >= 2) alpha_init[2:N] <- 0.3
        }
        alpha_init[1] <- 0
        alpha_init <- pmin(pmax(alpha_init, 0), alpha_ub)
    }
    if (order == 2) {
        if (!is.null(fit_init) && !is.null(fit_init$alpha)) {
            alpha_init <- as.matrix(fit_init$alpha)
        } else if (!is.null(init_alpha)) {
            if (is.matrix(init_alpha) && ncol(init_alpha) >= 2) {
                alpha_init <- init_alpha
                if (nrow(alpha_init) == 1L) alpha_init <- alpha_init[rep.int(1L, N), , drop = FALSE]
                if (nrow(alpha_init) != N) stop("for order=2, init_alpha matrix must have nrow 1 or ncol(y)")
                alpha_init <- alpha_init[, 1:2, drop = FALSE]
            } else if (is.list(init_alpha) && length(init_alpha) >= 2) {
                a1 <- as.numeric(init_alpha[[1]])
                a2 <- as.numeric(init_alpha[[2]])
                if (length(a1) == 1L) a1 <- rep(a1, N)
                if (length(a2) == 1L) a2 <- rep(a2, N)
                if (length(a1) != N || length(a2) != N) {
                    stop("init_alpha list entries must be length 1 or ncol(y)")
                }
                alpha_init <- cbind(a1, a2)
            } else {
                stop("for order=2, init_alpha must be matrix with 2 columns or list(alpha1, alpha2)")
            }
        } else {
            alpha_init <- matrix(0, nrow = N, ncol = 2)
            if (N >= 2) alpha_init[2:N, 1] <- 0.3
            if (N >= 3) alpha_init[3:N, 2] <- 0.1
        }
        alpha_init[1, 1] <- 0
        alpha_init[1, 2] <- 0
        if (N >= 2) alpha_init[2, 2] <- 0
        alpha_init <- pmin(pmax(alpha_init, 0), alpha_ub)
    }

    tau_init <- NULL
    if (B > 1L) {
        if (!is.null(fit_init) && !is.null(fit_init$tau)) {
            tau_init <- as.numeric(fit_init$tau)
        } else {
            tau_init <- as.numeric(init_tau)
        }
        if (length(tau_init) == 1L) {
            tau_init <- c(0, rep(tau_init, B - 1L))
        } else {
            if (length(tau_init) < B) stop("init_tau must be length 1 or at least max(blocks)")
            tau_init <- tau_init[seq_len(B)]
            tau_init[1] <- 0
        }
    }

    pack_init <- function() {
        out <- c(log(theta_init))
        if (innovation == "nbinom") out <- c(out, log(nb_init))
        if (order == 1 && N >= 2) {
            out <- c(out, stats::qlogis(pmin(pmax(alpha_init[2:N] / alpha_ub, eps_alpha), 1 - eps_alpha)))
        }
        if (order == 2 && N >= 2) {
            out <- c(out, stats::qlogis(pmin(pmax(alpha_init[2:N, 1] / alpha_ub, eps_alpha), 1 - eps_alpha)))
            if (N >= 3) {
                out <- c(out, stats::qlogis(pmin(pmax(alpha_init[3:N, 2] / alpha_ub, eps_alpha), 1 - eps_alpha)))
            }
        }
        if (B > 1L) out <- c(out, tau_init[2:B])
        out
    }

    unpack_par <- function(par) {
        idx <- 1L
        theta <- exp(par[idx:(idx + N - 1L)])
        idx <- idx + N

        nb <- NULL
        if (innovation == "nbinom") {
            nb <- exp(par[idx:(idx + N - 1L)])
            idx <- idx + N
        }

        alpha <- NULL
        if (order == 1) {
            alpha <- rep(0, N)
            if (N >= 2) {
                a_free <- par[idx:(idx + N - 2L)]
                alpha[2:N] <- alpha_ub * stats::plogis(a_free)
                idx <- idx + N - 1L
            }
        }
        if (order == 2) {
            alpha <- matrix(0, nrow = N, ncol = 2)
            if (N >= 2) {
                a1_free <- par[idx:(idx + N - 2L)]
                alpha[2:N, 1] <- alpha_ub * stats::plogis(a1_free)
                idx <- idx + N - 1L
            }
            if (N >= 3) {
                a2_free <- par[idx:(idx + N - 3L)]
                alpha[3:N, 2] <- alpha_ub * stats::plogis(a2_free)
                idx <- idx + N - 2L
            }
        }

        tau <- 0
        if (B > 1L) {
            tau <- c(0, par[idx:(idx + B - 2L)])
        }

        list(
            theta = pmax(theta, eps_pos),
            nb_inno_size = if (is.null(nb)) NULL else pmin(pmax(nb, eps_pos), nb_inno_size_ub),
            alpha = alpha,
            tau = tau
        )
    }

    obj <- function(par) {
        cur <- unpack_par(par)
        ll <- tryCatch(
            logL_inad(
                y = y,
                order = order,
                thinning = thinning,
                innovation = innovation,
                alpha = cur$alpha,
                theta = cur$theta,
                nb_inno_size = cur$nb_inno_size,
                blocks = blocks,
                tau = cur$tau,
                na_action = "marginalize"
            ),
            error = function(e) -Inf
        )
        if (!is.finite(ll)) return(1e10)
        -ll
    }

    par0 <- pack_init()
    control <- list(maxit = max(200L, 20L * max_iter), reltol = tol)
    opt <- optim(par = par0, fn = obj, method = "BFGS", control = control)
    if (!is.finite(opt$value)) {
        opt <- optim(par = par0, fn = obj, method = "Nelder-Mead", control = control)
    }

    cur <- unpack_par(opt$par)
    log_l <- logL_inad(
        y = y,
        order = order,
        thinning = thinning,
        innovation = innovation,
        alpha = cur$alpha,
        theta = cur$theta,
        nb_inno_size = cur$nb_inno_size,
        blocks = blocks,
        tau = cur$tau,
        na_action = "marginalize"
    )

    if (verbose) {
        message("Missing-data fit completed: logL=", round(log_l, 4), " convergence=", opt$convergence)
    }

    list(
        alpha = cur$alpha,
        theta = cur$theta,
        nb_inno_size = cur$nb_inno_size,
        tau = if (is.null(blocks)) NULL else cur$tau,
        log_l = log_l,
        loglik_i = NULL,
        n_obs = sum(!is.na(y)),
        n_missing = sum(is.na(y)),
        pct_missing = mean(is.na(y)) * 100,
        convergence = list(code = opt$convergence, message = opt$message, method = opt$method),
        settings = list(
            order = order,
            thinning = thinning,
            innovation = innovation,
            blocks = blocks,
            n_subjects = n,
            n_time = N,
            na_action = "marginalize"
        )
    )
}

fit_inad_no_fe <- function(
        y,
        order,
        thinning,
        innovation,
        init_alpha,
        init_theta,
        init_nb_inno_size,
        nb_inno_size_ub
) {
    N <- ncol(y)
    eps_alpha <- 1e-10
    eps_pos <- 1e-10
    alpha_ub <- 1 - eps_alpha

    mean_y <- colMeans(y)

    theta_from_inno_mean <- function(mu, innovation, eps_pos) {
        if (!is.finite(mu) || mu <= eps_pos) return(NA_real_)
        if (innovation == "bell") return(solve_theta_bell(mu))
        mu
    }

    default_theta0 <- function(innovation) {
        if (innovation == "bell") return(1)
        10
    }

    theta <- if (is.null(init_theta)) {
        rep(default_theta0(innovation), N)
    } else {
        th <- as.numeric(init_theta)
        if (length(th) == 1L) th <- rep(th, N)
        if (length(th) != N) stop("init_theta must be length 1 or ncol(y)")
        th
    }

    nb_inno_size <- NULL
    if (innovation == "nbinom") {
        nb_inno_size <- as.numeric(init_nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("init_nb_inno_size must be length 1 or ncol(y)")
        nb_inno_size <- pmin(pmax(nb_inno_size, eps_pos), nb_inno_size_ub)
    }

    init_alpha_order1 <- function(init_alpha, N) {
        if (is.null(init_alpha)) {
            a <- rep(0, N)
            if (N >= 2) a[2:N] <- 0.4
            return(a)
        }
        a <- as.numeric(init_alpha)
        if (length(a) == 1L) a <- rep(a, N)
        if (length(a) != N) stop("init_alpha must be length 1 or ncol(y) for order=1")
        a
    }

    init_alpha_order2 <- function(init_alpha, N) {
        if (is.null(init_alpha)) {
            a1 <- rep(0, N)
            a2 <- rep(0, N)
            if (N >= 2) a1[2:N] <- 0.4
            if (N >= 3) a2[3:N] <- 0.4
            return(list(a1 = a1, a2 = a2))
        }

        if (is.matrix(init_alpha) && ncol(init_alpha) >= 2) {
            a <- init_alpha
            if (nrow(a) == 1L) a <- a[rep.int(1L, N), , drop = FALSE]
            if (nrow(a) != N) stop("for order=2, init_alpha matrix must have nrow 1 or ncol(y)")
            return(list(a1 = as.numeric(a[, 1]), a2 = as.numeric(a[, 2])))
        }

        if (is.list(init_alpha) && length(init_alpha) >= 2) {
            a1 <- as.numeric(init_alpha[[1]])
            a2 <- as.numeric(init_alpha[[2]])
            if (length(a1) == 1L) a1 <- rep(a1, N)
            if (length(a2) == 1L) a2 <- rep(a2, N)
            if (length(a1) != N || length(a2) != N) stop("init_alpha list entries must be length 1 or ncol(y)")
            return(list(a1 = a1, a2 = a2))
        }

        stop("for order=2, init_alpha must be a matrix with 2 columns or a list(alpha1, alpha2)")
    }

    if (order == 0) {
        alpha_out <- NULL
    } else if (order == 1) {
        alpha1 <- init_alpha_order1(init_alpha, N)
        alpha1[1] <- 0
        alpha1 <- pmin(pmax(alpha1, 0), alpha_ub)
        alpha_out <- alpha1
    } else {
        ainit <- init_alpha_order2(init_alpha, N)
        alpha1 <- ainit$a1
        alpha2 <- ainit$a2
        alpha1[1] <- 0
        alpha2[1] <- 0
        if (N >= 2) alpha2[2] <- 0
        alpha1 <- pmin(pmax(alpha1, 0), alpha_ub)
        alpha2 <- pmin(pmax(alpha2, 0), alpha_ub)
        alpha_out <- cbind(alpha1, alpha2)
    }

    loglik_i <- rep(NA_real_, N)
    conv <- rep(NA_integer_, N)

    for (i in seq_len(N)) {
        if (order == 0) {
            if (innovation == "nbinom") {
                obj0_nb <- function(par) {
                    th <- pmax(par[1], eps_pos)
                    sz <- pmax(par[2], eps_pos)
                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = i,
                            order = 0,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = 0,
                            theta = replace(theta, i, th),
                            nb_inno_size = replace(nb_inno_size, i, sz)
                        ),
                        error = function(e) NA_real_
                    )
                    if (!is.finite(ll)) return(1e10)
                    -ll
                }

                par0 <- c(pmax(theta[i], eps_pos), pmax(nb_inno_size[i], eps_pos))
                fit0 <- optim(
                    par = par0,
                    fn = obj0_nb,
                    method = "L-BFGS-B",
                    lower = c(eps_pos, eps_pos),
                    upper = c(Inf, nb_inno_size_ub)
                )

                theta[i] <- pmax(fit0$par[1], eps_pos)
                nb_inno_size[i] <- pmin(pmax(fit0$par[2], eps_pos), nb_inno_size_ub)
                loglik_i[i] <- -fit0$value
                conv[i] <- fit0$convergence
            } else {
                obj0_th <- function(th) {
                    th <- pmax(th, eps_pos)
                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = i,
                            order = 0,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = 0,
                            theta = replace(theta, i, th),
                            nb_inno_size = NULL
                        ),
                        error = function(e) NA_real_
                    )
                    if (!is.finite(ll)) return(1e10)
                    -ll
                }

                fit0 <- optim(
                    par = pmax(theta[i], eps_pos),
                    fn = obj0_th,
                    method = "L-BFGS-B",
                    lower = eps_pos,
                    upper = Inf
                )

                theta[i] <- pmax(fit0$par, eps_pos)
                loglik_i[i] <- -fit0$value
                conv[i] <- fit0$convergence
            }
            next
        }

        if (i == 1) {
            if (innovation == "nbinom") {
                obj1_nb <- function(par) {
                    th <- pmax(par[1], eps_pos)
                    sz <- pmax(par[2], eps_pos)

                    th_try <- theta
                    sz_try <- nb_inno_size
                    th_try[1] <- th
                    sz_try[1] <- sz

                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = 1,
                            order = 0,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = 0,
                            theta = th_try,
                            nb_inno_size = sz_try
                        ),
                        error = function(e) NA_real_
                    )
                    if (!is.finite(ll)) return(1e10)
                    -ll
                }

                fit1 <- optim(
                    par = c(pmax(theta[1], eps_pos), pmax(nb_inno_size[1], eps_pos)),
                    fn = obj1_nb,
                    method = "L-BFGS-B",
                    lower = c(eps_pos, eps_pos),
                    upper = c(Inf, nb_inno_size_ub)
                )

                theta[1] <- pmax(fit1$par[1], eps_pos)
                nb_inno_size[1] <- pmin(pmax(fit1$par[2], eps_pos), nb_inno_size_ub)
                loglik_i[1] <- -fit1$value
                conv[1] <- fit1$convergence
            } else {
                obj1_th <- function(th) {
                    th <- pmax(th, eps_pos)
                    th_try <- theta
                    th_try[1] <- th

                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = 1,
                            order = 0,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = 0,
                            theta = th_try,
                            nb_inno_size = NULL
                        ),
                        error = function(e) NA_real_
                    )
                    if (!is.finite(ll)) return(1e10)
                    -ll
                }

                fit1 <- optim(
                    par = pmax(theta[1], eps_pos),
                    fn = obj1_th,
                    method = "L-BFGS-B",
                    lower = eps_pos,
                    upper = Inf
                )

                theta[1] <- pmax(fit1$par, eps_pos)
                loglik_i[1] <- -fit1$value
                conv[1] <- fit1$convergence
            }
            next
        }

        if (order == 1 || i == 2) {
            m1 <- mean_y[i - 1]

            obj_a1 <- function(par) {
                a <- pmin(pmax(par[1], 0), alpha_ub)
                mu_inno <- mean_y[i] - a * m1
                th <- theta_from_inno_mean(mu_inno, innovation, eps_pos)
                if (!is.finite(th)) return(1e10)

                alpha_try <- alpha_out
                theta_try <- theta

                if (order == 1) {
                    alpha_try[i] <- a
                } else {
                    alpha_try[i, 1] <- a
                }
                theta_try[i] <- th

                if (innovation == "nbinom") {
                    sz <- pmax(par[2], eps_pos)
                    nb_try <- nb_inno_size
                    nb_try[i] <- sz
                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = i,
                            order = order,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = alpha_try,
                            theta = theta_try,
                            nb_inno_size = nb_try
                        ),
                        error = function(e) NA_real_
                    )
                } else {
                    ll <- tryCatch(
                        logL_inad_i(
                            y = y,
                            i = i,
                            order = order,
                            thinning = thinning,
                            innovation = innovation,
                            alpha = alpha_try,
                            theta = theta_try,
                            nb_inno_size = NULL
                        ),
                        error = function(e) NA_real_
                    )
                }

                if (!is.finite(ll)) return(1e10)
                -ll
            }

            if (order == 1) {
                alpha_start <- alpha_out[i]
                if (innovation == "nbinom") {
                    fit_i <- optim(
                        par = c(alpha_start, nb_inno_size[i]),
                        fn = obj_a1,
                        method = "L-BFGS-B",
                        lower = c(0, eps_pos),
                        upper = c(alpha_ub, nb_inno_size_ub)
                    )
                    alpha_out[i] <- pmin(pmax(fit_i$par[1], 0), alpha_ub)
                    nb_inno_size[i] <- pmin(pmax(fit_i$par[2], eps_pos), nb_inno_size_ub)
                } else {
                    fit_i <- optim(
                        par = alpha_start,
                        fn = function(a) obj_a1(c(a)),
                        method = "L-BFGS-B",
                        lower = 0,
                        upper = alpha_ub
                    )
                    alpha_out[i] <- pmin(pmax(fit_i$par, 0), alpha_ub)
                }
            } else {
                alpha_start <- alpha_out[i, 1]
                if (innovation == "nbinom") {
                    fit_i <- optim(
                        par = c(alpha_start, nb_inno_size[i]),
                        fn = obj_a1,
                        method = "L-BFGS-B",
                        lower = c(0, eps_pos),
                        upper = c(alpha_ub, nb_inno_size_ub)
                    )
                    alpha_out[i, 1] <- pmin(pmax(fit_i$par[1], 0), alpha_ub)
                    nb_inno_size[i] <- pmin(pmax(fit_i$par[2], eps_pos), nb_inno_size_ub)
                } else {
                    fit_i <- optim(
                        par = alpha_start,
                        fn = function(a) obj_a1(c(a)),
                        method = "L-BFGS-B",
                        lower = 0,
                        upper = alpha_ub
                    )
                    alpha_out[i, 1] <- pmin(pmax(fit_i$par, 0), alpha_ub)
                }
            }

            if (order == 1) {
                mu_inno <- mean_y[i] - alpha_out[i] * m1
                th <- theta_from_inno_mean(mu_inno, innovation, eps_pos)
                if (!is.finite(th)) th <- eps_pos
                theta[i] <- th
            } else {
                theta[i] <- pmax(mean_y[i] - alpha_out[i, 1] * m1, eps_pos)
                alpha_out[i, 2] <- 0
            }

            loglik_i[i] <- -fit_i$value
            conv[i] <- fit_i$convergence
            next
        }

        m1 <- mean_y[i - 1]
        m2 <- mean_y[i - 2]

        obj_a12 <- function(par) {
            a1 <- pmin(pmax(par[1], 0), alpha_ub)
            a2 <- pmin(pmax(par[2], 0), alpha_ub)
            mu_inno <- mean_y[i] - a1 * m1 - a2 * m2
            th <- theta_from_inno_mean(mu_inno, innovation, eps_pos)
            if (!is.finite(th)) return(1e10)

            alpha_try <- alpha_out
            alpha_try[i, 1] <- a1
            alpha_try[i, 2] <- a2

            theta_try <- theta
            theta_try[i] <- th

            if (innovation == "nbinom") {
                sz <- pmax(par[3], eps_pos)
                nb_try <- nb_inno_size
                nb_try[i] <- sz
                ll <- tryCatch(
                    logL_inad_i(
                        y = y,
                        i = i,
                        order = order,
                        thinning = thinning,
                        innovation = innovation,
                        alpha = alpha_try,
                        theta = theta_try,
                        nb_inno_size = nb_try
                    ),
                    error = function(e) NA_real_
                )
            } else {
                ll <- tryCatch(
                    logL_inad_i(
                        y = y,
                        i = i,
                        order = order,
                        thinning = thinning,
                        innovation = innovation,
                        alpha = alpha_try,
                        theta = theta_try,
                        nb_inno_size = NULL
                    ),
                    error = function(e) NA_real_
                )
            }

            if (!is.finite(ll)) return(1e10)
            -ll
        }

        if (innovation == "nbinom") {
            fit_i <- optim(
                par = c(alpha_out[i, 1], alpha_out[i, 2], nb_inno_size[i]),
                fn = obj_a12,
                method = "L-BFGS-B",
                lower = c(0, 0, eps_pos),
                upper = c(alpha_ub, alpha_ub, nb_inno_size_ub)
            )
            alpha_out[i, 1] <- pmin(pmax(fit_i$par[1], 0), alpha_ub)
            alpha_out[i, 2] <- pmin(pmax(fit_i$par[2], 0), alpha_ub)
            nb_inno_size[i] <- pmin(pmax(fit_i$par[3], eps_pos), nb_inno_size_ub)
        } else {
            fit_i <- optim(
                par = c(alpha_out[i, 1], alpha_out[i, 2]),
                fn = function(p) obj_a12(p),
                method = "L-BFGS-B",
                lower = c(0, 0),
                upper = c(alpha_ub, alpha_ub)
            )
            alpha_out[i, 1] <- pmin(pmax(fit_i$par[1], 0), alpha_ub)
            alpha_out[i, 2] <- pmin(pmax(fit_i$par[2], 0), alpha_ub)
        }

        mu_inno <- mean_y[i] - alpha_out[i, 1] * m1 - alpha_out[i, 2] * m2
        th <- theta_from_inno_mean(mu_inno, innovation, eps_pos)
        if (!is.finite(th)) th <- eps_pos
        theta[i] <- th

        loglik_i[i] <- -fit_i$value
        conv[i] <- fit_i$convergence
    }

    if (!is.null(alpha_out)) alpha_out <- unname(alpha_out)
    theta <- unname(theta)
    if (!is.null(nb_inno_size)) nb_inno_size <- unname(nb_inno_size)

    list(
        alpha = alpha_out,
        theta = theta,
        nb_inno_size = nb_inno_size,
        tau = NULL,
        log_l = sum(loglik_i),
        loglik_i = loglik_i,
        convergence = list(per_time = conv),
        settings = list(
            order = order,
            thinning = thinning,
            innovation = innovation,
            blocks = NULL,
            n_subjects = nrow(y),
            n_time = ncol(y)
        )
    )
}

fit_inad_fe <- function(
        y,
        order,
        thinning,
        innovation,
        blocks,
        max_iter,
        tol,
        verbose,
        init_alpha,
        init_theta,
        init_tau,
        init_nb_inno_size,
        nb_inno_size_ub
) {
    if (!requireNamespace("nloptr", quietly = TRUE)) {
        stop("nloptr is required for fixed effect fitting")
    }

    n <- nrow(y)
    N <- ncol(y)

    blocks <- as.integer(blocks)
    if (length(blocks) != n) stop("length(blocks) must equal nrow(y)")
    if (any(blocks < 1)) stop("blocks must be positive integers starting at 1")
    B <- max(blocks)

    eps_alpha <- 1e-10
    eps_pos <- 1e-10
    alpha_ub <- 1 - eps_alpha

    default_theta0 <- function(innovation) {
        if (innovation == "bell") return(1)
        10
    }

    theta <- if (is.null(init_theta)) rep(default_theta0(innovation), N) else {
        th <- as.numeric(init_theta)
        if (length(th) == 1L) th <- rep(th, N)
        if (length(th) != N) stop("init_theta must be length 1 or ncol(y)")
        pmax(th, eps_pos)
    }

    nb_inno_size <- NULL
    if (innovation == "nbinom") {
        nb_inno_size <- as.numeric(init_nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("init_nb_inno_size must be length 1 or ncol(y)")
        nb_inno_size <- pmin(pmax(nb_inno_size, eps_pos), nb_inno_size_ub)
    }

    init_alpha_order1 <- function(init_alpha, N) {
        if (is.null(init_alpha)) {
            a <- rep(0, N)
            if (N >= 2) a[2:N] <- 0.4
            return(a)
        }
        a <- as.numeric(init_alpha)
        if (length(a) == 1L) a <- rep(a, N)
        if (length(a) != N) stop("init_alpha must be length 1 or ncol(y) for order=1")
        a
    }

    init_alpha_order2 <- function(init_alpha, N) {
        if (is.null(init_alpha)) {
            a1 <- rep(0, N)
            a2 <- rep(0, N)
            if (N >= 2) a1[2:N] <- 0.4
            if (N >= 3) a2[3:N] <- 0.4
            return(list(a1 = a1, a2 = a2))
        }

        if (is.matrix(init_alpha) && ncol(init_alpha) >= 2) {
            a <- init_alpha
            if (nrow(a) == 1L) a <- a[rep.int(1L, N), , drop = FALSE]
            if (nrow(a) != N) stop("for order=2, init_alpha matrix must have nrow 1 or ncol(y)")
            return(list(a1 = as.numeric(a[, 1]), a2 = as.numeric(a[, 2])))
        }

        if (is.list(init_alpha) && length(init_alpha) >= 2) {
            a1 <- as.numeric(init_alpha[[1]])
            a2 <- as.numeric(init_alpha[[2]])
            if (length(a1) == 1L) a1 <- rep(a1, N)
            if (length(a2) == 1L) a2 <- rep(a2, N)
            if (length(a1) != N || length(a2) != N) stop("init_alpha list entries must be length 1 or ncol(y)")
            return(list(a1 = a1, a2 = a2))
        }

        stop("for order=2, init_alpha must be a matrix with 2 columns or a list(alpha1, alpha2)")
    }

    alpha_full <- NULL
    if (order == 1) {
        alpha1 <- init_alpha_order1(init_alpha, N)
        alpha1[1] <- 0
        alpha1 <- pmin(pmax(alpha1, 0), alpha_ub)
        alpha_full <- alpha1
    }
    if (order == 2) {
        ainit <- init_alpha_order2(init_alpha, N)
        alpha1 <- ainit$a1
        alpha2 <- ainit$a2
        alpha1[1] <- 0
        alpha2[1] <- 0
        if (N >= 2) alpha2[2] <- 0
        alpha1 <- pmin(pmax(alpha1, 0), alpha_ub)
        alpha2 <- pmin(pmax(alpha2, 0), alpha_ub)
        alpha_full <- cbind(alpha1, alpha2)
    }

    tau <- as.numeric(init_tau)
    if (B <= 1L) {
        tau <- 0
    } else if (length(tau) == 1L) {
        tau <- c(0, rep(tau, B - 1L))
    } else {
        if (length(tau) < B) stop("init_tau must be length 1 or at least max(blocks)")
        tau <- tau[seq_len(B)]
        tau[1] <- 0
    }

    lambda_ok <- function(theta_vec, tau_vec) {
        mu <- if (innovation == "bell") {
            outer(rep(1, n), .bell_mean_from_theta(theta_vec)) + matrix(tau_vec[blocks], n, N)
        } else {
            outer(rep(1, n), theta_vec) + matrix(tau_vec[blocks], n, N)
        }
        all(is.finite(mu)) && all(mu > eps_pos)
    }

    log_old <- logL_inad(
        y = y,
        order = order,
        thinning = thinning,
        innovation = innovation,
        alpha = alpha_full,
        theta = theta,
        nb_inno_size = nb_inno_size,
        blocks = blocks,
        tau = tau
    )

    converged <- FALSE
    for (iter in seq_len(max_iter)) {
        if (verbose) message("iter ", iter, " logL=", log_old)

        if (B > 1L) {
            obj_tau <- function(tau_free) {
                tau_try <- c(0, tau_free)
                if (!lambda_ok(theta, tau_try)) return(1e8)
                -logL_inad(
                    y = y,
                    order = order,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = alpha_full,
                    theta = theta,
                    nb_inno_size = nb_inno_size,
                    blocks = blocks,
                    tau = tau_try
                )
            }

            opt_tau <- nloptr::nloptr(
                x0 = tau[2:B],
                eval_f = obj_tau,
                opts = list(
                    algorithm = "NLOPT_LN_BOBYQA",
                    xtol_rel = 1e-6,
                    maxeval = 1000
                )
            )
            tau <- c(0, opt_tau$solution)
        } else {
            tau <- 0
        }

        if (order == 1) {
            if (N >= 2) {
                obj_a <- function(a_free) {
                    a_try <- alpha_full
                    a_try[1] <- 0
                    a_try[2:N] <- pmin(pmax(a_free, 0), alpha_ub)
                    -logL_inad(
                        y = y,
                        order = order,
                        thinning = thinning,
                        innovation = innovation,
                        alpha = a_try,
                        theta = theta,
                        nb_inno_size = nb_inno_size,
                        blocks = blocks,
                        tau = tau
                    )
                }

                opt_a <- nloptr::nloptr(
                    x0 = alpha_full[2:N],
                    eval_f = obj_a,
                    lb = rep(0, N - 1),
                    ub = rep(alpha_ub, N - 1),
                    opts = list(
                        algorithm = "NLOPT_LN_BOBYQA",
                        xtol_rel = 1e-6,
                        maxeval = 1500
                    )
                )

                alpha_full[2:N] <- pmin(pmax(opt_a$solution, 0), alpha_ub)
                alpha_full[1] <- 0
            }
        }

        if (order == 2) {
            if (N >= 3) {
                obj_a <- function(a_free) {
                    a1 <- alpha_full[, 1]
                    a2 <- alpha_full[, 2]

                    a1[1] <- 0
                    a2[1] <- 0
                    a2[2] <- 0

                    a1[2:N] <- pmin(pmax(a_free[seq_len(N - 1)], 0), alpha_ub)
                    a2[3:N] <- pmin(pmax(a_free[(N):(N + (N - 3))], 0), alpha_ub)

                    a_try <- cbind(a1, a2)

                    -logL_inad(
                        y = y,
                        order = order,
                        thinning = thinning,
                        innovation = innovation,
                        alpha = a_try,
                        theta = theta,
                        nb_inno_size = nb_inno_size,
                        blocks = blocks,
                        tau = tau
                    )
                }

                x0 <- c(alpha_full[2:N, 1], alpha_full[3:N, 2])
                lb <- rep(0, length(x0))
                ub <- rep(alpha_ub, length(x0))

                opt_a <- nloptr::nloptr(
                    x0 = x0,
                    eval_f = obj_a,
                    lb = lb,
                    ub = ub,
                    opts = list(
                        algorithm = "NLOPT_LN_BOBYQA",
                        xtol_rel = 1e-6,
                        maxeval = 2000
                    )
                )

                sol <- opt_a$solution
                alpha_full[1, 1] <- 0
                alpha_full[1, 2] <- 0
                alpha_full[2, 2] <- 0
                alpha_full[2:N, 1] <- pmin(pmax(sol[seq_len(N - 1)], 0), alpha_ub)
                alpha_full[3:N, 2] <- pmin(pmax(sol[(N):(N + (N - 3))], 0), alpha_ub)
            }
        }

        obj_t <- function(t_try) {
            t_try <- pmax(t_try, eps_pos)
            if (!lambda_ok(t_try, tau)) return(1e8)
            -logL_inad(
                y = y,
                order = order,
                thinning = thinning,
                innovation = innovation,
                alpha = alpha_full,
                theta = t_try,
                nb_inno_size = nb_inno_size,
                blocks = blocks,
                tau = tau
            )
        }

        opt_t <- nloptr::nloptr(
            x0 = theta,
            eval_f = obj_t,
            lb = rep(eps_pos, N),
            opts = list(
                algorithm = "NLOPT_LN_BOBYQA",
                xtol_rel = 1e-6,
                maxeval = 2000
            )
        )
        theta <- pmax(opt_t$solution, eps_pos)

        if (innovation == "nbinom") {
            obj_sz <- function(sz_try) {
                sz_try <- pmax(sz_try, eps_pos)
                -logL_inad(
                    y = y,
                    order = order,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = alpha_full,
                    theta = theta,
                    nb_inno_size = sz_try,
                    blocks = blocks,
                    tau = tau
                )
            }

            opt_sz <- nloptr::nloptr(
                x0 = nb_inno_size,
                eval_f = obj_sz,
                lb = rep(eps_pos, N),
                ub = rep(nb_inno_size_ub, N),
                opts = list(
                    algorithm = "NLOPT_LN_BOBYQA",
                    xtol_rel = 1e-6,
                    maxeval = 2000
                )
            )
            nb_inno_size <- pmin(pmax(opt_sz$solution, eps_pos), nb_inno_size_ub)
        }

        log_new <- logL_inad(
            y = y,
            order = order,
            thinning = thinning,
            innovation = innovation,
            alpha = alpha_full,
            theta = theta,
            nb_inno_size = nb_inno_size,
            blocks = blocks,
            tau = tau
        )

        if (!is.finite(log_new)) break
        if (abs(log_new - log_old) < tol) {
            converged <- TRUE
            log_old <- log_new
            break
        }
        log_old <- log_new
    }

    list(
        alpha = alpha_full,
        theta = theta,
        tau = tau,
        nb_inno_size = nb_inno_size,
        log_l = log_old,
        convergence = if (converged) 0L else 1L,
        convergence_info = list(converged = converged, max_iter = max_iter),
        settings = list(
            order = order,
            thinning = thinning,
            innovation = innovation,
            blocks = blocks,
            n_subjects = n,
            n_time = N
        )
    )
}

solve_theta_bell <- function(mu, lower = 0, upper = 10) {
    if (!is.finite(mu) || mu <= 0) return(NA_real_)

    f <- function(theta) theta * exp(theta) - mu

    up <- upper
    fu <- f(up)
    while (is.finite(fu) && fu < 0 && up < 100) {
        up <- up * 2
        fu <- f(up)
    }
    if (!is.finite(fu) || fu < 0) return(NA_real_)

    uniroot(f, c(lower, up))$root
}
