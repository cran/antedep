# File: R/test_stationarity_gau.R
# Likelihood ratio tests for Gaussian AD stationarity constraints

#' Likelihood ratio test for stationarity (Gaussian AD data)
#'
#' Tests whether time-varying Gaussian AD covariance parameters can be
#' constrained to be constant over time.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param order Antedependence order (0, 1, or 2).
#' @param blocks Optional vector of block memberships (length n_subjects).
#' @param constrain Constraint to test:
#'   for order 0: "sigma" (or "all");
#'   for order 1: "phi", "sigma", or "both"/"all";
#'   for order 2: "phi1", "phi2", "phi", "sigma", or "all"/"both".
#' @param fit_unconstrained Optional pre-computed unconstrained fit from
#'   \code{\link{fit_gau}}.
#' @param verbose Logical; if TRUE, prints fitting progress.
#' @param max_iter Maximum number of optimization iterations for constrained fit.
#' @param rel_tol Relative tolerance for constrained optimization.
#' @param ... Additional arguments passed to \code{\link{fit_gau}} when
#'   \code{fit_unconstrained} is not provided.
#'
#' @return A list with class \code{"test_stationarity_gau"} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{fit_unconstrained}{Unconstrained Gaussian AD fit}
#'   \item{fit_constrained}{Constrained Gaussian AD fit}
#'   \item{constraint}{Human-readable null constraint description}
#'   \item{lrt_stat}{Likelihood-ratio statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{Chi-square p-value}
#'   \item{bic_unconstrained}{BIC of unconstrained model}
#'   \item{bic_constrained}{BIC of constrained model}
#'   \item{bic_selected}{Model selected by BIC}
#'   \item{table}{Two-row model summary table}
#' }
#'
#' @details
#' The mean structure is kept unrestricted in both models (time-specific means
#' plus optional block shifts), and the test constrains covariance parameters:
#' innovation standard deviations and/or antedependence coefficients.
#'
#' The likelihood-ratio statistic is:
#' \deqn{\lambda = 2(\ell_{alt} - \ell_{null})}
#' where \eqn{\ell_{null}} and \eqn{\ell_{alt}} are maximized log-likelihoods
#' under the constrained and unconstrained models, respectively.
#'
#' Degrees of freedom are computed from the number of constraints implied by
#' \code{constrain}.
#'
#' @references
#' Zimmerman, D.L. and Nunez-Anton, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Chapter 6.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 80, n_time = 6, order = 1, phi = 0.4, sigma = 1)
#'
#' # Test jointly constant phi and sigma (order 1)
#' out <- test_stationarity_gau(y, order = 1, constrain = "both")
#' out$p_value
#'
#' @seealso \code{\link{run_stationarity_tests_gau}}, \code{\link{test_order_gau}},
#'   \code{\link{fit_gau}}
#'
#' @export
test_stationarity_gau <- function(y, order = 1L, blocks = NULL, constrain = "both",
                                 fit_unconstrained = NULL, verbose = FALSE,
                                 max_iter = 2000L, rel_tol = 1e-8, ...) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "test_stationarity_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y))) stop("y must contain finite values", call. = FALSE)

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    order <- as.integer(order)
    if (!order %in% c(0L, 1L, 2L)) stop("order must be 0, 1, or 2", call. = FALSE)
    if (order >= n_time) stop("order must be less than n_time", call. = FALSE)

    block_info <- .normalize_blocks(blocks, n_subjects)
    blocks_id <- block_info$blocks_id
    n_blocks <- block_info$n_blocks
    block_levels <- block_info$block_levels
    blocks_fit <- if (is.null(blocks)) NULL else blocks_id

    cinfo <- .parse_stationarity_constraint_gau(constrain, order, n_time)

    fit_args <- list(...)
    if (is.null(fit_unconstrained)) {
        if (verbose) cat("Fitting unconstrained Gaussian AD model...\n")
        fit_unconstrained <- do.call(
            fit_gau,
            c(
                list(
                    y = y,
                    order = order,
                    blocks = blocks_fit,
                    na_action = "fail",
                    estimate_mu = TRUE
                ),
                fit_args
            )
        )
    } else {
        if (!is.list(fit_unconstrained) || is.null(fit_unconstrained$log_l)) {
            stop("fit_unconstrained must be a fit_gau-like object with log_l.", call. = FALSE)
        }
        fit_order <- fit_unconstrained$settings$order
        if (!is.null(fit_order) && as.integer(fit_order) != order) {
            stop("fit_unconstrained$settings$order does not match requested order.", call. = FALSE)
        }
    }

    if (verbose) {
        cat("Fitting constrained Gaussian AD model under H0:", cinfo$description, "\n")
    }
    fit_constrained <- .fit_constrained_stationarity_gau(
        y = y,
        order = order,
        blocks = blocks_fit,
        cinfo = cinfo,
        fit_init = fit_unconstrained,
        block_levels = if (is.null(blocks)) NULL else block_levels,
        max_iter = as.integer(max_iter),
        rel_tol = as.numeric(rel_tol),
        verbose = isTRUE(verbose)
    )

    logL_uncon <- fit_unconstrained$log_l
    logL_con <- fit_constrained$log_l

    lrt_stat <- 2 * (logL_uncon - logL_con)
    if (!is.finite(lrt_stat)) {
        lrt_stat <- NA_real_
    } else if (lrt_stat < 0) {
        if (lrt_stat > -1e-6) {
            lrt_stat <- 0
        } else {
            warning(
                "LRT statistic was negative. Returning zero (numerical tolerance).",
                call. = FALSE
            )
            lrt_stat <- 0
        }
    }

    df <- as.integer(cinfo$df)
    p_value <- if (is.finite(lrt_stat)) {
        stats::pchisq(lrt_stat, df = df, lower.tail = FALSE)
    } else {
        NA_real_
    }

    n_mean_params <- n_time + max(n_blocks - 1L, 0L)
    n_cov_uncon <- .count_params_gau(order, n_time, n_blocks = 1L, estimate_mu = FALSE)
    n_params_uncon <- n_mean_params + n_cov_uncon
    n_params_con <- n_params_uncon - df

    bic_uncon <- -2 * logL_uncon + n_params_uncon * log(n_subjects)
    bic_con <- -2 * logL_con + n_params_con * log(n_subjects)
    bic_selected <- if (bic_uncon <= bic_con) "unconstrained" else "constrained"

    table_df <- data.frame(
        model = c(
            "Unconstrained",
            paste0("Constrained: ", cinfo$description)
        ),
        logLik = c(logL_uncon, logL_con),
        n_params = c(n_params_uncon, n_params_con),
        BIC = c(bic_uncon, bic_con),
        stringsAsFactors = FALSE
    )

    out <- list(
        method = "lrt",
        fit_unconstrained = fit_unconstrained,
        fit_constrained = fit_constrained,
        constraint = cinfo$description,
        lrt_stat = lrt_stat,
        statistic = lrt_stat,
        df = df,
        p_value = p_value,
        bic_unconstrained = bic_uncon,
        bic_constrained = bic_con,
        bic_selected = bic_selected,
        table = table_df,
        settings = list(
            order = order,
            n_subjects = n_subjects,
            n_time = n_time,
            n_blocks = n_blocks,
            constrain = constrain,
            blocks = blocks_fit,
            block_levels = if (is.null(blocks)) NULL else block_levels,
            max_iter = as.integer(max_iter),
            rel_tol = as.numeric(rel_tol)
        )
    )

    class(out) <- "test_stationarity_gau"
    out
}


#' Run all stationarity-related tests for Gaussian AD
#'
#' Runs a standard set of stationarity constraints for Gaussian AD models.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param order Antedependence order (0, 1, or 2).
#' @param blocks Optional vector of block memberships (length n_subjects).
#' @param verbose Logical; if TRUE, prints progress.
#' @param max_iter Maximum number of optimization iterations for constrained fits.
#' @param rel_tol Relative tolerance for constrained optimization.
#' @param ... Additional arguments passed to \code{\link{fit_gau}} for the
#'   unconstrained fit.
#'
#' @return A list with class \code{"stationarity_tests_gau"} containing:
#' \describe{
#'   \item{fit_unconstrained}{Unconstrained Gaussian AD fit}
#'   \item{tests}{Named list of \code{\link{test_stationarity_gau}} results}
#'   \item{summary}{Summary table of all constraints}
#' }
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 80, n_time = 6, order = 1, phi = 0.4, sigma = 1)
#' out <- run_stationarity_tests_gau(y, order = 1, verbose = FALSE)
#' out$summary
#'
#' @seealso \code{\link{test_stationarity_gau}}
#'
#' @export
run_stationarity_tests_gau <- function(y, order = 1L, blocks = NULL, verbose = FALSE,
                                       max_iter = 2000L, rel_tol = 1e-8, ...) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "run_stationarity_tests_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    order <- as.integer(order)
    if (!order %in% c(0L, 1L, 2L)) stop("order must be 0, 1, or 2", call. = FALSE)
    if (order >= n_time) stop("order must be less than n_time", call. = FALSE)

    block_info <- .normalize_blocks(blocks, n_subjects)
    blocks_fit <- if (is.null(blocks)) NULL else block_info$blocks_id

    fit_args <- list(...)
    if (verbose) cat("Fitting unconstrained Gaussian AD model...\n")
    fit_uncon <- do.call(
        fit_gau,
        c(
            list(
                y = y,
                order = order,
                blocks = blocks_fit,
                na_action = "fail",
                estimate_mu = TRUE
            ),
            fit_args
        )
    )

    tests_to_run <- if (order == 0L) {
        c("sigma")
    } else if (order == 1L) {
        c("phi", "sigma", "both")
    } else {
        c("phi1", "phi2", "phi", "sigma", "all")
    }

    results <- vector("list", length(tests_to_run))
    names(results) <- tests_to_run

    for (name in tests_to_run) {
        if (verbose) cat("Testing constraint:", name, "\n")
        results[[name]] <- test_stationarity_gau(
            y = y,
            order = order,
            blocks = blocks_fit,
            constrain = name,
            fit_unconstrained = fit_uncon,
            verbose = FALSE,
            max_iter = max_iter,
            rel_tol = rel_tol
        )
    }

    summary_df <- data.frame(
        constraint = tests_to_run,
        method = vapply(results, function(x) x$method, character(1)),
        df = vapply(results, function(x) x$df, numeric(1)),
        LRT = vapply(results, function(x) x$lrt_stat, numeric(1)),
        p_value = vapply(results, function(x) x$p_value, numeric(1)),
        BIC_selected = vapply(results, function(x) x$bic_selected, character(1)),
        stringsAsFactors = FALSE
    )

    out <- list(
        method = "lrt",
        fit_unconstrained = fit_uncon,
        tests = results,
        summary = summary_df,
        settings = list(
            order = order,
            n_subjects = n_subjects,
            n_time = n_time,
            n_blocks = block_info$n_blocks,
            blocks = blocks_fit
        )
    )
    class(out) <- "stationarity_tests_gau"
    out
}


#' @export
print.test_stationarity_gau <- function(x, digits = 4, ...) {
    cat("\nLikelihood Ratio Test for Gaussian AD Stationarity\n")
    cat("==================================================\n\n")
    cat("Constraint: H0:", x$constraint, "\n\n")
    print(x$table, row.names = FALSE, digits = digits)
    cat(
        "\nLRT:", round(x$lrt_stat, digits),
        " df:", x$df,
        " p-value:", format.pval(x$p_value, digits = digits)
    )
    cat("\nBIC selects:", x$bic_selected, "\n")
    invisible(x)
}


#' @export
print.stationarity_tests_gau <- function(x, digits = 4, ...) {
    cat("\nStationarity Tests for Gaussian AD Model\n")
    cat("========================================\n\n")
    print(x$summary, row.names = FALSE, digits = digits)
    invisible(x)
}


#' @keywords internal
.parse_stationarity_constraint_gau <- function(constrain, order, n_time) {
    constrain <- tolower(as.character(constrain))
    if (length(constrain) == 0L || all(is.na(constrain))) {
        stop("constrain must be provided.", call. = FALSE)
    }

    if (order == 0L) {
        sigma_const <- any(constrain %in% c("sigma", "all", "both"))
        if (!sigma_const) {
            stop("For order 0, constrain must be 'sigma' (or 'all').", call. = FALSE)
        }
        df <- n_time - 1L
        return(list(
            order = 0L,
            sigma_const = TRUE,
            df = df,
            description = "sigma constant over time"
        ))
    }

    if (order == 1L) {
        if (length(constrain) == 1L) {
            if (constrain %in% c("both", "all")) {
                phi_const <- TRUE
                sigma_const <- TRUE
            } else if (constrain == "phi") {
                phi_const <- TRUE
                sigma_const <- FALSE
            } else if (constrain == "sigma") {
                phi_const <- FALSE
                sigma_const <- TRUE
            } else {
                stop("Invalid constraint for order 1. Use 'phi', 'sigma', or 'both'.", call. = FALSE)
            }
        } else {
            phi_const <- "phi" %in% constrain
            sigma_const <- "sigma" %in% constrain
            if ("both" %in% constrain || "all" %in% constrain) {
                phi_const <- TRUE
                sigma_const <- TRUE
            }
        }

        df <- 0L
        if (phi_const) df <- df + (n_time - 2L)
        if (sigma_const) df <- df + (n_time - 1L)
        if (df <= 0L) {
            stop("Selected constraint does not impose restrictions (df = 0).", call. = FALSE)
        }

        desc <- paste(
            c(
                if (phi_const) "phi constant over time",
                if (sigma_const) "sigma constant over time"
            ),
            collapse = ", "
        )
        return(list(
            order = 1L,
            phi_const = phi_const,
            sigma_const = sigma_const,
            df = df,
            description = desc
        ))
    }

    if (length(constrain) == 1L) {
        if (constrain %in% c("all", "both")) {
            phi1_const <- TRUE
            phi2_const <- TRUE
            sigma_const <- TRUE
        } else if (constrain == "phi") {
            phi1_const <- TRUE
            phi2_const <- TRUE
            sigma_const <- FALSE
        } else if (constrain == "phi1") {
            phi1_const <- TRUE
            phi2_const <- FALSE
            sigma_const <- FALSE
        } else if (constrain == "phi2") {
            phi1_const <- FALSE
            phi2_const <- TRUE
            sigma_const <- FALSE
        } else if (constrain == "sigma") {
            phi1_const <- FALSE
            phi2_const <- FALSE
            sigma_const <- TRUE
        } else {
            stop(
                "Invalid constraint for order 2. Use 'phi1', 'phi2', 'phi', 'sigma', or 'all'.",
                call. = FALSE
            )
        }
    } else {
        phi1_const <- "phi1" %in% constrain
        phi2_const <- "phi2" %in% constrain
        sigma_const <- "sigma" %in% constrain
        if ("phi" %in% constrain) {
            phi1_const <- TRUE
            phi2_const <- TRUE
        }
        if ("all" %in% constrain || "both" %in% constrain) {
            phi1_const <- TRUE
            phi2_const <- TRUE
            sigma_const <- TRUE
        }
    }

    df <- 0L
    if (phi1_const) df <- df + (n_time - 2L)
    if (phi2_const) df <- df + max(n_time - 3L, 0L)
    if (sigma_const) df <- df + (n_time - 1L)
    if (df <= 0L) {
        stop("Selected constraint does not impose restrictions (df = 0).", call. = FALSE)
    }

    desc <- paste(
        c(
            if (phi1_const) "phi1 constant over time",
            if (phi2_const) "phi2 constant over time",
            if (sigma_const) "sigma constant over time"
        ),
        collapse = ", "
    )

    list(
        order = 2L,
        phi1_const = phi1_const,
        phi2_const = phi2_const,
        sigma_const = sigma_const,
        df = df,
        description = desc
    )
}


#' @keywords internal
.fit_constrained_stationarity_gau <- function(y, order, blocks, cinfo, fit_init,
                                              block_levels = NULL, max_iter = 2000L,
                                              rel_tol = 1e-8, verbose = FALSE) {
    n_subjects <- nrow(y)
    n_time <- ncol(y)
    block_info <- .normalize_blocks(blocks, n_subjects)
    blocks_id <- block_info$blocks_id
    n_blocks <- block_info$n_blocks
    blocks_fit <- if (is.null(blocks)) NULL else blocks_id

    mu_init <- if (!is.null(fit_init$mu) && length(fit_init$mu) == n_time) {
        as.numeric(fit_init$mu)
    } else {
        colMeans(y)
    }

    tau_init <- if (n_blocks > 1L) {
        if (!is.null(fit_init$tau) && length(fit_init$tau) == n_blocks) {
            as.numeric(fit_init$tau)
        } else {
            rep(0, n_blocks)
        }
    } else {
        0
    }
    if (n_blocks > 1L) tau_init[1] <- 0

    sigma_full_init <- if (!is.null(fit_init$sigma) && length(fit_init$sigma) == n_time) {
        pmax(as.numeric(fit_init$sigma), 1e-8)
    } else {
        rep(max(stats::sd(as.numeric(y)), 1e-2), n_time)
    }

    phi1_init <- NULL
    phi2_init <- NULL
    if (order == 1L) {
        phi_full <- rep(0, n_time)
        if (!is.null(fit_init$phi)) {
            phi_raw <- as.numeric(fit_init$phi)
            if (length(phi_raw) == n_time) {
                phi_full <- phi_raw
            } else if (length(phi_raw) == n_time - 1L) {
                phi_full <- c(0, phi_raw)
            }
        }
        phi1_init <- if (n_time >= 2L) phi_full[2:n_time] else numeric(0)
    } else if (order == 2L) {
        phi_mat <- matrix(0, nrow = 2L, ncol = n_time)
        if (is.matrix(fit_init$phi) &&
            nrow(fit_init$phi) == 2L &&
            ncol(fit_init$phi) == n_time) {
            phi_mat <- fit_init$phi
        }
        if (n_time >= 2L) phi1_init <- phi_mat[1, 2:n_time]
        if (n_time >= 3L) phi2_init <- phi_mat[2, 3:n_time]
    }

    par0 <- numeric(0)
    par0 <- c(par0, mu_init)

    if (n_blocks > 1L) {
        par0 <- c(par0, tau_init[-1])
    }

    if (isTRUE(cinfo$sigma_const)) {
        par0 <- c(par0, log(exp(mean(log(sigma_full_init)))))
    } else {
        par0 <- c(par0, log(sigma_full_init))
    }

    if (order == 1L) {
        if (isTRUE(cinfo$phi_const)) {
            par0 <- c(par0, mean(phi1_init))
        } else {
            par0 <- c(par0, phi1_init)
        }
    } else if (order == 2L) {
        if (isTRUE(cinfo$phi1_const)) {
            par0 <- c(par0, mean(phi1_init))
        } else {
            par0 <- c(par0, phi1_init)
        }

        if (isTRUE(cinfo$phi2_const)) {
            par0 <- c(par0, mean(phi2_init))
        } else if (n_time >= 3L) {
            par0 <- c(par0, phi2_init)
        }
    }

    unpack <- function(par) {
        idx <- 1L

        mu <- par[idx:(idx + n_time - 1L)]
        idx <- idx + n_time

        if (n_blocks > 1L) {
            tau <- c(0, par[idx:(idx + n_blocks - 2L)])
            idx <- idx + n_blocks - 1L
        } else {
            tau <- 0
        }

        if (isTRUE(cinfo$sigma_const)) {
            sigma_scalar <- exp(par[idx])
            sigma <- rep(sigma_scalar, n_time)
            idx <- idx + 1L
        } else {
            sigma <- exp(par[idx:(idx + n_time - 1L)])
            idx <- idx + n_time
        }

        phi <- NULL

        if (order == 1L) {
            if (n_time >= 2L) {
                if (isTRUE(cinfo$phi_const)) {
                    phi_tail <- rep(par[idx], n_time - 1L)
                    idx <- idx + 1L
                } else {
                    phi_tail <- par[idx:(idx + n_time - 2L)]
                    idx <- idx + n_time - 1L
                }
                phi <- c(0, phi_tail)
            } else {
                phi <- numeric(0)
            }
        } else if (order == 2L) {
            phi_mat <- matrix(0, nrow = 2L, ncol = n_time)

            if (n_time >= 2L) {
                if (isTRUE(cinfo$phi1_const)) {
                    phi1 <- rep(par[idx], n_time - 1L)
                    idx <- idx + 1L
                } else {
                    phi1 <- par[idx:(idx + n_time - 2L)]
                    idx <- idx + n_time - 1L
                }
                phi_mat[1, 2:n_time] <- phi1
            }

            if (n_time >= 3L) {
                if (isTRUE(cinfo$phi2_const)) {
                    phi2 <- rep(par[idx], n_time - 2L)
                    idx <- idx + 1L
                } else {
                    phi2 <- par[idx:(idx + n_time - 3L)]
                    idx <- idx + n_time - 2L
                }
                phi_mat[2, 3:n_time] <- phi2
            }

            phi <- phi_mat
        }

        list(mu = mu, tau = tau, sigma = sigma, phi = phi)
    }

    objective <- function(par) {
        params <- unpack(par)
        log_l <- logL_gau(
            y = y,
            order = order,
            mu = params$mu,
            phi = params$phi,
            sigma = params$sigma,
            blocks = blocks_fit,
            tau = params$tau,
            na_action = "fail"
        )
        if (!is.finite(log_l)) return(1e20)
        -log_l
    }

    control <- list(maxit = as.integer(max_iter), reltol = as.numeric(rel_tol))
    if (isTRUE(verbose)) control$trace <- 1

    opt <- stats::optim(par = par0, fn = objective, method = "BFGS", control = control)
    params_hat <- unpack(opt$par)
    log_l <- -opt$value

    n_mean_params <- n_time + max(n_blocks - 1L, 0L)
    n_cov_uncon <- .count_params_gau(order, n_time, n_blocks = 1L, estimate_mu = FALSE)
    n_params <- n_mean_params + (n_cov_uncon - cinfo$df)

    out <- list(
        mu = params_hat$mu,
        phi = params_hat$phi,
        sigma = params_hat$sigma,
        tau = params_hat$tau,
        log_l = log_l,
        aic = -2 * log_l + 2 * n_params,
        bic = -2 * log_l + n_params * log(n_subjects),
        n_params = as.integer(n_params),
        convergence = opt$convergence,
        n_obs = n_subjects * n_time,
        n_missing = 0L,
        pct_missing = 0,
        missing_pattern = "complete",
        settings = list(
            order = order,
            n_time = n_time,
            n_subjects = n_subjects,
            blocks = blocks_fit,
            block_levels = block_levels,
            estimate_mu = TRUE,
            na_action = "fail",
            stationary = TRUE,
            constraints = cinfo
        )
    )

    class(out) <- c("gau_fit", "list")
    out
}
