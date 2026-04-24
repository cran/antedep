#' Fit Gaussian antedependence model by maximum likelihood
#'
#' Fits an AD(0), AD(1), or AD(2) model for Gaussian longitudinal data
#' by maximum likelihood. Missing values can be handled by complete-case
#' deletion or by EM (see \code{\link{em_gau}} for an explicit EM wrapper).
#'
#' @details
#' For missing data with \code{na_action = "em"}, AD orders 0 and 1 are the
#' primary production path. AD order 2 is available, but the current EM
#' implementation uses simplified second-order updates and should be treated as
#' provisional for high-stakes inference.
#'
#' For observed-data likelihood evaluation under MAR without fitting, use
#' \code{\link{logL_gau}} with \code{na_action = "marginalize"}. In contrast,
#' \code{fit_gau} handles missingness via complete-case fitting
#' (\code{na_action = "complete"}) or EM (\code{na_action = "em"}).
#'
#' @param y Numeric matrix (n_subjects x n_time). May contain NA.
#' @param order Integer 0, 1, or 2.
#' @param blocks Optional vector of block membership (length n_subjects).
#' @param na_action One of "fail", "complete", or "em".
#' @param estimate_mu Logical, whether to estimate mu (default TRUE).
#' @param em_max_iter Maximum EM iterations (only used when na_action = "em").
#' @param em_tol EM convergence tolerance (only used when na_action = "em").
#' @param em_verbose Logical, print EM progress (only used when na_action = "em").
#' @param ... Passed through to the EM fitter.
#'
#' @return A list with components including mu, phi, sigma, tau, log_l, n_obs, n_missing.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 30, n_time = 5, order = 1, phi = 0.3)
#' fit <- fit_gau(y, order = 1)
#' fit$log_l
#'
#' y_miss <- y
#' y_miss[1, 2] <- NA
#' fit_em <- fit_gau(y_miss, order = 1, na_action = "em", em_max_iter = 20)
#' fit_em$settings$na_action
#'
#' @seealso \code{\link{em_gau}},
#'   \code{\link{fit_cat}}, \code{\link{fit_inad}}
#' @export
fit_gau <- function(y, order = 1, blocks = NULL,
                   na_action = c("fail", "complete", "em"),
                   estimate_mu = TRUE,
                   em_max_iter = 100, em_tol = 1e-6, em_verbose = FALSE, ...) {
    na_action <- match.arg(na_action)
    blocks_input <- blocks

    if (!is.matrix(y)) y <- as.matrix(y)
    storage.mode(y) <- "double"

    if (!order %in% c(0, 1, 2)) {
        stop("order must be 0, 1, or 2.", call. = FALSE)
    }

    if (!is.null(blocks) && length(blocks) != nrow(y)) {
        stop("blocks must have length n_subjects.", call. = FALSE)
    }

    has_missing <- any(is.na(y))

    if (has_missing) {
        if (na_action == "fail") {
            stop("y contains NA values. Use na_action = 'complete' or 'em'.", call. = FALSE)
        }

        if (na_action == "complete") {
            cc <- .extract_complete_cases(y, blocks)
            y <- cc$y
            blocks <- cc$blocks
            has_missing <- FALSE
        }

        if (na_action == "em") {
            block_info <- .normalize_blocks(blocks, nrow(y))
            blocks_id <- if (is.null(blocks_input)) NULL else block_info$blocks_id
            if (order == 2L) {
                warning(
                    "na_action = 'em' with order = 2 uses a provisional implementation; validate results with care.",
                    call. = FALSE
                )
            }
            fit <- .fit_gau_em(
                y = y,
                order = order,
                blocks = blocks_id,
                estimate_mu = estimate_mu,
                max_iter = em_max_iter,
                tol = em_tol,
                verbose = em_verbose,
                ...
            )
            if (!is.null(blocks_input)) {
                fit$settings$block_levels <- block_info$block_levels
            }
            return(fit)
        }
    }

    if (!is.matrix(y)) y <- as.matrix(y)
    n_subjects <- nrow(y)
    n_time <- ncol(y)
    block_levels <- NULL

    if (n_subjects < 2) {
        stop("Need at least 2 subjects to fit the model.", call. = FALSE)
    }

    if (!is.null(blocks)) {
        block_info <- .normalize_blocks(blocks, n_subjects)
        blocks <- block_info$blocks_id
        block_levels <- block_info$block_levels

        k <- length(unique(blocks))
        tau <- numeric(k)
        tau[1] <- 0

        idx_ref <- which(blocks == 1)
        if (length(idx_ref) == 0) {
            stop("Reference block (level 1) has no subjects.", call. = FALSE)
        }

        mu <- if (isTRUE(estimate_mu)) {
            colMeans(y[idx_ref, , drop = FALSE])
        } else {
            rep(0, n_time)
        }

        if (k > 1) {
            for (b in 2:k) {
                idx_b <- which(blocks == b)
                if (length(idx_b) == 0) next
                diffs <- y[idx_b, , drop = FALSE] - matrix(mu, nrow = length(idx_b), ncol = n_time, byrow = TRUE)
                tau[b] <- mean(diffs)
            }
        }

        mu_mat <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
        tau_mat <- matrix(tau[blocks], nrow = n_subjects, ncol = n_time)
        y_center <- y - mu_mat - tau_mat
    } else {
        tau <- 0
        mu <- if (isTRUE(estimate_mu)) colMeans(y) else rep(0, n_time)
        y_center <- y - matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
    }

    if (order == 0) {
        phi <- numeric(0)
        sigma <- apply(y_center, 2, function(z) sqrt(max(mean(z^2), 1e-8)))
    }

    if (order == 1) {
        phi <- numeric(n_time)
        sigma <- numeric(n_time)
        phi[1] <- 0

        sigma[1] <- sqrt(max(mean(y_center[, 1]^2), 1e-8))

        for (t in 2:n_time) {
            x <- y_center[, t - 1]
            yy <- y_center[, t]
            denom <- sum(x^2)

            if (!is.finite(denom) || denom <= 1e-12) {
                phi[t] <- 0
                resid <- yy
            } else {
                phi[t] <- sum(x * yy) / denom
                resid <- yy - phi[t] * x
            }

            sigma[t] <- sqrt(max(mean(resid^2), 1e-8))
        }
    }

    if (order == 2) {
        phi <- matrix(0, nrow = 2, ncol = n_time)
        sigma <- numeric(n_time)

        sigma[1] <- sqrt(max(mean(y_center[, 1]^2), 1e-8))

        if (n_time >= 2) {
            x <- y_center[, 1]
            yy <- y_center[, 2]
            denom <- sum(x^2)

            if (!is.finite(denom) || denom <= 1e-12) {
                phi[1, 2] <- 0
                resid <- yy
            } else {
                phi[1, 2] <- sum(x * yy) / denom
                resid <- yy - phi[1, 2] * x
            }
            sigma[2] <- sqrt(max(mean(resid^2), 1e-8))
        }

        if (n_time >= 3) {
            for (t in 3:n_time) {
                X <- cbind(y_center[, t - 1], y_center[, t - 2])
                yy <- y_center[, t]

                fit <- tryCatch(
                    stats::lm.fit(x = X, y = yy),
                    error = function(e) NULL
                )

                if (is.null(fit) || any(!is.finite(fit$coefficients))) {
                    phi[, t] <- c(0, 0)
                    resid <- yy
                } else {
                    phi[, t] <- fit$coefficients
                    resid <- fit$residuals
                }

                sigma[t] <- sqrt(max(mean(resid^2), 1e-8))
            }
        }
    }

    log_l <- logL_gau(
        y = y,
        order = order,
        mu = mu,
        phi = phi,
        sigma = sigma,
        blocks = blocks,
        tau = tau,
        na_action = "fail"
    )

    n_blocks <- if (is.null(blocks)) 1 else length(unique(blocks))
    n_params <- .count_params_gau(order, n_time, n_blocks)

    result <- list(
        mu = mu,
        phi = phi,
        sigma = sigma,
        tau = tau,
        log_l = log_l,
        aic = -2 * log_l + 2 * n_params,
        bic = -2 * log_l + n_params * log(n_subjects),
        n_params = as.integer(n_params),
        convergence = 0L,
        n_obs = sum(!is.na(y)),
        n_missing = sum(is.na(y)),
        pct_missing = mean(is.na(y)) * 100,
        missing_pattern = "complete",
        settings = list(
            order = order,
            n_time = n_time,
            n_subjects = n_subjects,
            blocks = blocks,
            block_levels = block_levels,
            estimate_mu = estimate_mu,
            na_action = na_action
        )
    )

    class(result) <- c("gau_fit", "list")
    result
}

#' Print method for gau_fit objects
#'
#' @param x A \code{gau_fit} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.gau_fit <- function(x, ...) {
    cat("Gaussian Antedependence Model Fit\n")
    cat("=================================\n\n")

    cat("Order:", x$settings$order, "\n")
    cat("Time points:", x$settings$n_time, "\n")
    cat("Subjects:", x$settings$n_subjects, "\n")
    cat("Missing values:", x$n_missing, "(", round(x$pct_missing, 2), "%)\n")

    cat("\nLog-likelihood:", round(x$log_l, 4), "\n")
    if (!is.null(x$aic)) cat("AIC:", round(x$aic, 4), "\n")
    if (!is.null(x$bic)) cat("BIC:", round(x$bic, 4), "\n")
    if (!is.null(x$n_params)) cat("Parameters:", x$n_params, "\n")
    if (!is.null(x$convergence)) cat("Convergence code:", x$convergence, "\n")

    invisible(x)
}
