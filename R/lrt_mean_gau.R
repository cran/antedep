# File: R/lrt_mean_gau.R
# Likelihood ratio tests for the mean structure (Chapter 7)

#' One-sample test for mean structure under antedependence
#'
#' Tests the null hypothesis that the mean vector equals a specified value
#' mu = mu_0 against the alternative mu != mu_0, under an AD(p) covariance
#' structure. This implements Theorem 7.1 of Zimmerman & Núñez-Antón (2009).
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param mu0 Hypothesized mean vector under the null (length n_time).
#' @param p Antedependence order of the covariance structure. This is the same
#'   order parameter named \code{order} in \code{\link{fit_gau}}.
#' @param order Optional alias for \code{p}. Supply only one of \code{p} or
#'   \code{order}.
#' @param use_modified Logical. If TRUE (default), use the modified test statistic
#'   (formula 7.7) for better small-sample approximation.
#'
#' @return A list with class \code{gau_mean_test} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{test_type}{"one-sample"}
#'   \item{mu0}{Hypothesized mean under null}
#'   \item{mu_hat}{MLE of mean (sample mean)}
#'   \item{statistic}{Test statistic value}
#'   \item{statistic_modified}{Modified test statistic (if use_modified = TRUE)}
#'   \item{df}{Degrees of freedom (n_time)}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{p_value_modified}{P-value from modified test}
#'   \item{order}{Antedependence order used}
#' }
#'
#' @details
#' The test exploits the AD structure to gain power over tests that don't
#' assume any covariance structure. The likelihood ratio test statistic
#' (Theorem 7.1) is:
#' \deqn{N \sum_{i=1}^{n} [\log RSS_i(\mu_0) - \log RSS_i(\hat{\mu})]}
#'
#' where RSS_i(mu) is the residual sum of squares from the regression of
#' Y_i - mu_i on its p predecessors Y_(i-1) - mu_(i-1), ..., Y_(i-p) - mu_(i-p).
#'
#' The test has n degrees of freedom (one for each component of mu).
#'
#' @references
#' Zimmerman, D.L. and Núñez-Antón, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Chapter 7.
#'
#' @seealso \code{\link{test_two_sample_gau}}, \code{\link{test_order_gau}}
#'
#' @examples
#' \donttest{
#' # Simulate data with known mean
#' mu_true <- c(10, 11, 12, 13, 14, 15)
#' y <- simulate_gau(n_subjects = 50, n_time = 6, order = 1, mu = mu_true)
#'
#' # Test if mean is zero
#' test <- test_one_sample_gau(y, mu0 = rep(0, 6), p = 1)
#' print(test)
#'
#' # Test if mean equals true value (should not reject)
#' test2 <- test_one_sample_gau(y, mu0 = mu_true, p = 1)
#' print(test2)
#' }
#'
#' @export
test_one_sample_gau <- function(y, mu0, p = 1L, order = NULL, use_modified = TRUE) {

    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "test_one_sample_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y))) stop("y must contain finite values")

    n <- nrow(y)  # N in the book
    n_time <- ncol(y)  # n in the book

    if (!is.null(order)) {
        if (!missing(p)) {
            stop("Specify only one of 'p' or 'order'.", call. = FALSE)
        }
        p <- order
    }
    p <- as.integer(p)
    if (p < 0 || p >= n_time) stop("p must be in [0, n_time - 1)")

    mu0 <- as.numeric(mu0)
    if (length(mu0) != n_time) stop("mu0 must have length n_time")

    # MLE of mu is the sample mean
    mu_hat <- colMeans(y)

    # Compute RSS under null and alternative
    rss_null <- .rss_vector_gau(y, p, mu = mu0)
    rss_alt <- .rss_vector_gau(y, p, mu = mu_hat)

    # Standard LRT statistic (Theorem 7.1)
    # N * sum [log RSS_i(mu0) - log RSS_i(mu_hat)]
    log_ratios <- log(rss_null) - log(rss_alt)
    statistic <- n * sum(log_ratios)

    # Degrees of freedom = n_time
    df <- n_time

    # P-value
    p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)

    # Modified test statistic (formula 7.7)
    statistic_modified <- NULL
    p_value_modified <- NULL

    if (use_modified) {
        # Modified criterion uses psi weights
        numerator <- 0
        denominator <- 0

        for (i in 1:n_time) {
            p_i <- min(p, i - 1)
            if (rss_null[i] > 0 && rss_alt[i] > 0) {
                log_ratio <- log(rss_null[i]) - log(rss_alt[i])
                numerator <- numerator + log_ratio

                psi_val <- .psi_kenward(1, n - 1 - p_i)
                denominator <- denominator + psi_val
            }
        }

        if (denominator > 0) {
            # Eq. 7.7: leading factor is n_time (number of time points), not N_subjects.
            # Analogous to Eq. 6.9 where Σ diff_i = n_time for the mean test (each diff = 1).
            statistic_modified <- n_time * numerator / denominator
            p_value_modified <- stats::pchisq(statistic_modified, df = df, lower.tail = FALSE)
        }
    }

    out <- list(
        method = "lrt",
        test_type = "one-sample",
        mu0 = mu0,
        mu_hat = mu_hat,
        statistic = statistic,
        statistic_modified = statistic_modified,
        df = df,
        p_value = p_value,
        p_value_modified = p_value_modified,
        order = p,
        n_subjects = n,
        n_time = n_time
    )

    class(out) <- "gau_mean_test"
    out
}


#' Two-sample test for equality of mean profiles under antedependence
#'
#' Tests the null hypothesis that two groups have equal mean profiles
#' mu_1 = mu_2 against the alternative mu_1 != mu_2, assuming a common
#' AD(p) covariance structure. This implements Theorem 7.3 of
#' Zimmerman & Núñez-Antón (2009).
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param blocks Integer vector of length n_subjects indicating group membership
#'   (must contain exactly two unique values, typically 1 and 2).
#' @param p Antedependence order of the common covariance structure. This is the
#'   same order parameter named \code{order} in \code{\link{fit_gau}}.
#' @param order Optional alias for \code{p}. Supply only one of \code{p} or
#'   \code{order}.
#' @param use_modified Logical. If TRUE (default), use modified test statistic.
#'
#' @return A list with class \code{gau_mean_test} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{test_type}{"two-sample"}
#'   \item{mu1_hat}{Estimated mean for group 1}
#'   \item{mu2_hat}{Estimated mean for group 2}
#'   \item{mu_pooled}{Pooled mean estimate under H0}
#'   \item{statistic}{Test statistic value}
#'   \item{statistic_modified}{Modified test statistic}
#'   \item{df}{Degrees of freedom (n_time)}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{p_value_modified}{P-value from modified test}
#'   \item{order}{Antedependence order used}
#' }
#'
#' @details
#' This test is also known as a "profile comparison" test. The likelihood
#' ratio test statistic (Theorem 7.3) compares the pooled RSS (under H0:
#' common mean) to the sum of within-group RSS (under H1: separate means):
#'
#' \deqn{N \sum_{i=1}^{n} [\log RSS_i(\mu) - \log RSS_i(\mu_1, \mu_2)]}
#'
#' where RSS_i(mu) uses a common mean and RSS_i(mu_1, mu_2) uses group-specific
#' means.
#'
#' @references
#' Zimmerman, D.L. and Núñez-Antón, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Chapter 7.
#'
#' @examples
#' \donttest{
#' # Simulate data from two groups with different means
#' n1 <- 30
#' n2 <- 35
#' y1 <- simulate_gau(n1, n_time = 6, order = 1, mu = rep(10, 6))
#' y2 <- simulate_gau(n2, n_time = 6, order = 1, mu = rep(12, 6))
#' y <- rbind(y1, y2)
#' blocks <- c(rep(1, n1), rep(2, n2))
#'
#' # Test equality of profiles
#' test <- test_two_sample_gau(y, blocks, p = 1)
#' print(test)
#' }
#'
#' @export
test_two_sample_gau <- function(y, blocks, p = 1L, order = NULL, use_modified = TRUE) {

    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "test_two_sample_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y))) stop("y must contain finite values")

    n <- nrow(y)
    n_time <- ncol(y)

    blocks <- as.integer(blocks)
    if (length(blocks) != n) stop("blocks must have length nrow(y)")

    unique_blocks <- sort(unique(blocks))
    if (length(unique_blocks) != 2) {
        stop("blocks must contain exactly two unique values for two-sample test")
    }

    # Relabel to 1 and 2
    blocks <- ifelse(blocks == unique_blocks[1], 1L, 2L)

    g1 <- which(blocks == 1)
    g2 <- which(blocks == 2)
    n1 <- length(g1)
    n2 <- length(g2)

    if (!is.null(order)) {
        if (!missing(p)) {
            stop("Specify only one of 'p' or 'order'.", call. = FALSE)
        }
        p <- order
    }
    p <- as.integer(p)
    if (p < 0 || p >= n_time) stop("p must be in [0, n_time - 1)")

    # Group means
    mu1_hat <- colMeans(y[g1, , drop = FALSE])
    mu2_hat <- colMeans(y[g2, , drop = FALSE])

    # Pooled mean (under H0)
    mu_pooled <- colMeans(y)

    # RSS under null (common mean)
    rss_null <- numeric(n_time)
    for (i in 1:n_time) {
        rss_null[i] <- .rss_gau(y, i, p, mu = mu_pooled)
    }

    # RSS under alternative (group-specific means)
    # Sum of within-group RSS
    rss_alt <- numeric(n_time)
    for (i in 1:n_time) {
        rss1 <- .rss_gau(y[g1, , drop = FALSE], i, p, mu = mu1_hat)
        rss2 <- .rss_gau(y[g2, , drop = FALSE], i, p, mu = mu2_hat)
        rss_alt[i] <- rss1 + rss2
    }

    # Standard LRT statistic (Theorem 7.3)
    log_ratios <- log(rss_null) - log(rss_alt)
    statistic <- n * sum(log_ratios)

    # Degrees of freedom = n_time
    df <- n_time

    # P-value
    p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)

    # Modified test statistic (formula 7.11)
    statistic_modified <- NULL
    p_value_modified <- NULL

    if (use_modified) {
        numerator <- 0
        denominator <- 0

        for (i in 1:n_time) {
            p_i <- min(p, i - 1)
            if (rss_null[i] > 0 && rss_alt[i] > 0) {
                log_ratio <- log(rss_null[i]) - log(rss_alt[i])
                numerator <- numerator + log_ratio

                psi_val <- .psi_kenward(1, n - 2 - p_i)
                denominator <- denominator + psi_val
            }
        }

        if (denominator > 0) {
            # Eq. 7.11: leading factor is n_time (number of time points), not N_subjects.
            # Analogous to Eq. 6.9 where Σ diff_i = n_time for the mean test (each diff = 1).
            statistic_modified <- n_time * numerator / denominator
            p_value_modified <- stats::pchisq(statistic_modified, df = df, lower.tail = FALSE)
        }
    }

    out <- list(
        method = "lrt",
        test_type = "two-sample",
        mu1_hat = mu1_hat,
        mu2_hat = mu2_hat,
        mu_pooled = mu_pooled,
        statistic = statistic,
        statistic_modified = statistic_modified,
        df = df,
        p_value = p_value,
        p_value_modified = p_value_modified,
        order = p,
        n_subjects = n,
        n1 = n1,
        n2 = n2,
        n_time = n_time
    )

    class(out) <- "gau_mean_test"
    out
}


#' Test linear hypotheses on the mean under antedependence
#'
#' Tests the null hypothesis C * mu = c for a specified contrast matrix C
#' and vector c, under an AD(p) covariance structure. This implements
#' Theorem 7.2 of Zimmerman & Núñez-Antón (2009).
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param C Contrast matrix with c rows and n_time columns, where c is the
#'   number of contrasts being tested. Rows must be linearly independent.
#' @param c Right-hand side vector of length equal to nrow(C). Default is
#'   a vector of zeros.
#' @param p Antedependence order of the covariance structure. This is the same
#'   order parameter named \code{order} in \code{\link{fit_gau}}.
#'
#' @return A list with class \code{gau_contrast_test} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"wald"}).}
#'   \item{C}{Contrast matrix}
#'   \item{c}{Right-hand side vector}
#'   \item{mu_hat}{Estimated mean vector}
#'   \item{contrast_est}{Estimated value of C * mu}
#'   \item{statistic}{Wald test statistic}
#'   \item{df}{Degrees of freedom (number of contrasts)}
#'   \item{p_value}{P-value from chi-square distribution}
#' }
#'
#' @details
#' The Wald test statistic (Theorem 7.2) is:
#' \deqn{(C\bar{Y} - c)^T (C \hat{\Sigma} C^T)^{-1} (C\bar{Y} - c)}
#'
#' where \eqn{\hat{\Sigma}} is the REML estimator of the covariance matrix
#' under the AD(p) model.
#'
#' Common examples include:
#' \itemize{
#'   \item Testing if mean is constant: C is the first-difference matrix
#'   \item Testing for linear trend: C tests deviations from linearity
#' }
#'
#' @references
#' Zimmerman, D.L. and Núñez-Antón, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Chapter 7.
#'
#' @examples
#' \donttest{
#' y <- simulate_gau(n_subjects = 50, n_time = 5, order = 1)
#'
#' # Test if mean is constant (all differences = 0)
#' # C is 4x5 matrix of first differences
#' C <- matrix(0, nrow = 4, ncol = 5)
#' for (i in 1:4) {
#'   C[i, i] <- 1
#'   C[i, i+1] <- -1
#' }
#' test <- test_contrast_gau(y, C = C, p = 1)
#' print(test)
#' }
#'
#' @export
test_contrast_gau <- function(y, C, c = NULL, p = 1L) {

    if (!is.matrix(y)) y <- as.matrix(y)
    if (!is.matrix(C)) C <- as.matrix(C)
    if (anyNA(y)) {
        stop(
            "test_contrast_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y))) stop("y must contain finite values")

    n <- nrow(y)
    n_time <- ncol(y)

    if (ncol(C) != n_time) {
        stop("C must have n_time columns")
    }

    n_contrasts <- nrow(C)

    if (is.null(c)) {
        c <- rep(0, n_contrasts)
    }
    c <- as.numeric(c)
    if (length(c) != n_contrasts) {
        stop("c must have length equal to nrow(C)")
    }

    p <- as.integer(p)

    # MLE of mean
    mu_hat <- colMeans(y)

    # Estimate covariance matrix under AD(p)
    # For this we need the phi and sigma estimates
    fit <- fit_gau(y, order = p, estimate_mu = FALSE)

    # Construct the AD covariance matrix from fit
    Sigma_hat <- .construct_gau_covariance(fit, n_time)

    # Wald test
    contrast_est <- as.numeric(C %*% mu_hat)
    diff <- contrast_est - c

    # (C Sigma C')^{-1}
    CSigmaC <- C %*% Sigma_hat %*% t(C)
    CSigmaC_inv <- tryCatch(
        solve(CSigmaC),
        error = function(e) {
            warning("Singular contrast covariance matrix")
            NULL
        }
    )

    if (is.null(CSigmaC_inv)) {
        statistic <- NA
        p_value <- NA
    } else {
        # Wald statistic: N * diff' (C Sigma C')^{-1} diff
        statistic <- n * as.numeric(t(diff) %*% CSigmaC_inv %*% diff)
        p_value <- stats::pchisq(statistic, df = n_contrasts, lower.tail = FALSE)
    }

    out <- list(
        method = "wald",
        C = C,
        c = c,
        mu_hat = mu_hat,
        contrast_est = contrast_est,
        statistic = statistic,
        df = n_contrasts,
        p_value = p_value,
        order = p
    )

    class(out) <- "gau_contrast_test"
    out
}


#' Construct AD covariance matrix from fitted model
#'
#' @param fit Fitted AD model from fit_gau().
#' @param n_time Number of time points.
#'
#' @return Covariance matrix.
#'
#' @keywords internal
.construct_gau_covariance <- function(fit, n_time) {
    ord <- fit$settings$order
    sigma <- fit$sigma
    phi <- fit$phi

    # Innovation variances (delta_i^2 in the book)
    delta2 <- sigma^2

    # Build the modified Cholesky factor T and D
    # Sigma^{-1} = T' D^{-1} T
    # where T is unit lower triangular and D is diagonal

    T_mat <- diag(n_time)
    D_vec <- delta2

    if (ord >= 1 && !is.null(phi)) {
        if (is.matrix(phi)) {
            # Order 2
            for (i in 2:n_time) {
                T_mat[i, i-1] <- -phi[1, i]
                if (i >= 3 && ord >= 2) {
                    T_mat[i, i-2] <- -phi[2, i]
                }
            }
        } else {
            # Order 1
            for (i in 2:n_time) {
                T_mat[i, i-1] <- -phi[i]
            }
        }
    }

    # Sigma = T^{-1} D T^{-T}
    T_inv <- solve(T_mat)
    Sigma <- T_inv %*% diag(D_vec) %*% t(T_inv)

    Sigma
}


#' Print method for AD mean test
#'
#' @param x Object of class \code{gau_mean_test}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gau_mean_test <- function(x, ...) {
    if (x$test_type == "one-sample") {
        cat("One-Sample Test for Mean Under AD(", x$order, ")\n", sep = "")
        cat("================================================\n\n")
        cat("H0: mu = mu_0\n")
        cat("H1: mu != mu_0\n\n")
    } else {
        cat("Two-Sample Test for Mean Profiles Under AD(", x$order, ")\n", sep = "")
        cat("======================================================\n\n")
        cat("H0: mu_1 = mu_2\n")
        cat("H1: mu_1 != mu_2\n\n")
        cat("Group sizes: n1 =", x$n1, ", n2 =", x$n2, "\n")
    }

    cat("Sample size:", x$n_subjects, "subjects,", x$n_time, "time points\n\n")

    cat("Standard LRT:\n")
    cat("  Statistic:", round(x$statistic, 4), "\n")
    cat("  df:", x$df, "\n")
    cat("  p-value:", format.pval(x$p_value, digits = 4), "\n")

    if (!is.null(x$statistic_modified)) {
        cat("\nModified LRT:\n")
        cat("  Statistic:", round(x$statistic_modified, 4), "\n")
        cat("  p-value:", format.pval(x$p_value_modified, digits = 4), "\n")
    }

    invisible(x)
}


#' Print method for AD contrast test
#'
#' @param x Object of class \code{gau_contrast_test}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gau_contrast_test <- function(x, ...) {
    cat("Wald Test for Linear Contrasts Under AD(", x$order, ")\n", sep = "")
    cat("===================================================\n\n")

    cat("H0: C * mu = c\n")
    cat("Number of contrasts:", x$df, "\n\n")

    cat("Wald statistic:", round(x$statistic, 4), "\n")
    cat("df:", x$df, "\n")
    cat("p-value:", format.pval(x$p_value, digits = 4), "\n")

    invisible(x)
}
