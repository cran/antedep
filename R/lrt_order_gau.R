# File: R/test_order_gau.R
# Likelihood ratio tests for the order of antedependence (Chapter 6.2)

#' Likelihood ratio test for antedependence order (Gaussian AD data)
#'
#' Tests the null hypothesis that the data follow an AD(p) model against the
#' alternative that they follow an AD(p+q) model, using the likelihood ratio
#' test described in Theorem 6.4 and 6.5 of Zimmerman & Núñez-Antón (2009).
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param p Order under the null hypothesis (default 0). This is the same
#'   antedependence order parameter named \code{order} in \code{\link{fit_gau}}.
#' @param q Order increment under the alternative (default 1, so alternative is AD(p+q)).
#' @param order_null Optional alias for \code{p} (null order).
#' @param order_alt Optional absolute order under the alternative. When supplied,
#'   \code{q} is computed as \code{order_alt - p}.
#' @param mu Optional mean vector. If NULL, the saturated mean (sample means) is used.
#' @param use_modified Logical. If TRUE (default), use the modified test statistic
#'   (formula 6.9) which has better small-sample properties.
#'
#' @return A list with class \code{gau_order_test} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{p}{Order under null hypothesis}
#'   \item{q}{Order increment}
#'   \item{statistic}{Test statistic value}
#'   \item{statistic_modified}{Modified test statistic (if use_modified = TRUE)}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{p_value_modified}{P-value from modified test (if use_modified = TRUE)}
#'   \item{n_subjects}{Number of subjects}
#'   \item{n_time}{Number of time points}
#' }
#'
#' @details
#' The test is based on the intervenor-adjusted sample partial correlations.
#' Under the null hypothesis AD(p), the partial correlations r_(i,i-k|(i-k+1:i-1))
#' should be zero for k > p.
#'
#' The likelihood ratio test statistic (Theorem 6.4) is:
#' \deqn{-N \sum_{j=1}^{q} \sum_{i=p+j+1}^{n} \log(1 - r^2_{i,i-p-j\cdot(i-p-j+1:i-1)})}
#'
#' which is asymptotically chi-square with (2n - 2p - q - 1)(q/2) degrees of freedom.
#'
#' The modified test (formula 6.9) adjusts for small-sample bias using Kenward's
#' (1987) correction.
#'
#' @references
#' Zimmerman, D.L. and Núñez-Antón, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Chapter 6.
#'
#' Kenward, M.G. (1987). A method for comparing profiles of repeated measurements.
#' Applied Statistics, 36, 296-308.
#'
#' @seealso \code{\link{test_one_sample_gau}}, \code{\link{test_homogeneity_gau}}
#'
#' @examples
#' \donttest{
#' # Simulate AD(1) data
#' y <- simulate_gau(n_subjects = 50, n_time = 6, order = 1, phi = 0.5)
#'
#' # Test AD(0) vs AD(1)
#' test01 <- test_order_gau(y, p = 0, q = 1)
#' print(test01)
#'
#' # Test AD(1) vs AD(2)
#' test12 <- test_order_gau(y, p = 1, q = 1)
#' print(test12)
#' }
#'
#' @export
test_order_gau <- function(y, p = 0L, q = 1L, mu = NULL, use_modified = TRUE,
                           order_null = NULL, order_alt = NULL) {

    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) {
        stop(
            "test_order_gau currently supports complete data only. Missing-data AD likelihood-ratio tests are not implemented yet.",
            call. = FALSE
        )
    }
    if (any(!is.finite(y))) stop("y must contain finite values")

    n <- nrow(y)  # Number of subjects (N in the book)
    n_time <- ncol(y)  # Number of time points (n in the book)

    if (!is.null(order_null)) p <- order_null
    if (!is.null(order_alt)) q <- as.integer(order_alt) - as.integer(p)

    p <- as.integer(p)
    q <- as.integer(q)

    if (p < 0) stop("p must be non-negative")
    if (q < 1) stop("q must be at least 1")
    if (p + q >= n_time) stop("p + q must be less than n_time")

    # Use saturated mean if not provided
    if (is.null(mu)) {
        mu <- colMeans(y)
    }
    if (length(mu) != n_time) stop("mu must have length n_time")

    # Compute intervenor-adjusted partial correlations and test statistic
    # Formula (6.5): -N * sum_{j=1}^{q} sum_{i=p+j+1}^{n} log(1 - r^2)

    log_terms <- numeric(0)
    psi_terms <- numeric(0)

    for (j in 1:q) {
        for (i in (p + j + 1):n_time) {
            # Partial correlation r_{i, i-p-j · (i-p-j+1 : i-1)}
            # This is the correlation between Y_i and Y_{i-p-j} adjusting for intervenors
            lag <- p + j
            r <- .partial_corr_gau(y, i, i - lag, mu = mu)

            log_terms <- c(log_terms, log(1 - r^2))

            # For modified test, compute psi weights
            if (use_modified) {
                # psi(p+q)_i from the book for time i
                p_i <- min(p + q, i - 1)
                psi_val <- .psi_kenward(p_i, n - 2 - p_i)
                psi_terms <- c(psi_terms, psi_val)
            }
        }
    }

    # Standard LRT statistic
    statistic <- -n * sum(log_terms)

    # Degrees of freedom: (2n - 2p - q - 1)(q/2)
    # But this is for the book's "n" which is n_time
    df <- (2 * n_time - 2 * p - q - 1) * q / 2

    # P-value
    p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)

    # Modified test statistic (formula 6.9)
    statistic_modified <- NULL
    p_value_modified <- NULL

    if (use_modified) {
        # The modified statistic uses a weighted sum
        # For simplicity, we use the RSS-based formulation (Theorem 6.5)
        rss_p <- .rss_vector_gau(y, p, mu = mu)
        rss_pq <- .rss_vector_gau(y, p + q, mu = mu)

        # Modified statistic (6.9):
        # G~^2 = [sum_i diff_i] * [sum_i log(RSS_i(p)/RSS_i(p+q))] / sum_i psi_i
        # where diff_i = (p+q)_i - p_i  (additional parameters at time i)
        # This is a product of two separate sums (not an inner product).

        sum_diff <- 0
        sum_lr   <- 0
        denominator <- 0

        for (i in (p + 2):n_time) {
            p_i_null <- min(p, i - 1)
            p_i_alt <- min(p + q, i - 1)
            diff_order <- p_i_alt - p_i_null

            if (diff_order > 0 && rss_p[i] > 0 && rss_pq[i] > 0) {
                log_ratio <- log(rss_p[i]) - log(rss_pq[i])
                sum_diff <- sum_diff + diff_order
                sum_lr   <- sum_lr   + log_ratio

                # m = number of mean parameters (1 for saturated within each regression)
                m <- 1
                psi_val <- .psi_kenward(diff_order, n - m - p_i_alt)
                denominator <- denominator + psi_val
            }
        }

        if (denominator > 0) {
            # Eq. 6.9: [Σ diff_i] × [Σ log-ratio_i] / Σ psi_i
            statistic_modified <- sum_diff * sum_lr / denominator
            p_value_modified <- stats::pchisq(statistic_modified, df = df, lower.tail = FALSE)
        }
    }

    out <- list(
        method = "lrt",
        p = p,
        q = q,
        order_null = p,
        order_alt = p + q,
        statistic = statistic,
        statistic_modified = statistic_modified,
        df = df,
        p_value = p_value,
        p_value_modified = p_value_modified,
        n_subjects = n,
        n_time = n_time
    )

    class(out) <- "gau_order_test"
    out
}


#' BIC-based order selection for Gaussian AD models
#'
#' Fits AD models of increasing orders and selects the best by BIC.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param max_order Maximum order to consider.
#' @param ... Additional arguments passed to \code{\link{fit_gau}}.
#'
#' @return A list with class \code{gau_bic_order} containing:
#' \describe{
#'   \item{fits}{List of fitted models}
#'   \item{bic}{BIC values for each order}
#'   \item{best_order}{Order with lowest BIC}
#'   \item{table}{Summary table}
#' }
#'
#' @seealso \code{\link{bic_order_cat}}, \code{\link{bic_order_inad}},
#'   \code{\link{fit_gau}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 80, n_time = 6, order = 1, phi = 0.4)
#' ord <- bic_order_gau(y, max_order = 2)
#' ord$best_order
#' ord$table
#' }
#'
#' @export
bic_order_gau <- function(y, max_order = 2L, ...) {

    if (!is.matrix(y)) y <- as.matrix(y)

    n <- nrow(y)
    n_time <- ncol(y)

    max_order <- min(as.integer(max_order), n_time - 1L)

    fits <- list()
    bic_vals <- numeric(max_order + 1)

    for (ord in 0:max_order) {
        fits[[ord + 1]] <- fit_gau(y, order = ord, ...)
        bic_vals[ord + 1] <- bic_gau(fits[[ord + 1]], n)
    }

    names(fits) <- paste0("order", 0:max_order)
    names(bic_vals) <- paste0("order", 0:max_order)

    best_order <- which.min(bic_vals) - 1L

    tbl <- data.frame(
        order = 0:max_order,
        log_l = sapply(fits, function(f) f$log_l),
        bic = bic_vals,
        selected = (0:max_order) == best_order
    )

    out <- list(
        fits = fits,
        bic = bic_vals,
        best_order = best_order,
        table = tbl
    )

    class(out) <- "gau_bic_order"
    out
}


#' Print method for AD order test
#'
#' @param x Object of class \code{gau_order_test}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gau_order_test <- function(x, ...) {
    cat("Order Test for Antedependence\n")
    cat("==================================================\n\n")

    cat("H0: AD(", x$p, ")\n", sep = "")
    cat("H1: AD(", x$p + x$q, ")\n\n", sep = "")

    cat("Sample size: n =", x$n_subjects, "subjects,", x$n_time, "time points\n\n")

    cat("Standard statistic (LRT):\n")
    cat("  Statistic:", round(x$statistic, 4), "\n")
    cat("  df:", x$df, "\n")
    cat("  p-value:", format.pval(x$p_value, digits = 4), "\n")

    if (!is.null(x$statistic_modified)) {
        cat("\nModified LRT (Kenward correction):\n")
        cat("  Statistic:", round(x$statistic_modified, 4), "\n")
        cat("  p-value:", format.pval(x$p_value_modified, digits = 4), "\n")
    }

    invisible(x)
}


#' Print method for BIC order selection
#'
#' @param x Object of class \code{gau_bic_order}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gau_bic_order <- function(x, ...) {
    cat("BIC-based Order Selection for Gaussian AD\n")
    cat("==========================================\n\n")

    print(x$table, row.names = FALSE)

    cat("\nBest order by BIC:", x$best_order, "\n")

    invisible(x)
}
