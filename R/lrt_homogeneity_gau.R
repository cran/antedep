# File: R/test_homogeneity_gau.R
# Likelihood ratio test for homogeneity of covariance across groups (Section 6.4)

#' Likelihood ratio test for homogeneity across groups (Gaussian AD data)
#'
#' Tests the null hypothesis that G groups have the same AD(p) covariance
#' structure against the alternative that they have different AD(p) structures.
#' This implements Theorem 6.6 of Zimmerman & Núñez-Antón (2009).
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param blocks Integer vector of length n_subjects indicating group membership.
#' @param p Antedependence order. This is the same order parameter named
#'   \code{order} in \code{\link{fit_gau}}.
#' @param use_modified Logical. If TRUE (default), use modified test statistic
#'   for better small-sample approximation.
#'
#' @return A list with class \code{gau_homogeneity_test} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{statistic}{Test statistic value}
#'   \item{statistic_modified}{Modified test statistic (if use_modified = TRUE)}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{p_value_modified}{P-value from modified test}
#'   \item{G}{Number of groups}
#'   \item{group_sizes}{Sample sizes for each group}
#'   \item{order}{Antedependence order}
#' }
#'
#' @details
#' The test compares:
#' \itemize{
#'   \item H0: All G groups have the same AD(p) covariance matrix Sigma(theta)
#'   \item H1: Groups have different AD(p) covariance matrices Sigma(theta_g)
#' }
#'
#' The likelihood ratio test statistic (Theorem 6.6) involves comparing
#' pooled and within-group RSS values. The degrees of freedom are
#' (G-1)(2n - p)(p + 1)/2.
#'
#' This test is useful for determining whether a common covariance structure
#' can be assumed across treatment groups before performing mean comparisons.
#'
#' @references
#' Zimmerman, D.L. and Núñez-Antón, V. (2009). Antedependence Models for
#' Longitudinal Data. Chapman & Hall/CRC. Section 6.4.
#'
#' Kenward, M.G. (1987). A method for comparing profiles of repeated measurements.
#' Applied Statistics, 36, 296-308.
#'
#' @seealso \code{\link{test_order_gau}}, \code{\link{test_two_sample_gau}}
#'
#' @importFrom stats setNames
#'
#' @examples
#' \donttest{
#' # Simulate data from two groups with same covariance
#' n1 <- 30
#' n2 <- 35
#' y1 <- simulate_gau(n1, n_time = 6, order = 1, phi = 0.5, sigma = 1)
#' y2 <- simulate_gau(n2, n_time = 6, order = 1, phi = 0.5, sigma = 1)
#' y <- rbind(y1, y2)
#' blocks <- c(rep(1, n1), rep(2, n2))
#'
#' # Test homogeneity
#' test <- test_homogeneity_gau(y, blocks, p = 1)
#' print(test)
#' }
#'
#' @export
test_homogeneity_gau <- function(y, blocks, p = 1L, use_modified = TRUE) {

    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y[!is.na(y)]))) stop("y must contain finite values")

    n <- nrow(y)  # Total N
    n_time <- ncol(y)  # n in the book

    blocks <- as.integer(blocks)
    if (length(blocks) != n) stop("blocks must have length nrow(y)")

    unique_blocks <- sort(unique(blocks))
    G <- length(unique_blocks)

    if (G < 2) stop("Need at least 2 groups for homogeneity test")

    p <- as.integer(p)
    if (p < 0 || p >= n_time) stop("p must be in [0, n_time - 1)")

    missing_info <- .validate_missing(y)
    has_missing <- isTRUE(missing_info$has_missing)

    if (has_missing) {
        allowed_patterns <- c("complete", "dropout", "all_missing")
        bad_patterns <- setdiff(unique(missing_info$patterns), allowed_patterns)
        if (length(bad_patterns) > 0) {
            stop(
                paste0(
                    "test_homogeneity_gau with missing data currently supports monotone dropout only ",
                    "(missing values allowed only at the end of each subject profile). ",
                    "Found unsupported patterns: ",
                    paste(bad_patterns, collapse = ", "),
                    "."
                ),
                call. = FALSE
            )
        }
        if (sum(missing_info$patterns != "all_missing") == 0L) {
            stop("At least one subject must have observed data.", call. = FALSE)
        }
    }

    # Relabel blocks to 1, 2, ..., G
    block_map <- setNames(seq_along(unique_blocks), unique_blocks)
    blocks <- as.integer(block_map[as.character(blocks)])

    # Group indices and sizes
    group_indices <- lapply(1:G, function(g) which(blocks == g))
    group_sizes <- sapply(group_indices, length)
    names(group_sizes) <- paste0("g", 1:G)

    # Check minimum group size
    if (any(group_sizes < p + 2)) {
        warning("Some groups have very small sample sizes relative to order p")
    }

    # Compute pooled and within-group RSS
    # Under H0: common covariance, use pooled RSS
    # Under H1: separate covariances, sum within-group RSS

    # For the pooled case, we assume common mean structure (saturated within groups)
    # Actually for homogeneity of covariance, we typically assume means can differ

    # Group-specific means
    group_means <- lapply(1:G, function(g) {
        colMeans(y[group_indices[[g]], , drop = FALSE])
    })

    # Within-group innovation variance estimates
    delta2_g <- matrix(NA_real_, nrow = G, ncol = n_time)
    if (!has_missing) {
        for (g in 1:G) {
            y_g <- y[group_indices[[g]], , drop = FALSE]
            mu_g <- group_means[[g]]
            rss_g <- .rss_vector_gau(y_g, p, mu = mu_g)
            delta2_g[g, ] <- rss_g / group_sizes[g]
        }
    }

    # Pooled innovation variance estimates under H0 (common covariance structure).
    # For complete data: center by group means across all subjects and estimate common phi.
    # For monotone dropout: apply theorem adaptation using subjects complete through time i.
    delta2_pooled <- rep(NA_real_, n_time)
    Nig_mat <- matrix(NA_integer_, nrow = G, ncol = n_time)

    for (i in 1:n_time) {
        p_i <- min(p, i - 1L)

        if (has_missing) {
            keep_i <- rowSums(is.na(y[, 1:i, drop = FALSE])) == 0L
        } else {
            keep_i <- rep(TRUE, n)
        }

        if (!any(keep_i)) {
            next
        }

        y_i <- y[keep_i, 1:i, drop = FALSE]
        blocks_i <- blocks[keep_i]
        N_i <- nrow(y_i)
        # Group means at time i based on subjects complete through i.
        group_means_i <- vector("list", G)
        for (g in 1:G) {
            idx_g <- which(blocks_i == g)
            N_ig <- length(idx_g)
            Nig_mat[g, i] <- N_ig
            if (N_ig > 0) {
                group_means_i[[g]] <- colMeans(y_i[idx_g, , drop = FALSE])
            } else {
                group_means_i[[g]] <- rep(NA_real_, i)
            }
        }

        # Pooled RSS under H0: center by group means and fit one common regression.
        y_centered_i <- y_i
        for (g in 1:G) {
            idx_g <- which(blocks_i == g)
            if (length(idx_g) > 0) {
                y_centered_i[idx_g, ] <- sweep(
                    y_i[idx_g, , drop = FALSE],
                    2,
                    group_means_i[[g]]
                )
            }
        }

        rss_pooled_i <- .rss_gau(y_centered_i, i = i, p = p_i, mu = rep(0, i))
        delta2_pooled[i] <- rss_pooled_i / N_i

        # Within-group RSS under H1.
        for (g in 1:G) {
            idx_g <- which(blocks_i == g)
            N_ig <- length(idx_g)
            if (N_ig <= 0) {
                next
            }
            y_g <- y_i[idx_g, , drop = FALSE]
            mu_g <- group_means_i[[g]]
            rss_g <- .rss_gau(y_g, i = i, p = p_i, mu = mu_g)
            delta2_g[g, i] <- rss_g / N_ig
        }
    }

    # LRT statistic (Theorem 6.6 part a)
    # complete: sum_g N(g) * sum_i [log(RSS_i(p)/N) - log(RSS_ig(p)/N(g))]
    # dropout: replace N and N(g) by N_i and N_i(g) at each time i

    statistic <- 0
    for (i in 1:n_time) {
        if (!is.finite(delta2_pooled[i]) || delta2_pooled[i] <= 0) {
            next
        }
        for (g in 1:G) {
            if (!is.finite(delta2_g[g, i]) || delta2_g[g, i] <= 0) {
                next
            }
            weight <- if (has_missing) Nig_mat[g, i] else group_sizes[g]
            if (!is.finite(weight) || weight <= 0) {
                next
            }
            statistic <- statistic + weight * (log(delta2_pooled[i]) - log(delta2_g[g, i]))
        }
    }

    # Degrees of freedom: (G - 1)(2n - p)(p + 1)/2
    df <- (G - 1) * (2 * n_time - p) * (p + 1) / 2

    # P-value
    p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)

    # Modified test statistic
    statistic_modified <- NULL
    p_value_modified <- NULL

    if (use_modified) {
        # Modified statistic — Eq. (6.13) of Zimmerman & Núñez-Antón (2009):
        #
        #   G̃² = df × G² / Σ_i Σ_g N_ig [ψ(N_i - N_ig - G + 1, N_ig - p_i - 1) - log(N_i / N_ig)]
        #
        # where df = (G-1)(2n-p)(p+1)/2, G² is the standard statistic,
        # N_i  = total subjects observed through time i (= N for complete data),
        # N_ig = subjects in group g observed through time i (= N_g for complete data).
        # With monotone dropout, N_i and N_ig are the time-varying effective counts
        # already stored in Nig_mat.

        denominator <- 0

        for (i in 1:n_time) {
            p_i    <- min(p, i - 1L)
            N_i    <- if (has_missing) sum(Nig_mat[, i], na.rm = TRUE) else n

            for (g in 1:G) {
                N_ig <- if (has_missing) Nig_mat[g, i] else group_sizes[g]
                if (!is.finite(N_ig) || N_ig <= 0L) next

                arg1 <- N_i - N_ig - G + 1L
                arg2 <- N_ig - p_i - 1L
                if (arg1 > 1L && arg2 > 1L) {
                    psi_val  <- .psi_kenward(arg1, arg2)
                    log_term <- log(N_i / N_ig)
                    denominator <- denominator + N_ig * (psi_val - log_term)
                }
            }
        }

        if (denominator > 0) {
            statistic_modified <- df * statistic / denominator
            p_value_modified <- stats::pchisq(statistic_modified, df = df, lower.tail = FALSE)
        }
    }

    out <- list(
        method = "lrt",
        statistic = statistic,
        statistic_modified = statistic_modified,
        df = df,
        p_value = p_value,
        p_value_modified = p_value_modified,
        G = G,
        group_sizes = group_sizes,
        order = p,
        n_time = n_time
    )

    class(out) <- "gau_homogeneity_test"
    out
}


#' Print method for AD homogeneity test
#'
#' @param x Object of class \code{gau_homogeneity_test}.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#' @export
print.gau_homogeneity_test <- function(x, ...) {
    cat("Test for Homogeneity of AD(", x$order, ") Covariance Across Groups\n", sep = "")
    cat("================================================================\n\n")

    cat("H0: Sigma_1 = Sigma_2 = ... = Sigma_G (common AD covariance)\n")
    cat("H1: At least one Sigma_g differs\n\n")

    cat("Number of groups:", x$G, "\n")
    cat("Group sizes:", paste(names(x$group_sizes), "=", x$group_sizes, collapse = ", "), "\n")
    cat("Time points:", x$n_time, "\n\n")

    cat("Standard LRT:\n")
    cat("  Statistic:", round(x$statistic, 4), "\n")
    cat("  df:", round(x$df, 2), "\n")
    cat("  p-value:", format.pval(x$p_value, digits = 4), "\n")

    if (!is.null(x$statistic_modified)) {
        cat("\nModified LRT:\n")
        cat("  Statistic:", round(x$statistic_modified, 4), "\n")
        cat("  p-value:", format.pval(x$p_value_modified, digits = 4), "\n")
    }

    invisible(x)
}
