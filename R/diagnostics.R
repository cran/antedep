# diagnostics.R - Diagnostic and visualization functions for antedependence models
#
# Functions:
#   - partial_corr()      : Compute intervenor-adjusted partial correlation matrix
#   - plot_prism()        : PRISM (Partial Residual Intervenor Scatterplot Matrix)
#   - plot_profile()      : Profile/spaghetti plot with mean and SD bands
#
# Internal helpers:
#   - .prism_residuals()  : Compute partial residuals for PRISM
# Suppress R CMD check notes for ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("time", "value", "subject", "mean_val", "sd_val"))
}

# ==============================================================================
# PARTIAL CORRELATION MATRIX
# ==============================================================================

#' Compute intervenor-adjusted partial correlation matrix
#'
#' Computes the partial correlation between `Y[i]` and `Y[j]` adjusting for the
#' "intervenor" variables `Y[i+1], ..., Y[j-1]`. Under an antedependence model
#' of order p, partial correlations for |i-j| > p should be approximately zero.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param test Logical; if TRUE, returns significance flags based on approximate
#'   threshold 2/sqrt(n_eff) where n_eff = n_subjects - (lag - 1). Default FALSE.
#' @param n_digits Integer; number of decimal places for rounding. Default 3.
#'
#' @importFrom stats var cor cov lm residuals sd
#' @importFrom graphics par mtext
#'
#' @return A list with components:
#'   \item{correlation}{Matrix with correlations (upper triangle) and variances (diagonal)}
#'   \item{partial_correlation}{Matrix with partial correlations (lower triangle) and variances (diagonal)}
#'   \item{significant}{(If test=TRUE) Matrix flagging significant partial correlations (1 = significant)}
#'   \item{n_subjects}{Number of subjects}
#'   \item{n_time}{Number of time points}
#'
#' @details
#' The intervenor-adjusted partial correlation between `Y[i]` and `Y[j]` (i < j) is
#' computed as the correlation between the residuals from regressing `Y[i]` and `Y[j]`
#' on the intervenor set `Y[i+1], ..., Y[j-1]`.
#'
#' For adjacent time points (|i-j| = 1), the partial correlation equals the
#' ordinary correlation since there are no intervenors.
#'
#' The diagonal of both returned matrices contains variances (not correlations).
#' This keeps scale information available alongside correlation structure.
#'
#' The significance test uses an approximate threshold of 2/sqrt(n_eff), which
#' corresponds roughly to a 95% confidence bound under normality. This is a
#' rough screening tool, not a formal hypothesis test.
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#' pc <- partial_corr(bolus_inad$y, test = TRUE)
#'
#' # View partial correlations (lower triangle)
#' pc$partial_correlation
#'
#' # Extract variances from the diagonal
#' variances <- diag(pc$partial_correlation)
#'
#' # Check which are "significant" (rough screen for AD order)
#' pc$significant
#' }
#'
#' @references
#' Zimmerman, D. L. and Nunez-Anton, V. (2009). Antedependence Models for
#' Longitudinal Data. CRC Press.
#'
#' @seealso \code{\link{plot_prism}} for visual diagnostics
#'
#' @export
partial_corr <- function(y, test = FALSE, n_digits = 3) {

    # Input validation
    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y))) stop("y must contain finite values only")

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    if (n_time < 2) stop("y must have at least 2 time points")
    if (n_subjects < 3) stop("y must have at least 3 subjects")

    # Compute variance and correlation matrices
    var_y <- var(y)
    cor_y <- cor(y)

    # Upper triangle = correlations, diagonal = variances
    CV <- cor_y
    CV[!upper.tri(cor_y, diag = TRUE)] <- 0
    diag(CV) <- diag(var_y)

    # Lower triangle = partial correlations, diagonal = variances
    PC <- diag(diag(var_y))

    for (i in 2:n_time) {
        for (j in 1:(i - 1)) {
            if (i - j == 1) {
                # Adjacent: partial correlation = ordinary correlation
                PC[i, j] <- cor_y[i, j]
            } else {
                # Non-adjacent: condition on intervenors
                y_A <- y[, c(j, i), drop = FALSE]  # Variables of interest
                y_B <- y[, (j + 1):(i - 1), drop = FALSE]  # Intervenors

                s_AA <- var(y_A)
                s_BB <- var(y_B)
                s_AB <- cov(y_A, y_B)

                # Conditional variance: Var(A|B) = Var(A) - Cov(A,B) Var(B)^{-1} Cov(B,A)
                s_BB_inv <- tryCatch(
                    solve(s_BB),
                    error = function(e) {
                        warning("Singular intervenor covariance at lag ", i - j,
                                "; using generalized inverse")
                        MASS::ginv(s_BB)
                    }
                )
                s_AA_B <- s_AA - s_AB %*% s_BB_inv %*% t(s_AB)

                # Partial correlation from conditional covariance
                denom <- sqrt(s_AA_B[1, 1] * s_AA_B[2, 2])
                if (denom > 1e-10) {
                    PC[i, j] <- s_AA_B[1, 2] / denom
                } else {
                    PC[i, j] <- NA
                    warning("Near-zero conditional variance at (", i, ",", j, ")")
                }
            }
        }
    }

    # Round for display
    CV_rounded <- round(CV, n_digits)
    PC_rounded <- round(PC, n_digits)

    result <- list(
        correlation = CV_rounded,
        partial_correlation = PC_rounded,
        n_subjects = n_subjects,
        n_time = n_time
    )

    # Optional significance testing
    if (test) {
        sig_mat <- diag(n_time) * NA  # NA on diagonal
        for (i in 2:n_time) {
            for (j in 1:(i - 1)) {
                lag <- i - j
                n_eff <- n_subjects - (lag - 1)  # Effective sample size
                if (n_eff > 2) {
                    threshold <- 2 / sqrt(n_eff)
                    sig_mat[i, j] <- as.integer(abs(PC[i, j]) >= threshold)
                } else {
                    sig_mat[i, j] <- NA
                }
            }
        }
        result$significant <- sig_mat
    }

    class(result) <- "partial_corr"
    return(result)
}

#' @export
print.partial_corr <- function(x, ...) {
    cat("Partial Correlation Analysis\n")
    cat("----------------------------\n")
    cat("Subjects:", x$n_subjects, "| Time points:", x$n_time, "\n\n")

    cat("Matrix layout:\n")
    cat("  - Upper triangle: ordinary correlations\n")
    cat("  - Lower triangle: intervenor-adjusted partial correlations\n")
    cat("  - Diagonal: variances\n\n")

    # Combine into single display matrix
    display <- x$correlation
    display[lower.tri(display)] <- x$partial_correlation[lower.tri(x$partial_correlation)]
    print(display)

    if (!is.null(x$significant)) {
        cat("\nSignificance flags (lower triangle, 1 = |r| >= 2/sqrt(n_eff)):\n")
        print(x$significant)
    }

    invisible(x)
}


# ==============================================================================
# PRISM PLOTS
# ==============================================================================

#' Compute partial residuals for PRISM plot
#'
#' @param y Data matrix (n_subjects x n_time).
#' @param i First time index.
#' @param j Second time index (j > i + 1).
#'
#' @return List with resid_i and resid_j (partial residuals).
#'
#' @keywords internal
.prism_residuals <- function(y, i, j) {
    if (j <= i + 1) {
        stop("PRISM requires j > i + 1 (need at least one intervenor)")
    }

    y_i <- y[, i]
    y_j <- y[, j]
    intervenors <- y[, (i + 1):(j - 1), drop = FALSE]

    fit_i <- lm(y_i ~ intervenors)
    fit_j <- lm(y_j ~ intervenors)

    list(
        resid_i = residuals(fit_i),
        resid_j = residuals(fit_j)
    )
}


#' PRISM plot (Partial Residual Intervenor Scatterplot Matrix)
#'
#' Creates a matrix of scatterplots for diagnosing antedependence structure.
#' The upper triangle shows ordinary scatterplots of `Y[i]` vs `Y[j]`.
#' The lower triangle shows PRISM plots: residuals from regressing `Y[i]` and `Y[j]`
#' on the intervenor variables `Y[i+1], ..., Y[j-1]`.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns.
#' @param time_labels Optional character vector of time point labels.
#'   Default uses column names or "T1", "T2", etc.
#' @param pch Point character for scatterplots. Default 20 (filled circle).
#' @param cex Point size. Default 0.6.
#' @param col_upper Color for upper triangle plots. Default "steelblue".
#' @param col_lower Color for lower triangle (PRISM) plots. Default "firebrick".
#' @param main Overall title. Default "PRISM Diagnostic Plot".
#'
#' @return Invisibly returns NULL. Called for side effect (plotting).
#'
#' @details
#' Under an antedependence model of order p, the partial correlation between
#' `Y[i]` and `Y[j]` given the intervenors should be zero when |i-j| > p.
#' This means PRISM plots in the lower triangle should show no association
#' for lags greater than p.
#'
#' Interpretation:
#' \itemize{
#'   \item Upper triangle: Shows marginal associations between time points
#'   \item Lower triangle (PRISM): Shows conditional associations after
#'     removing effects of intervenor variables
#'   \item If AD(1) holds: Only the first sub-diagonal of lower triangle
#'     should show association
#'   \item If AD(2) holds: First two sub-diagonals should show association
#' }
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#' plot_prism(bolus_inad$y)
#'
#' # With custom labels
#' plot_prism(bolus_inad$y, time_labels = paste0("Hour ", seq(0, 44, by = 4)))
#' }
#'
#' @references
#' Zimmerman, D. L. and Nunez-Anton, V. (2009). Antedependence Models for
#' Longitudinal Data. CRC Press. Chapter 2.
#'
#' @seealso \code{\link{partial_corr}} for numerical partial correlations
#'
#' @export
plot_prism <- function(y,
                       time_labels = NULL,
                       pch = 20,
                       cex = 0.6,
                       col_upper = "steelblue",
                       col_lower = "firebrick",
                       main = "PRISM Diagnostic Plot") {

    # Input validation
    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y))) stop("y must contain finite values only")

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    if (n_time < 3) stop("PRISM requires at least 3 time points")

    # Time labels
    if (is.null(time_labels)) {
        if (!is.null(colnames(y))) {
            time_labels <- colnames(y)
        } else {
            time_labels <- paste0("T", 1:n_time)
        }
    }
    if (length(time_labels) != n_time) {
        stop("time_labels must have length equal to number of columns in y")
    }

    # Set up plot grid
    n_panels <- n_time - 1
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mfrow = c(n_panels, n_panels),
        mar = c(2, 2, 2, 1),
        oma = c(0, 0, 3, 0))

    for (row in 1:n_panels) {
        for (col in 1:n_panels) {

            if (row <= col) {
                # Upper triangle + diagonal: ordinary scatterplot
                # Plot Y[row] vs Y[col + 1]
                i <- row
                j <- col + 1

                plot(y[, i], y[, j],
                     pch = pch, cex = cex, col = col_upper,
                     xlab = "", ylab = "",
                     main = if (row == 1) time_labels[j] else "")

                if (col == 1) {
                    mtext(time_labels[i], side = 2, line = 2.5, cex = 0.7)
                }

            } else {
                # Lower triangle: PRISM plot
                # Partial residuals of Y[col] and Y[row + 1] given intervenors
                i <- col
                j <- row + 1

                resids <- .prism_residuals(y, i, j)

                plot(resids$resid_i, resids$resid_j,
                     pch = pch, cex = cex, col = col_lower,
                     xlab = "", ylab = "",
                     main = "")

                # Add correlation to plot
                r <- cor(resids$resid_i, resids$resid_j)
                mtext(sprintf("r=%.2f", r), side = 3, line = -1.5,
                      cex = 0.6, col = col_lower)
            }
        }
    }

    # Overall title
    mtext(main, outer = TRUE, cex = 1.2, font = 2)

    invisible(NULL)
}


# ==============================================================================
# PROFILE PLOTS
# ==============================================================================

#' Profile plot (spaghetti plot) for longitudinal data
#'
#' Creates a profile plot showing individual subject trajectories with
#' overlaid mean trajectory and standard deviation bands.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns,
#'   or a data frame with measurements.
#' @param time_points Optional numeric vector of time points for x-axis.
#'   Default uses 1:n_time or attempts to extract from column names.
#' @param blocks Optional integer vector of block memberships for stratified plotting.
#'   If provided, creates separate panels for each block.
#' @param block_labels Optional character vector of labels for blocks.
#' @param title Plot title. Default "Profile Plot".
#' @param xlab X-axis label. Default "Time".
#' @param ylab Y-axis label. Default "Measurement".
#' @param ylim Optional y-axis limits as c(min, max).
#' @param show_sd Logical; if TRUE (default), show +/- 1 SD error bars.
#' @param individual_alpha Alpha (transparency) for individual trajectories. Default 0.3.
#' @param individual_color Color for individual trajectories. Default "grey50".
#' @param mean_color Color for mean trajectory. Default "blue".
#' @param sd_color Color for SD error bars. Default "red".
#' @param mean_lwd Line width for mean trajectory. Default 2.
#'
#' @return A ggplot2 object (invisibly). Called primarily for side effect (plotting).
#'
#' @details
#' This function provides a quick visual summary of longitudinal data showing:
#' \itemize{
#'   \item Individual subject trajectories (light grey lines)
#'   \item Mean trajectory across subjects (bold colored line)
#'   \item +/- 1 standard deviation bands (error bars)
#' }
#'
#' When `blocks` is provided, the plot is faceted by block membership,
#' allowing comparison of trajectories across treatment groups or other strata.
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#'
#' # Basic profile plot
#' plot_profile(bolus_inad$y)
#'
#' # With block stratification
#' plot_profile(bolus_inad$y, blocks = bolus_inad$blocks,
#'              block_labels = c("2mg", "1mg"))
#'
#' # Customized
#' plot_profile(bolus_inad$y,
#'              time_points = seq(0, 44, by = 4),
#'              title = "Bolus Counts Over Time",
#'              xlab = "Hours", ylab = "Count")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_errorbar labs theme_minimal
#'   theme element_text facet_wrap coord_cartesian
#'
#' @export
plot_profile <- function(y,
                         time_points = NULL,
                         blocks = NULL,
                         block_labels = NULL,
                         title = "Profile Plot",
                         xlab = "Time",
                         ylab = "Measurement",
                         ylim = NULL,
                         show_sd = TRUE,
                         individual_alpha = 0.3,
                         individual_color = "grey50",
                         mean_color = "blue",
                         sd_color = "red",
                         mean_lwd = 2) {

    # Check ggplot2 availability
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for plot_profile()")
    }

    # Input validation
    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y[!is.na(y)]))) stop("y must contain finite values (NA allowed)")

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    # Time points
    if (is.null(time_points)) {
        if (!is.null(colnames(y))) {
            # Try to extract numbers from column names
            time_points <- suppressWarnings(
                as.numeric(gsub("[^0-9.-]", "", colnames(y)))
            )
            if (any(is.na(time_points))) {
                time_points <- 1:n_time
            }
        } else {
            time_points <- 1:n_time
        }
    }
    if (length(time_points) != n_time) {
        stop("time_points must have length equal to number of columns in y")
    }

    # Blocks validation
    if (!is.null(blocks)) {
        if (length(blocks) != n_subjects) {
            stop("blocks must have length equal to number of rows in y")
        }
        blocks <- as.factor(blocks)
        if (!is.null(block_labels)) {
            if (length(block_labels) != nlevels(blocks)) {
                stop("block_labels must match number of unique blocks")
            }
            levels(blocks) <- block_labels
        }
    }

    # Convert to long format
    df <- data.frame(
        subject = rep(1:n_subjects, times = n_time),
        time = rep(time_points, each = n_subjects),
        value = as.vector(y)
    )

    if (!is.null(blocks)) {
        df$block <- rep(blocks, times = n_time)
    }

    # Compute summary statistics
    if (!is.null(blocks)) {
        stats <- do.call(rbind, lapply(split(df, list(df$time, df$block)), function(x) {
            data.frame(
                time = x$time[1],
                block = x$block[1],
                mean_val = mean(x$value, na.rm = TRUE),
                sd_val = sd(x$value, na.rm = TRUE)
            )
        }))
    } else {
        stats <- do.call(rbind, lapply(split(df, df$time), function(x) {
            data.frame(
                time = x$time[1],
                mean_val = mean(x$value, na.rm = TRUE),
                sd_val = sd(x$value, na.rm = TRUE)
            )
        }))
    }
    rownames(stats) <- NULL

    # Build plot
    p <- ggplot2::ggplot()

    # Individual trajectories
    if (!is.null(blocks)) {
        p <- p + ggplot2::geom_line(
            data = df,
            ggplot2::aes(x = time, y = value, group = subject),
            color = individual_color,
            alpha = individual_alpha
        )
    } else {
        p <- p + ggplot2::geom_line(
            data = df,
            ggplot2::aes(x = time, y = value, group = subject),
            color = individual_color,
            alpha = individual_alpha
        )
    }

    # Mean trajectory: use version-compatible line width argument.
    mean_line_args <- list(
        data = stats,
        mapping = ggplot2::aes(x = time, y = mean_val),
        color = mean_color
    )
    if (utils::packageVersion("ggplot2") >= "3.4.0") {
        mean_line_args$linewidth <- mean_lwd
    } else {
        mean_line_args$size <- mean_lwd
    }
    p <- p + do.call(ggplot2::geom_line, mean_line_args)

    # SD error bars
    if (show_sd) {
        eb_width <- diff(range(time_points)) * 0.03
        p <- p + ggplot2::geom_errorbar(
            data = stats,
            ggplot2::aes(x = time,
                         ymin = mean_val - sd_val,
                         ymax = mean_val + sd_val),
            width = eb_width,
            color = sd_color
        )
    }

    # Labels and theme
    p <- p +
        ggplot2::labs(title = title, x = xlab, y = ylab) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
        )

    # Y-axis limits
    if (!is.null(ylim)) {
        p <- p + ggplot2::coord_cartesian(ylim = ylim)
    }

    # Facet by block
    if (!is.null(blocks)) {
        p <- p + ggplot2::facet_wrap(~ block)
    }

    print(p)
    invisible(p)
}
