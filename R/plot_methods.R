# R/plot_methods.R
# plot() S3 methods for antedependence model fit objects and diagnostics.
#
# Provides:
#   plot.gau_fit       - parameter panel: mean / phi / sigma over time
#   plot.cat_fit       - marginal probability bar chart over time
#   plot.inad_fit      - alpha and theta parameter traces over time
#   plot.partial_corr  - heatmap of the partial correlation matrix

# Suppress R CMD check note for ggplot2 aesthetics
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("x", "y", "fill", "label"))
}

# ===========================================================================
# plot.gau_fit
# ===========================================================================

#' Plot estimated parameters of a Gaussian AD fit
#'
#' Produces a multi-panel base-graphics plot showing the estimated mean
#' (\eqn{\mu}), innovation standard deviation (\eqn{\sigma}), and (for
#' AD(1) or AD(2)) the antedependence coefficients (\eqn{\phi}) over time.
#'
#' @param x A \code{gau_fit} object from \code{\link{fit_gau}}.
#' @param ... Unused; present for S3 consistency.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 40, n_time = 6, order = 1, phi = 0.4)
#' fit <- fit_gau(y, order = 1)
#' plot(fit)
#' }
#'
#' @export
plot.gau_fit <- function(x, ...) {
    ord  <- x$settings$order
    T    <- x$settings$n_time
    tpts <- seq_len(T)

    has_phi <- !is.null(x$phi) && ord >= 1

    n_panels <- if (has_phi) 3L else 2L
    old_par  <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1L, n_panels), mar = c(4, 4, 3, 1))

    # Panel 1: mean
    mu_vals <- as.numeric(x$mu)
    graphics::plot(tpts, mu_vals, type = "b", pch = 16,
                   xlab = "Time", ylab = expression(mu[t]),
                   main = "Estimated Mean", xaxt = "n")
    graphics::axis(1, at = tpts)

    # Panel 2: sigma
    sig_vals <- as.numeric(x$sigma)
    graphics::plot(tpts, sig_vals, type = "b", pch = 16, col = "steelblue",
                   xlab = "Time", ylab = expression(sigma[t]),
                   main = "Innovation Std Dev", xaxt = "n",
                   ylim = c(0, max(sig_vals) * 1.1))
    graphics::axis(1, at = tpts)
    graphics::abline(h = 0, lty = 2, col = "grey70")

    # Panel 3: phi (if any)
    if (has_phi) {
        if (ord == 1) {
            phi_vals <- as.numeric(x$phi)
            # phi[1] is always 0 (reference), plot from t=2
            graphics::plot(2:T, phi_vals[2:T], type = "b", pch = 16,
                           col = "firebrick",
                           xlab = "Time", ylab = expression(phi[t]),
                           main = expression(paste("Antedep. Coeff. (", phi, ")")),
                           xaxt = "n",
                           ylim = c(-1, 1))
            graphics::axis(1, at = 2:T)
            graphics::abline(h = 0, lty = 2, col = "grey70")
        } else {
            phi_mat <- as.matrix(x$phi)   # 2 x T
            graphics::plot(2:T, phi_mat[1, 2:T], type = "b", pch = 16,
                           col = "firebrick",
                           xlab = "Time",
                           ylab = expression(phi[t]),
                           main = expression(paste("Antedep. Coeff. (", phi, ")")),
                           xaxt = "n",
                           ylim = c(-1, 1))
            if (T >= 3)
                graphics::lines(3:T, phi_mat[2, 3:T], type = "b", pch = 17,
                                col = "darkorange")
            graphics::axis(1, at = 2:T)
            graphics::abline(h = 0, lty = 2, col = "grey70")
            graphics::legend("topright",
                             legend = c(expression(phi[1][t]), expression(phi[2][t])),
                             col = c("firebrick", "darkorange"), pch = c(16, 17),
                             lty = 1, cex = 0.8)
        }
    }

    invisible(x)
}


# ===========================================================================
# plot.cat_fit
# ===========================================================================

#' Plot marginal probabilities of a categorical AD fit
#'
#' Displays a stacked bar chart of the estimated marginal category
#' probabilities at each time point.
#'
#' @param x A \code{cat_fit} object from \code{\link{fit_cat}}.
#' @param ... Unused.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_cat(n_subjects = 80, n_time = 5, order = 1, n_categories = 2)
#' fit <- fit_cat(y, order = 1, n_categories = 2)
#' plot(fit)
#' }
#'
#' @export
plot.cat_fit <- function(x, ...) {
    ord    <- x$settings$order
    T      <- x$settings$n_time
    n_cat  <- x$settings$n_categories

    # Extract marginal probabilities at each time point
    marg <- x$marginal

    # Build a matrix: rows = categories, cols = time points
    prob_mat <- matrix(NA_real_, nrow = n_cat, ncol = T)

    get_t1_probs <- function(marg) {
        if (!is.null(marg[["t1"]])) return(as.numeric(marg[["t1"]]))
        if (!is.null(marg[[1]]))    return(as.numeric(marg[[1]]))
        NULL
    }

    t1_probs <- get_t1_probs(marg)
    if (!is.null(t1_probs) && length(t1_probs) == n_cat)
        prob_mat[, 1] <- t1_probs

    for (t in 2:T) {
        key <- paste0("t", t)
        if (!is.null(marg[[key]])) {
            prob_mat[, t] <- as.numeric(marg[[key]])
        } else {
            # Marginalise from transition probs using t1 marginal
            if (t == 2 && !is.null(x$marginal[["t2_given_1to1"]]) &&
                    !is.null(t1_probs)) {
                trans <- x$marginal[["t2_given_1to1"]]
                prob_mat[, 2] <- as.numeric(crossprod(trans, t1_probs))
            } else {
                trans_key <- paste0("t", t)
                if (!is.null(x$transition[[trans_key]]) && !is.null(prob_mat[, t - 1])) {
                    trans <- x$transition[[trans_key]]
                    if (is.matrix(trans))
                        prob_mat[, t] <- as.numeric(crossprod(trans, prob_mat[, t - 1]))
                }
            }
        }
    }

    # If any column is all NA fall back to uniform
    for (t in seq_len(T)) {
        if (any(is.na(prob_mat[, t])))
            prob_mat[, t] <- rep(1 / n_cat, n_cat)
    }

    cat_labels <- paste0("Cat ", seq_len(n_cat))
    col_pal    <- grDevices::hcl.colors(n_cat, palette = "Set2")

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mar = c(4, 4, 3, 7))

    graphics::barplot(prob_mat,
                      names.arg = paste0("T", seq_len(T)),
                      col       = col_pal,
                      xlab      = "Time",
                      ylab      = "Probability",
                      main      = "Estimated Marginal Probabilities",
                      legend.text = cat_labels,
                      args.legend = list(x = "topright", bty = "n",
                                         inset = c(-0.22, 0), xpd = TRUE))

    invisible(x)
}


# ===========================================================================
# plot.inad_fit
# ===========================================================================

#' Plot estimated parameters of an INAD fit
#'
#' Produces a two-panel base-graphics plot of the estimated thinning
#' parameters (\eqn{\alpha}) and innovation means (\eqn{\theta}) over time.
#'
#' @param x An \code{inad_fit} object from \code{\link{fit_inad}}.
#' @param ... Unused.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_inad(n_subjects = 60, n_time = 5, order = 1,
#'                    thinning = "binom", innovation = "pois",
#'                    alpha = 0.3, theta = 2)
#' fit <- fit_inad(y, order = 1, thinning = "binom",
#'                 innovation = "pois", max_iter = 20)
#' plot(fit)
#' }
#'
#' @export
plot.inad_fit <- function(x, ...) {
    T    <- x$settings$n_time
    tpts <- seq_len(T)
    ord  <- x$settings$order

    has_alpha <- !is.null(x$alpha) && ord >= 1

    n_panels <- if (has_alpha) 2L else 1L
    old_par  <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = c(1L, n_panels), mar = c(4, 4, 3, 1))

    # Panel 1: theta
    theta_vals <- as.numeric(x$theta)
    graphics::plot(tpts, theta_vals, type = "b", pch = 16, col = "steelblue",
                   xlab = "Time", ylab = expression(theta[t]),
                   main = "Innovation Mean (theta)", xaxt = "n",
                   ylim = c(0, max(theta_vals) * 1.15))
    graphics::axis(1, at = tpts)
    graphics::abline(h = 0, lty = 2, col = "grey70")

    # Panel 2: alpha
    if (has_alpha) {
        if (ord == 1) {
            alpha_vals <- as.numeric(x$alpha)
            graphics::plot(2:T, alpha_vals[2:T], type = "b", pch = 16,
                           col = "firebrick",
                           xlab = "Time", ylab = expression(alpha[t]),
                           main = "Thinning Parameter (alpha)", xaxt = "n",
                           ylim = c(0, 1))
            graphics::axis(1, at = 2:T)
            graphics::abline(h = c(0, 1), lty = 2, col = "grey70")
        } else {
            alpha_mat <- as.matrix(x$alpha)   # T x 2
            graphics::plot(2:T, alpha_mat[2:T, 1], type = "b", pch = 16,
                           col = "firebrick",
                           xlab = "Time", ylab = expression(alpha[t]),
                           main = "Thinning Parameters (alpha)", xaxt = "n",
                           ylim = c(0, 1))
            if (T >= 3)
                graphics::lines(3:T, alpha_mat[3:T, 2], type = "b", pch = 17,
                                col = "darkorange")
            graphics::axis(1, at = 2:T)
            graphics::abline(h = c(0, 1), lty = 2, col = "grey70")
            graphics::legend("topright",
                             legend = c(expression(alpha[1][t]), expression(alpha[2][t])),
                             col = c("firebrick", "darkorange"), pch = c(16, 17),
                             lty = 1, cex = 0.8)
        }
    }

    invisible(x)
}


# ===========================================================================
# plot.partial_corr
# ===========================================================================

#' Heatmap plot for a partial_corr object
#'
#' Displays the partial correlation matrix as a colour heatmap using base
#' graphics.  Upper triangle shows ordinary correlations; lower triangle
#' shows intervenor-adjusted partial correlations; the diagonal shows
#' standardised variances.
#'
#' @param x A \code{partial_corr} object from \code{\link{partial_corr}}.
#' @param ... Unused.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#' pc <- partial_corr(bolus_inad$y)
#' plot(pc)
#' }
#'
#' @export
plot.partial_corr <- function(x, ...) {
    T <- x$n_time

    # Build combined display matrix (upper = corr, lower = partial corr)
    mat <- x$correlation
    mat[lower.tri(mat)] <- x$partial_correlation[lower.tri(x$partial_correlation)]

    # Standardise diagonal to [-1, 1] range for colour mapping
    diag(mat) <- 0

    # Colour palette centred at 0
    pal <- grDevices::colorRampPalette(c("royalblue", "white", "firebrick"))(101)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mar = c(4, 4, 4, 5))

    # image() expects matrix[col, row] with col increasing left-to-right
    # and row increasing bottom-to-top, so we flip/transpose
    graphics::image(seq_len(T), seq_len(T), t(mat[T:1, ]),
                    zlim = c(-1, 1),
                    col  = pal,
                    xlab = "Time index j",
                    ylab = "Time index i",
                    main = "Correlation / Partial Correlation Matrix",
                    axes = FALSE)
    graphics::axis(1, at = seq_len(T), labels = seq_len(T))
    graphics::axis(2, at = seq_len(T), labels = T:1)

    # Cell labels
    for (i in seq_len(T)) {
        for (j in seq_len(T)) {
            val <- mat[i, j]
            if (is.finite(val))
                graphics::text(j, T + 1 - i,
                               labels = sprintf("%.2f", val),
                               cex = 0.65,
                               col = if (abs(val) > 0.5) "white" else "black")
        }
    }

    # Colour bar legend (approximate)
    legend_vals <- seq(-1, 1, length.out = 101)
    legend_cols <- pal
    graphics::image(x = c(T + 0.6, T + 0.9),
                    y = seq(1, T, length.out = 101),
                    z = matrix(legend_vals, nrow = 1),
                    col = legend_cols, add = TRUE)
    graphics::mtext("Corr.", side = 4, line = 1, cex = 0.8)

    graphics::box()
    invisible(x)
}
