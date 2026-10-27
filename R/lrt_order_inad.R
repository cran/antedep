#' Likelihood Ratio Test for Antedependence Order (INAD Data)
#'
#' Compares INAD models of consecutive orders by their likelihood ratio
#' statistic and BIC.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param order_null Order under null hypothesis (0 or 1).
#' @param order_alt Order under alternative hypothesis (1 or 2). Must be order_null + 1.
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional integer vector for block effects.
#' @param reference Reference distribution for a p-value. \code{"none"} (the
#'   default) reports no p-value, because under the null hypothesis the added
#'   thinning coefficients lie on the boundary of the parameter space and no
#'   calibrated reference has been established in general. \code{"chibar"}
#'   requests a working p-value from the binomial-weight chi-bar-square mixture,
#'   which assumes orthogonal added parameters with block-diagonal efficient
#'   information; see Details.
#' @param weights Optional weights for the chi-bar-square mixture (only with
#'   \code{reference = "chibar"}).
#' @param fit_null Optional pre-computed null fit.
#' @param fit_alt Optional pre-computed alternative fit.
#' @param use_chibar Deprecated. \code{use_chibar = TRUE} is equivalent to
#'   \code{reference = "chibar"}.
#' @param ... Additional arguments passed to fit_inad.
#'
#' @return A list with class \code{c("test_order_inad", "htest")} containing:
#' \describe{
#'   \item{statistic, lrt_stat}{Likelihood ratio statistic}
#'   \item{df}{Number of added thinning parameters}
#'   \item{p.value, p_value}{Working chi-bar-square p-value when
#'     \code{reference = "chibar"}; \code{NA} otherwise}
#'   \item{p_value_chibar}{Same as \code{p_value} (kept for compatibility)}
#'   \item{reference}{The reference used; \code{reference_note} describes its
#'     assumptions}
#'   \item{fit_null, fit_alt}{The two fits}
#'   \item{bic_null, bic_alt, bic_selected}{BIC values and the order BIC prefers}
#'   \item{table}{Two-row model comparison table}
#'   \item{settings}{Input and derived settings for the test}
#' }
#'
#' @details
#' The statistic is \eqn{\lambda = 2(\ell_{alt} - \ell_{null})}, where
#' \eqn{\ell_{null}} and \eqn{\ell_{alt}} are the maximized log-likelihoods.
#' Under the null hypothesis the \code{df} added thinning coefficients are zero,
#' on the boundary of their parameter space, so the ordinary chi-square
#' reference is not valid. The limiting distribution is a chi-bar-square
#' mixture whose weights depend on the efficient information of the added
#' parameters; the binomial weights \eqn{\binom{r}{j}2^{-r}} apply only when
#' these are orthogonal with block-diagonal efficient information. For models
#' with block effects (INADFE) the shared effects couple the time points and
#' these weights have not been established; requesting
#' \code{reference = "chibar"} then gives a warning. Calibrated p-values in
#' general require further work, for example a parametric bootstrap, which is
#' not implemented here. BIC comparisons are always reported.
#'
#' The fits are computed with \code{\link{fit_inad}} unless supplied. For order
#' 2, fits from \code{\link{em_inad}} with \code{newton_polish = TRUE} are
#' usually more precise and can be passed through \code{fit_null} and
#' \code{fit_alt}.
#'
#' Missing-data inputs are supported through the same \code{na_action} options
#' available in \code{\link{fit_inad}}. If \code{y} has missing values and
#' \code{na_action} is not supplied via \code{...}, this function defaults to
#' \code{na_action = "marginalize"}.
#'
#' @references
#' Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. \emph{Journal of the American Statistical Association}, 82,
#' 605-610.
#'
#' @seealso \code{\link{fit_inad}}, \code{\link{em_inad}}, \code{\link{bic_order_inad}},
#'   \code{\link{test_stationarity_inad}}
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 15,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' out <- test_order_inad(
#'   y,
#'   order_null = 0,
#'   order_alt = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   max_iter = 8
#' )
#' out$statistic
#' out$bic_selected
#' # Working chi-bar-square p-value for a baseline model:
#' test_order_inad(y, order_null = 0, order_alt = 1, reference = "chibar",
#'                 fit_null = out$fit_null, fit_alt = out$fit_alt)$p_value
#' @export
test_order_inad <- function(y, order_null = 0, order_alt = 1,
                           thinning = "binom", innovation = "pois",
                           blocks = NULL, reference = c("none", "chibar"), weights = NULL,
                           fit_null = NULL, fit_alt = NULL, use_chibar = NULL, ...) {
  data_name <- deparse(substitute(y))
  if (!is.null(use_chibar)) {
    .Deprecated(msg = "use_chibar is deprecated; use reference = \"chibar\" or \"none\".")
    reference <- if (isTRUE(use_chibar)) "chibar" else "none"
  }
  reference <- match.arg(reference)
  if (!is.matrix(y)) y <- as.matrix(y)
  y_obs <- y[!is.na(y)]
  if (any(y_obs < 0) || any(y_obs != floor(y_obs))) stop("y must contain non-negative integers")
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (!order_null %in% c(0, 1)) stop("order_null must be 0 or 1")
  if (!order_alt %in% c(1, 2)) stop("order_alt must be 1 or 2")
  if (order_alt != order_null + 1) stop("order_alt must be order_null + 1")
  if (order_alt >= n_time) stop("order_alt must be less than n_time")
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  n_blocks <- length(unique(blocks))
  fit_args <- list(...)
  if (anyNA(y) && is.null(fit_args$na_action)) {
    fit_args$na_action <- "marginalize"
  }

  if (is.null(fit_null)) {
    fit_null <- do.call(
      fit_inad,
      c(
        list(y = y, order = order_null, thinning = thinning, innovation = innovation, blocks = blocks),
        fit_args
      )
    )
  }
  if (is.null(fit_alt)) {
    fit_alt <- do.call(
      fit_inad,
      c(
        list(y = y, order = order_alt, thinning = thinning, innovation = innovation, blocks = blocks),
        fit_args
      )
    )
  }

  logL_null <- fit_null$log_l; logL_alt <- fit_alt$log_l
  lrt_stat <- 2 * (logL_alt - logL_null)
  # Nested fits can differ by rounding when the larger model's maximum is on
  # the boundary; only a material shortfall indicates a fitting problem.
  if (lrt_stat < -1e-6 * max(1, abs(logL_null))) {
    warning("The higher-order fit has a smaller log-likelihood; LRT set to 0")
  }
  lrt_stat <- max(lrt_stat, 0)

  df <- if (order_null == 0) n_time - 1 else n_time - 2
  if (reference == "chibar") {
    if (n_blocks > 1L) {
      warning("The binomial chi-bar-square weights have not been established for models ",
              "with block effects; the p-value is a working approximation only.", call. = FALSE)
    }
    p_value <- if (df > 0) .chi_bar_pvalue(lrt_stat, df, weights) else NA_real_
    reference_note <- paste0(
      "working chi-bar-square reference with ", if (is.null(weights)) "binomial" else "user-supplied",
      " weights; assumes orthogonal added parameters with block-diagonal efficient information",
      if (n_blocks > 1L) "; weights not established for models with block effects" else "")
  } else {
    p_value <- NA_real_
    reference_note <- paste0(
      "no calibrated p-value: the added thinning coefficients are on the boundary under H0; ",
      "compare BIC, or use reference = \"chibar\" for a working reference")
  }

  bic_null <- bic_inad(fit_null, n_subjects); bic_alt <- bic_inad(fit_alt, n_subjects)
  bic_selected <- if (bic_null <= bic_alt) "null" else "alt"

  table <- data.frame(model = c(paste0("Order ", order_null), paste0("Order ", order_alt)),
                      order = c(order_null, order_alt), logLik = c(logL_null, logL_alt),
                      n_params = c(.count_params_inad_fit(fit_null), .count_params_inad_fit(fit_alt)),
                      BIC = c(bic_null, bic_alt), stringsAsFactors = FALSE)

  result <- list(
    # htest required slots
    statistic  = stats::setNames(lrt_stat, "LR"),
    parameter  = stats::setNames(df, "added parameters"),
    p.value    = p_value,
    method     = paste0("LRT for INAD Order (H0: Order ", order_null,
                        " vs H1: Order ", order_alt, "); ", reference_note),
    data.name  = data_name,
    # additional slots
    inference  = "lrt",
    fit_null = fit_null, fit_alt = fit_alt,
    lrt_stat = lrt_stat, df = df,
    p_value = p_value, p_value_chibar = p_value,
    reference = reference, reference_note = reference_note,
    bic_null = bic_null, bic_alt = bic_alt, bic_selected = bic_selected, table = table,
    settings = list(order_null = order_null, order_alt = order_alt, thinning = thinning,
                    innovation = innovation, n_subjects = n_subjects, n_time = n_time,
                    n_blocks = n_blocks, reference = reference))
  class(result) <- c("test_order_inad", "htest")
  result
}

#' @export
print.test_order_inad <- function(x, digits = 4, ...) {
  stat <- if (!is.null(x$statistic)) x$statistic else x$lrt_stat
  cat("\nLikelihood Ratio Test for INAD Model Order\n\n")
  cat("H0: Order", x$settings$order_null, " vs H1: Order", x$settings$order_alt, "\n\n")
  print(x$table, row.names = FALSE, digits = digits)
  cat("\nLR statistic:", round(stat, digits), " added parameters:", x$df, "\n")
  if (identical(x$reference, "chibar")) {
    cat("Working p-value (chi-bar-square):", format.pval(x$p_value, digits = digits), "\n")
  }
  if (!is.null(x$reference_note)) cat("Note:", x$reference_note, "\n")
  cat("BIC selects: Order", ifelse(x$bic_selected == "null", x$settings$order_null, x$settings$order_alt), "\n")
  invisible(x)
}

#' Upper-tail probability of a chi-bar-square mixture
#'
#' The j = 0 component is a point mass at zero, so its upper-tail probability
#' P(chi^2_0 >= stat) is 1 when stat <= 0 and 0 when stat > 0.
#'
#' @keywords internal
.chi_bar_pvalue <- function(stat, r, weights = NULL) {
  tail0 <- as.numeric(stat <= 0)
  if (r == 0) return(tail0)
  if (is.null(weights)) weights <- choose(r, 0:r) * 2^(-r)
  tails <- vapply(0:r, function(j) {
    if (j == 0) tail0 else stats::pchisq(stat, df = j, lower.tail = FALSE)
  }, numeric(1))
  sum(weights * tails)
}

#' Free-parameter count of an unconstrained INAD model
#'
#' Matches \code{.count_params_inad_fit()}: order 2 has \code{2 * n_time - 3}
#' thinning parameters (\eqn{\alpha_{2,1}}, then two per time point), and negative
#' binomial innovations carry one size parameter per time point, as fitted by
#' \code{fit_inad()}.
#'
#' @keywords internal
.count_params_inad <- function(order, n_time, n_blocks, thinning, innovation) {
  n_alpha <- switch(as.character(order), "0" = 0, "1" = n_time - 1, "2" = 2 * n_time - 3)
  n_theta <- n_time
  n_tau <- max(0, n_blocks - 1)
  n_nb_inno <- if (innovation == "nbinom") n_time else 0
  n_alpha + n_theta + n_tau + n_nb_inno
}

#' BIC-Based Order Selection for INAD Models
#'
#' Fits INAD models across candidate orders and reports BIC-based selection.
#' @param y Integer matrix.
#' @param max_order Maximum order (1 or 2).
#' @param thinning Thinning operator.
#' @param innovation Innovation distribution.
#' @param blocks Optional block assignments.
#' @param ... Additional arguments.
#' @return A list with class \code{"bic_order_inad"} containing:
#' \describe{
#'   \item{fits}{List of fitted INAD models by candidate order}
#'   \item{bic}{Named numeric vector of BIC values by order}
#'   \item{best_order}{Order with minimum BIC}
#'   \item{table}{Data frame with order, logLik, n_params, and BIC}
#'   \item{settings}{Input and derived settings used for selection}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 30,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.3,
#'   theta = 2
#' )
#' ord <- bic_order_inad(y, max_order = 2, thinning = "binom", innovation = "pois", max_iter = 10)
#' ord$best_order
#' }
#' @export
bic_order_inad <- function(y, max_order = 2, thinning = "binom", innovation = "pois",
                           blocks = NULL, ...) {
  if (!is.matrix(y)) y <- as.matrix(y)
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  n_blocks <- length(unique(blocks))
  max_order <- min(max_order, n_time - 1); orders <- 0:max_order
  fits <- lapply(orders, function(ord) fit_inad(y, order = ord, thinning = thinning, 
                                                 innovation = innovation, blocks = blocks, ...))
  names(fits) <- paste0("order_", orders)
  table <- data.frame(order = orders, logLik = sapply(fits, function(f) f$log_l),
                      n_params = sapply(fits, .count_params_inad_fit),
                      stringsAsFactors = FALSE)
  table$BIC <- -2 * table$logLik + table$n_params * log(n_subjects)
  bic_vals <- stats::setNames(table$BIC, paste0("order_", orders))
  result <- list(fits = fits, bic = bic_vals, table = table, best_order = orders[which.min(table$BIC)],
                 settings = list(thinning = thinning, innovation = innovation, n_subjects = n_subjects, 
                                 n_time = n_time, n_blocks = n_blocks))
  class(result) <- "bic_order_inad"
  result
}

#' @export
print.bic_order_inad <- function(x, digits = 4, ...) {
  cat("\nBIC Model Order Comparison for INAD\n====================================\n\n")
  print(x$table, row.names = FALSE, digits = digits)
  cat("\nBest order by BIC:", x$best_order, "\n")
  invisible(x)
}
