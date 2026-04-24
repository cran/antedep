#' Likelihood ratio test for antedependence order (INAD data)
#'
#' Performs a likelihood ratio test comparing INAD models of different orders.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param order_null Order under null hypothesis (0 or 1).
#' @param order_alt Order under alternative hypothesis (1 or 2). Must be order_null + 1.
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional integer vector for block effects.
#' @param use_chibar Logical; if TRUE, use chi-bar-square for boundary test.
#' @param weights Optional weights for chi-bar-square mixture.
#' @param fit_null Optional pre-computed null fit.
#' @param fit_alt Optional pre-computed alternative fit.
#' @param ... Additional arguments passed to fit_inad.
#'
#' @return A list with class \code{"test_order_inad"} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{fit_null}{Fitted model under H0}
#'   \item{fit_alt}{Fitted model under H1}
#'   \item{statistic}{Test statistic value}
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{Chi-square p-value}
#'   \item{p_value_chibar}{Chi-bar-square p-value (if \code{use_chibar = TRUE})}
#'   \item{bic_null}{BIC under H0}
#'   \item{bic_alt}{BIC under H1}
#'   \item{bic_selected}{Which model BIC prefers}
#'   \item{table}{Two-row model comparison table}
#'   \item{settings}{Input and derived settings for the test}
#' }
#'
#' @details
#' The test compares nested INAD models of orders \code{order_null} and
#' \code{order_alt = order_null + 1} using:
#' \deqn{\lambda = 2(\ell_{alt} - \ell_{null})}
#' where \eqn{\ell_{null}} and \eqn{\ell_{alt}} are maximized log-likelihoods
#' under the null and alternative models.
#'
#' The default p-value uses the chi-square approximation with degrees of freedom
#' matching the number of additional dependence parameters introduced under the
#' higher-order model. When \code{use_chibar = TRUE}, a chi-bar-square mixture
#' p-value is also reported for boundary-aware inference.
#'
#' Missing-data inputs are supported through the same \code{na_action} options
#' available in \code{\link{fit_inad}}. If \code{y} has missing values and
#' \code{na_action} is not supplied via \code{...}, this function defaults to
#' \code{na_action = "marginalize"}.
#'
#' @references
#' Li, C. and Zimmerman, D.L. (2026). Integer-valued antedependence models for
#' longitudinal count data. \emph{Biostatistics}.
#'
#' @seealso \code{\link{fit_inad}}, \code{\link{bic_order_inad}},
#'   \code{\link{test_stationarity_inad}}
#'
#' @importFrom stats pchisq
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 40,
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
#'   max_iter = 20
#' )
#' out$statistic
#' @export
test_order_inad <- function(y, order_null = 0, order_alt = 1,
                           thinning = "binom", innovation = "pois",
                           blocks = NULL, use_chibar = TRUE, weights = NULL,
                           fit_null = NULL, fit_alt = NULL, ...) {
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
  if (lrt_stat < 0) { warning("LRT negative, setting to 0"); lrt_stat <- 0 }
  
  df <- if (order_null == 0) n_time - 1 else n_time - 2
  p_value <- 1 - pchisq(lrt_stat, df = df)
  p_value_chibar <- if (use_chibar && df > 0) .chi_bar_pvalue(lrt_stat, df, weights) else NA_real_
  
  bic_null <- bic_inad(fit_null, n_subjects); bic_alt <- bic_inad(fit_alt, n_subjects)
  bic_selected <- if (bic_null <= bic_alt) "null" else "alt"
  
  table <- data.frame(model = c(paste0("Order ", order_null), paste0("Order ", order_alt)),
                      order = c(order_null, order_alt), logLik = c(logL_null, logL_alt),
                      n_params = c(.count_params_inad(order_null, n_time, n_blocks, thinning, innovation),
                                   .count_params_inad(order_alt, n_time, n_blocks, thinning, innovation)),
                      BIC = c(bic_null, bic_alt), stringsAsFactors = FALSE)
  
  result <- list(method = "lrt", fit_null = fit_null, fit_alt = fit_alt,
                 statistic = lrt_stat, lrt_stat = lrt_stat, df = df,
                 p_value = p_value, p_value_chibar = p_value_chibar,
                 bic_null = bic_null, bic_alt = bic_alt, bic_selected = bic_selected, table = table,
                 settings = list(order_null = order_null, order_alt = order_alt, thinning = thinning,
                                 innovation = innovation, n_subjects = n_subjects, n_time = n_time,
                                 n_blocks = n_blocks, use_chibar = use_chibar))
  class(result) <- "test_order_inad"
  result
}

#' @export
print.test_order_inad <- function(x, digits = 4, ...) {
  stat <- if (!is.null(x$statistic)) x$statistic else x$lrt_stat
  cat("\nLikelihood Ratio Test for INAD Model Order\n")
  cat("===========================================\n\n")
  cat("H0: Order", x$settings$order_null, " vs H1: Order", x$settings$order_alt, "\n\n")
  print(x$table, row.names = FALSE, digits = digits)
  cat("\nLRT:", round(stat, digits), " df:", x$df, " p-value:", format.pval(x$p_value, digits = digits))
  if (!is.na(x$p_value_chibar)) cat(" p-value(chi-bar):", format.pval(x$p_value_chibar, digits = digits))
  cat("\nBIC selects: Order", ifelse(x$bic_selected == "null", x$settings$order_null, x$settings$order_alt), "\n")
  invisible(x)
}

#' @keywords internal
.chi_bar_pvalue <- function(stat, r, weights = NULL) {
  if (r == 0) return(as.numeric(stat > 0))
  if (is.null(weights)) weights <- choose(r, 0:r) * 2^(-r)
  tails <- vapply(0:r, function(j) if (j == 0) as.numeric(stat > 0) else 1 - pchisq(stat, df = j), numeric(1))
  sum(weights * tails)
}

#' @keywords internal
.count_params_inad <- function(order, n_time, n_blocks, thinning, innovation) {
  n_alpha <- switch(as.character(order), "0" = 0, "1" = n_time - 1, "2" = 2 * (n_time - 2))
  n_theta <- n_time
  n_tau <- max(0, n_blocks - 1)
  n_nb_inno <- if (innovation == "nbinom") 1 else 0
  n_alpha + n_theta + n_tau + n_nb_inno
}

#' BIC-based order selection for INAD models
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
                      n_params = sapply(orders, function(ord) .count_params_inad(ord, n_time, n_blocks, thinning, innovation)),
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
