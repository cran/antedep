# bic_cat.R - BIC computation for categorical antedependence models

#' @rdname bic_gau
#' @export
bic_cat <- function(fit, n_subjects = NULL) {
  if (!inherits(fit, "cat_fit")) {
    stop("fit must be a cat_fit object from fit_cat()")
  }
  
  if (is.null(n_subjects)) {
    n_subjects <- fit$settings$n_subjects
  }
  
  if (is.null(n_subjects) || n_subjects < 1) {
    stop("n_subjects must be a positive integer")
  }
  
  # BIC = -2 * log_l + k * log(N)
  -2 * fit$log_l + fit$n_params * log(n_subjects)
}


#' @rdname bic_gau
#' @export
aic_cat <- function(fit) {
  if (!inherits(fit, "cat_fit")) {
    stop("fit must be a cat_fit object from fit_cat()")
  }
  
  # AIC = -2 * log_l + 2 * k
  -2 * fit$log_l + 2 * fit$n_params
}


#' BIC-Based Order Selection for Categorical AD Models
#'
#' Fits AD models of increasing orders and selects the best by BIC.
#'
#' @param y Integer matrix of categorical data (n_subjects x n_time).
#' @param max_order Maximum order to consider. Default is 2.
#' @param blocks Optional block membership vector.
#' @param homogeneous Whether to use homogeneous parameters across blocks.
#' @param n_categories Number of categories (inferred if NULL).
#' @param criterion Which criterion to use: "bic" (default) or "aic".
#'
#' @return A list containing:
#'   \item{table}{Data frame with order, log_l, n_params, aic, bic}
#'   \item{bic}{Named numeric vector of BIC values by order}
#'   \item{best_order}{Order with lowest criterion value}
#'   \item{criterion}{Criterion used for order selection (\code{"bic"} or \code{"aic"})}
#'   \item{fits}{List of fitted models}
#'
#' @examples
#' \donttest{
#' y <- simulate_cat(100, 5, order = 1, n_categories = 2)
#' result <- bic_order_cat(y, max_order = 2)
#' print(result)
#' print(result$best_order)
#' }
#'
#' @export
bic_order_cat <- function(y, max_order = 2, blocks = NULL, homogeneous = TRUE,
                          n_categories = NULL, criterion = "bic") {
  
  # Validate criterion
  criterion <- match.arg(criterion, c("bic", "aic"))
  
  # Validate data
  validated <- .validate_y_cat(y, n_categories)
  y <- validated$y
  n_time <- ncol(y)
  
  # Determine maximum allowable order
  max_possible <- n_time - 1
  if (max_order > max_possible) {
    warning("max_order reduced to ", max_possible, " (n_time - 1)")
    max_order <- max_possible
  }
  if (max_order > 2) {
    warning("max_order > 2 not supported; setting to 2")
    max_order <- 2
  }
  
  # Fit models for each order
  orders <- 0:max_order
  fits <- vector("list", length(orders))
  names(fits) <- paste0("order_", orders)
  
  for (i in seq_along(orders)) {
    p <- orders[i]
    fits[[i]] <- fit_cat(y, order = p, blocks = blocks, homogeneous = homogeneous,
                         n_categories = n_categories)
  }
  
  # Build comparison table
  table_df <- data.frame(
    order = orders,
    log_l = sapply(fits, function(f) f$log_l),
    n_params = sapply(fits, function(f) f$n_params),
    aic = sapply(fits, function(f) f$aic),
    bic = sapply(fits, function(f) f$bic)
  )
  
  # Find best order
  if (criterion == "bic") {
    best_idx <- which.min(table_df$bic)
  } else {
    best_idx <- which.min(table_df$aic)
  }
  best_order <- orders[best_idx]
  bic_vals <- stats::setNames(table_df$bic, paste0("order_", orders))
  
  out <- list(
    table = table_df,
    bic = bic_vals,
    best_order = best_order,
    criterion = criterion,
    fits = fits
  )
  class(out) <- "bic_order_cat"
  out
}

#' Print Method for Categorical BIC Order Selection
#'
#' @param x Object returned by \code{\link{bic_order_cat}}.
#' @param digits Number of digits used for printed numeric columns.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.bic_order_cat <- function(x, digits = 4, ...) {
  cat("BIC-Based Order Selection for Categorical AD\n")
  cat("============================================\n\n")
  if (!is.null(x$table)) {
    tbl <- x$table
    num_cols <- vapply(tbl, is.numeric, logical(1))
    tbl[num_cols] <- lapply(tbl[num_cols], round, digits = digits)
    print(tbl, row.names = FALSE)
  }
  criterion <- if (!is.null(x$criterion)) toupper(x$criterion) else "BIC"
  cat("\nBest order by", criterion, ":", x$best_order, "\n")
  invisible(x)
}
