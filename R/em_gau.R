#' EM algorithm for Gaussian AD model estimation
#'
#' Convenience wrapper around \code{\link{fit_gau}} with
#' \code{na_action = "em"} to provide a parallel entry point to
#' \code{\link{em_inad}}.
#'
#' @param y Numeric matrix (n_subjects x n_time), may contain NA.
#' @param order Integer 0, 1, or 2.
#' @param blocks Optional vector of block membership (length n_subjects).
#' @param estimate_mu Logical, whether to estimate mu (default TRUE).
#' @param max_iter Maximum EM iterations.
#' @param tol EM convergence tolerance.
#' @param verbose Logical, print EM progress.
#' @param ... Additional arguments passed to \code{\link{fit_gau}}.
#'
#' @return An \code{gau_fit} object as returned by \code{\link{fit_gau}}.
#'
#' @details
#' This is an alias-style helper for users who prefer explicit \code{em_*}
#' entry points across model families.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_gau(n_subjects = 35, n_time = 5, order = 1, phi = 0.3)
#' y[sample(length(y), 8)] <- NA
#' fit <- em_gau(y, order = 1, max_iter = 20, tol = 1e-5)
#' fit$settings$na_action
#'
#' @seealso \code{\link{fit_gau}}, \code{\link{em_inad}}, \code{\link{em_cat}},
#'   \code{\link{fit_cat}}
#' @export
em_gau <- function(y, order = 1, blocks = NULL, estimate_mu = TRUE,
                  max_iter = 100, tol = 1e-6, verbose = FALSE, ...) {
  fit_gau(
    y = y,
    order = order,
    blocks = blocks,
    na_action = "em",
    estimate_mu = estimate_mu,
    em_max_iter = max_iter,
    em_tol = tol,
    em_verbose = verbose,
    ...
  )
}
