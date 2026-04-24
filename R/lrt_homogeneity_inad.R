# File: R/test_homogeneity_inad.R
# Likelihood ratio tests for homogeneity across groups in INAD models
# Based on Section 3.7 of Li & Zimmerman (2026) Biostatistics

#' Likelihood ratio test for homogeneity across groups (INAD data)
#'
#' Tests hypotheses about parameter equality across treatment or grouping
#' factors in integer-valued antedependence models. Implements the homogeneity
#' testing framework from Section 3.7 of Li & Zimmerman (2026).
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param blocks Integer vector of length n_subjects specifying group membership.
#' @param order Antedependence order (0, 1, or 2).
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param test Type of homogeneity test:
#'   \itemize{
#'     \item \code{"all"}: Tests M1 (pooled) vs M3 (fully heterogeneous)
#'     \item \code{"mean"}: Tests M1 (pooled) vs M2 (shared dependence, different means)
#'     \item \code{"dependence"}: Tests M2 (INADFE) vs M3 (fully heterogeneous)
#'   }
#' @param fit_pooled Optional pre-computed pooled fit (M1).
#' @param fit_inadfe Optional pre-computed INADFE fit (M2).
#' @param fit_hetero Optional pre-computed heterogeneous fit (M3).
#' @param ... Additional arguments passed to \code{fit_inad}.
#'
#' @return A list with class \code{"test_homogeneity_inad"} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{statistic}{Test statistic value}
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{test}{Type of test performed}
#'   \item{fit_null}{Fitted model under H0}
#'   \item{fit_alt}{Fitted model under H1}
#'   \item{bic_null}{BIC under H0}
#'   \item{bic_alt}{BIC under H1}
#'   \item{bic_selected}{Which model BIC prefers}
#'   \item{table}{Summary data frame}
#' }
#'
#' @details
#' The function supports three nested model comparisons as described in the paper:
#'
#' \strong{M1 (Pooled)}: All parameters are common across groups. This corresponds
#' to fitting \code{fit_inad(y, blocks = NULL)}.
#'
#' \strong{M2 (INADFE)}: The thinning parameters \eqn{\alpha} are shared across
#' groups, but innovation means differ via block effects \eqn{\tau}. This is the
#' standard INADFE model fitted via \code{fit_inad(y, blocks = blocks)}.
#'
#' \strong{M3 (Fully Heterogeneous)}: Both \eqn{\alpha} and \eqn{\theta}
#' parameters can differ across groups. This is fitted by running separate
#' \code{fit_inad} calls for each group.
#'
#' The three test types correspond to:
#' \itemize{
#'   \item \code{"all"}: H0: M1 vs H1: M3 (complete homogeneity vs complete heterogeneity)
#'   \item \code{"mean"}: H0: M1 vs H1: M2 (test for group differences in means only)
#'   \item \code{"dependence"}: H0: M2 vs H1: M3 (test for group differences in dependence)
#' }
#'
#' Degrees of freedom are computed as the difference in free parameters between
#' the null and alternative models.
#'
#' @references
#' Li, C. and Zimmerman, D.L. (2026). Integer-valued antedependence models for
#' longitudinal count data. \emph{Biostatistics}. Section 3.7.
#'
#' @seealso \code{\link{fit_inad}}, \code{\link{test_order_inad}},
#'   \code{\link{test_stationarity_inad}}
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#' y <- bolus_inad$y
#' blocks <- bolus_inad$blocks
#'
#' # Test for any group differences (M1 vs M3)
#' test_all <- test_homogeneity_inad(y, blocks, order = 1,
#'                                   thinning = "nbinom", innovation = "bell",
#'                                   test = "all")
#' print(test_all)
#'
#' # Test only for mean differences (M1 vs M2)
#' test_mean <- test_homogeneity_inad(y, blocks, order = 1,
#'                                    thinning = "nbinom", innovation = "bell",
#'                                    test = "mean")
#' print(test_mean)
#'
#' # Test for dependence differences given different means (M2 vs M3)
#' test_dep <- test_homogeneity_inad(y, blocks, order = 1,
#'                                   thinning = "nbinom", innovation = "bell",
#'                                   test = "dependence")
#' print(test_dep)
#' }
#'
#' @importFrom stats pchisq
#' @export
test_homogeneity_inad <- function(y, blocks, order = 1,
                                  thinning = "binom", innovation = "pois",
                                  test = c("all", "mean", "dependence"),
                                  fit_pooled = NULL, fit_inadfe = NULL,
                                  fit_hetero = NULL, ...) {
  
  # Validate inputs
  if (!is.matrix(y)) y <- as.matrix(y)
  y_obs <- y[!is.na(y)]
  if (any(y_obs < 0) || any(y_obs != floor(y_obs))) {
    stop("y must contain non-negative integers")
  }
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Validate blocks
  if (missing(blocks) || is.null(blocks)) {
    stop("blocks must be provided for homogeneity testing")
  }
  blocks <- as.integer(blocks)
  if (length(blocks) != n_subjects) {
    stop("blocks must have length equal to nrow(y)")
  }
  
  unique_blocks <- sort(unique(blocks))
  n_blocks <- length(unique_blocks)
  
  if (n_blocks < 2) {
    stop("Need at least 2 groups for homogeneity testing")
  }
  
  # Normalize block labels to 1, 2, ..., G
  block_map <- setNames(seq_along(unique_blocks), unique_blocks)
  blocks_norm <- as.integer(block_map[as.character(blocks)])
  
  # Validate other arguments
  test <- match.arg(test)
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  order <- as.integer(order)
  if (!order %in% c(0L, 1L, 2L)) stop("order must be 0, 1, or 2")
  fit_args <- list(...)
  if (anyNA(y) && is.null(fit_args$na_action)) {
    fit_args$na_action <- "marginalize"
  }
  
  # Group sizes
  group_sizes <- table(blocks_norm)
  group_indices <- lapply(1:n_blocks, function(g) which(blocks_norm == g))
  
  # Fit required models
  # M1: Pooled model (no block effects)
  if (is.null(fit_pooled) && test %in% c("all", "mean")) {
    fit_pooled <- do.call(
      fit_inad,
      c(
        list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = NULL),
        fit_args
      )
    )
  }
  
  # M2: INADFE model (shared dependence, different means via tau)
  if (is.null(fit_inadfe) && test %in% c("mean", "dependence")) {
    fit_inadfe <- do.call(
      fit_inad,
      c(
        list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = blocks_norm),
        fit_args
      )
    )
  }
  
  # M3: Fully heterogeneous (separate fits per group)
  if (is.null(fit_hetero) && test %in% c("all", "dependence")) {
    fit_hetero <- do.call(
      .fit_inad_heterogeneous,
      c(
        list(
          y = y,
          blocks = blocks_norm,
          n_blocks = n_blocks,
          group_indices = group_indices,
          order = order,
          thinning = thinning,
          innovation = innovation
        ),
        fit_args
      )
    )
  }
  
  # Determine null and alternative models based on test type
  if (test == "all") {
    fit_null <- fit_pooled
    fit_alt <- fit_hetero
    null_name <- "M1 (Pooled)"
    alt_name <- "M3 (Heterogeneous)"
  } else if (test == "mean") {
    fit_null <- fit_pooled
    fit_alt <- fit_inadfe
    null_name <- "M1 (Pooled)"
    alt_name <- "M2 (INADFE)"
  } else {  # test == "dependence"
    fit_null <- fit_inadfe
    fit_alt <- fit_hetero
    null_name <- "M2 (INADFE)"
    alt_name <- "M3 (Heterogeneous)"
  }
  
  # Compute LRT statistic
  logL_null <- fit_null$log_l
  logL_alt <- fit_alt$log_l
  lrt_stat <- 2 * (logL_alt - logL_null)
  
  if (lrt_stat < 0) {
    warning("LRT statistic is negative; setting to 0")
    lrt_stat <- 0
  }
  
  # Count parameters
  n_params_null <- .count_params_homogeneity(fit_null, test, n_blocks, 
                                              order, n_time, thinning, innovation,
                                              is_null = TRUE)
  n_params_alt <- .count_params_homogeneity(fit_alt, test, n_blocks, 
                                             order, n_time, thinning, innovation,
                                             is_null = FALSE)
  
  df <- n_params_alt - n_params_null
  
  if (df <= 0) {
    warning("Degrees of freedom <= 0; p-value may not be meaningful")
    df <- max(1, df)
  }
  
  # P-value from chi-square distribution
  p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
  
  # Compute BIC
  bic_null <- -2 * logL_null + n_params_null * log(n_subjects)
  bic_alt <- -2 * logL_alt + n_params_alt * log(n_subjects)
  bic_selected <- if (bic_null <= bic_alt) "null" else "alt"
  
  # Summary table
  table_df <- data.frame(
    model = c(null_name, alt_name),
    logLik = c(logL_null, logL_alt),
    n_params = c(n_params_null, n_params_alt),
    BIC = c(bic_null, bic_alt),
    stringsAsFactors = FALSE
  )
  
  # Assemble output
  result <- list(
    method = "lrt",
    statistic = lrt_stat,
    lrt_stat = lrt_stat,
    df = df,
    p_value = p_value,
    test = test,
    fit_null = fit_null,
    fit_alt = fit_alt,
    bic_null = bic_null,
    bic_alt = bic_alt,
    bic_selected = bic_selected,
    table = table_df,
    settings = list(
      order = order,
      thinning = thinning,
      innovation = innovation,
      n_subjects = n_subjects,
      n_time = n_time,
      n_blocks = n_blocks,
      group_sizes = as.integer(group_sizes)
    )
  )
  
  class(result) <- "test_homogeneity_inad"
  result
}


#' Fit fully heterogeneous INAD model
#'
#' Internal function that fits separate INAD models for each group and
#' combines the log-likelihoods.
#'
#' @param y Data matrix
#' @param blocks Normalized block vector
#' @param n_blocks Number of groups
#' @param group_indices List of subject indices per group
#' @param order AD order
#' @param thinning Thinning type
#' @param innovation Innovation type
#' @param ... Additional arguments passed to fit_inad
#'
#' @return A list mimicking fit_inad output with combined log-likelihood
#'
#' @keywords internal
.fit_inad_heterogeneous <- function(y, blocks, n_blocks, group_indices,
                                     order, thinning, innovation, ...) {
  
  fits <- vector("list", n_blocks)
  total_log_l <- 0
  
  for (g in seq_len(n_blocks)) {
    y_g <- y[group_indices[[g]], , drop = FALSE]
    fits[[g]] <- fit_inad(y_g, order = order, thinning = thinning,
                          innovation = innovation, blocks = NULL, ...)
    total_log_l <- total_log_l + fits[[g]]$log_l
  }
  
  names(fits) <- paste0("group_", seq_len(n_blocks))
  
  # Return a list with combined info
  list(
    log_l = total_log_l,
    fits_by_group = fits,
    settings = list(
      order = order,
      thinning = thinning,
      innovation = innovation,
      n_blocks = n_blocks,
      heterogeneous = TRUE
    )
  )
}


#' Count parameters for homogeneity test
#'
#' @param fit Fitted model
#' @param test Test type
#' @param n_blocks Number of groups
#' @param order AD order
#' @param n_time Number of time points
#' @param thinning Thinning type
#' @param innovation Innovation type
#' @param is_null Whether this is the null model
#'
#' @return Number of free parameters
#'
#' @keywords internal
.count_params_homogeneity <- function(fit, test, n_blocks, order, n_time,
                                       thinning, innovation, is_null) {

  # Base parameter counts (per group, excluding NB innovation size)
  n_alpha <- switch(as.character(order),
                    "0" = 0,
                    "1" = n_time - 1,
                    "2" = 2 * (n_time - 2))
  n_theta <- n_time

  # NB innovation size params: count directly from the fit so this matches
  # what fit_inad() actually produced (length 1 or n_time per group). For
  # heterogeneous (M3) fits, sum across per-group fits.
  count_nb <- function(f) {
    if (innovation != "nbinom") return(0L)
    if (!is.null(f$nb_inno_size)) return(length(f$nb_inno_size))
    if (!is.null(f$fits_by_group)) {
      return(sum(vapply(f$fits_by_group,
                        function(g) length(g$nb_inno_size), integer(1))))
    }
    n_time  # conservative fallback
  }

  n_nb_inno <- count_nb(fit)
  k_single  <- n_alpha + n_theta + n_nb_inno

  if (test == "all") {
    if (is_null) {
      # M1: Pooled - single set of parameters (one shared NB size vector)
      return(k_single)
    } else {
      # M3: Heterogeneous - G sets of alpha+theta, NB sizes already counted
      # across groups via count_nb(fit)
      return(n_blocks * (n_alpha + n_theta) + n_nb_inno)
    }
  } else if (test == "mean") {
    if (is_null) {
      # M1: Pooled
      return(k_single)
    } else {
      # M2: INADFE - shared alpha, theta, and NB size; group-specific tau
      return(n_alpha + n_theta + (n_blocks - 1) + n_nb_inno)
    }
  } else {  # test == "dependence"
    if (is_null) {
      # M2: INADFE
      return(n_alpha + n_theta + (n_blocks - 1) + n_nb_inno)
    } else {
      # M3: Heterogeneous
      return(n_blocks * (n_alpha + n_theta) + n_nb_inno)
    }
  }
}


#' Run all homogeneity tests for INAD
#'
#' Convenience function to run all three homogeneity tests at once and
#' return a summary.
#'
#' @inheritParams test_homogeneity_inad
#'
#' @return A list with class \code{"homogeneity_tests_inad"} containing results
#'   for all three tests and a summary table.
#'
#' @examples
#' \donttest{
#' data("bolus_inad")
#' tests <- run_homogeneity_tests_inad(bolus_inad$y, bolus_inad$blocks,
#'                                      order = 1, thinning = "nbinom",
#'                                      innovation = "bell")
#' print(tests)
#' }
#'
#' @export
run_homogeneity_tests_inad <- function(y, blocks, order = 1,
                                        thinning = "binom", innovation = "pois",
                                        ...) {
  
  # Fit all three models once
  if (!is.matrix(y)) y <- as.matrix(y)
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  fit_args <- list(...)
  if (anyNA(y) && is.null(fit_args$na_action)) {
    fit_args$na_action <- "marginalize"
  }
  
  blocks <- as.integer(blocks)
  unique_blocks <- sort(unique(blocks))
  n_blocks <- length(unique_blocks)
  block_map <- setNames(seq_along(unique_blocks), unique_blocks)
  blocks_norm <- as.integer(block_map[as.character(blocks)])
  group_indices <- lapply(1:n_blocks, function(g) which(blocks_norm == g))
  
  # Fit M1 (Pooled)
  fit_pooled <- do.call(
    fit_inad,
    c(
      list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = NULL),
      fit_args
    )
  )
  
  # Fit M2 (INADFE)
  fit_inadfe <- do.call(
    fit_inad,
    c(
      list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = blocks_norm),
      fit_args
    )
  )
  
  # Fit M3 (Heterogeneous)
  fit_hetero <- do.call(
    .fit_inad_heterogeneous,
    c(
      list(
        y = y,
        blocks = blocks_norm,
        n_blocks = n_blocks,
        group_indices = group_indices,
        order = order,
        thinning = thinning,
        innovation = innovation
      ),
      fit_args
    )
  )
  
  # Run all three tests
  test_all <- test_homogeneity_inad(y, blocks, order = order,
                                    thinning = thinning, innovation = innovation,
                                    test = "all",
                                    fit_pooled = fit_pooled,
                                    fit_hetero = fit_hetero,
                                    na_action = fit_args$na_action)
  
  test_mean <- test_homogeneity_inad(y, blocks, order = order,
                                     thinning = thinning, innovation = innovation,
                                     test = "mean",
                                     fit_pooled = fit_pooled,
                                     fit_inadfe = fit_inadfe,
                                     na_action = fit_args$na_action)
  
  test_dep <- test_homogeneity_inad(y, blocks, order = order,
                                    thinning = thinning, innovation = innovation,
                                    test = "dependence",
                                    fit_inadfe = fit_inadfe,
                                    fit_hetero = fit_hetero,
                                    na_action = fit_args$na_action)
  
  # Create summary table
  summary_table <- data.frame(
    test = c("M1 vs M3 (all)", "M1 vs M2 (mean)", "M2 vs M3 (dependence)"),
    method = c(test_all$method, test_mean$method, test_dep$method),
    statistic = c(test_all$statistic, test_mean$statistic, test_dep$statistic),
    lrt_stat = c(test_all$lrt_stat, test_mean$lrt_stat, test_dep$lrt_stat),
    df = c(test_all$df, test_mean$df, test_dep$df),
    p_value = c(test_all$p_value, test_mean$p_value, test_dep$p_value),
    bic_selected = c(test_all$bic_selected, test_mean$bic_selected, test_dep$bic_selected),
    stringsAsFactors = FALSE
  )
  
  # Model comparison table
  model_table <- data.frame(
    model = c("M1 (Pooled)", "M2 (INADFE)", "M3 (Heterogeneous)"),
    logLik = c(fit_pooled$log_l, fit_inadfe$log_l, fit_hetero$log_l),
    n_params = c(
      .count_params_homogeneity(fit_pooled, "mean", n_blocks, order, n_time, 
                                 thinning, innovation, is_null = TRUE),
      .count_params_homogeneity(fit_inadfe, "dependence", n_blocks, order, n_time, 
                                 thinning, innovation, is_null = TRUE),
      .count_params_homogeneity(fit_hetero, "all", n_blocks, order, n_time, 
                                 thinning, innovation, is_null = FALSE)
    ),
    stringsAsFactors = FALSE
  )
  model_table$BIC <- -2 * model_table$logLik + model_table$n_params * log(n_subjects)
  
  result <- list(
    method = "lrt",
    test_all = test_all,
    test_mean = test_mean,
    test_dependence = test_dep,
    summary_table = summary_table,
    model_table = model_table,
    fits = list(
      pooled = fit_pooled,
      inadfe = fit_inadfe,
      heterogeneous = fit_hetero
    ),
    best_model_bic = model_table$model[which.min(model_table$BIC)],
    settings = list(
      order = order,
      thinning = thinning,
      innovation = innovation,
      n_subjects = n_subjects,
      n_time = n_time,
      n_blocks = n_blocks
    )
  )
  
  class(result) <- "homogeneity_tests_inad"
  result
}


#' Print method for test_homogeneity_inad
#'
#' @param x Object of class \code{test_homogeneity_inad}
#' @param digits Number of digits for printing
#' @param ... Unused
#' @return Invisibly returns \code{x}.
#' @export
print.test_homogeneity_inad <- function(x, digits = 4, ...) {
  stat <- if (!is.null(x$statistic)) x$statistic else x$lrt_stat
  cat("\nLikelihood Ratio Test for Homogeneity (INAD)\n")
  cat("=============================================\n\n")
  
  test_desc <- switch(x$test,
    "all" = "H0: M1 (Pooled) vs H1: M3 (Heterogeneous)",
    "mean" = "H0: M1 (Pooled) vs H1: M2 (INADFE - different means)",
    "dependence" = "H0: M2 (INADFE) vs H1: M3 (Heterogeneous - different dependence)"
  )
  cat(test_desc, "\n\n")
  
  cat("Number of groups:", x$settings$n_blocks, "\n")
  cat("Group sizes:", paste(x$settings$group_sizes, collapse = ", "), "\n")
  cat("Order:", x$settings$order, "\n")
  cat("Thinning:", x$settings$thinning, " Innovation:", x$settings$innovation, "\n\n")
  
  print(x$table, row.names = FALSE, digits = digits)
  
  cat("\nLRT statistic:", round(stat, digits), "\n")
  cat("Degrees of freedom:", x$df, "\n")
  cat("P-value:", format.pval(x$p_value, digits = digits), "\n")
  cat("BIC selects:", ifelse(x$bic_selected == "null", 
                              x$table$model[1], x$table$model[2]), "\n")
  
  invisible(x)
}


#' Print method for homogeneity_tests_inad
#'
#' @param x Object of class \code{homogeneity_tests_inad}
#' @param digits Number of digits for printing
#' @param ... Unused
#' @return Invisibly returns \code{x}.
#' @export
print.homogeneity_tests_inad <- function(x, digits = 4, ...) {
  cat("\nHomogeneity Tests for INAD Models\n")
  cat("==================================\n\n")
  
  cat("Settings: Order", x$settings$order, ",", 
      x$settings$thinning, "thinning,", x$settings$innovation, "innovations\n")
  cat("Groups:", x$settings$n_blocks, "  N:", x$settings$n_subjects, 
      "  T:", x$settings$n_time, "\n\n")
  
  cat("Model Comparison:\n")
  print(x$model_table, row.names = FALSE, digits = digits)
  cat("\n")
  
  cat("Test Results:\n")
  print(x$summary_table, row.names = FALSE, digits = digits)
  cat("\n")
  
  cat("Best model by BIC:", x$best_model_bic, "\n")
  
  invisible(x)
}
