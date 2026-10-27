#' Likelihood Ratio Test for Time-Invariance (INAD Data)
#'
#' Tests whether time-varying INAD parameters can be constrained to be constant
#' over time.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param order Model order (1 or 2).
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional integer vector for block effects.
#' @param constrain Which parameters to constrain to be constant over time. A
#'   character vector built from \code{"alpha"}, \code{"alpha1"},
#'   \code{"alpha2"} (order 2 only), \code{"theta"}, and, for
#'   \code{innovation = "nbinom"}, \code{"kappa"} (alias
#'   \code{"nb_inno_size"}). The shortcuts \code{"innovation"} (\code{theta}
#'   and, for negative binomial innovations, \code{kappa}), \code{"both"}
#'   (\code{alpha} and \code{theta}), and \code{"all"} (every time-varying
#'   parameter, including \code{kappa}) are also accepted.
#' @param fit_unconstrained Optional pre-computed unconstrained fit.
#' @param verbose Logical; if TRUE, print progress.
#' @param ... Additional arguments passed to \code{\link{fit_inad}} for the
#'   unconstrained fit (for example \code{max_iter}, \code{tol},
#'   \code{na_action}, \code{nb_inno_size_ub}).
#'
#' @return A list with class \code{"test_stationarity_inad"} containing:
#' \describe{
#'   \item{method}{Inference method used (\code{"lrt"}).}
#'   \item{fit_unconstrained}{Unconstrained INAD fit}
#'   \item{fit_constrained}{Constrained INAD fit, an \code{inad_fit} object}
#'   \item{constraint}{Human-readable null constraint description}
#'   \item{statistic}{Test statistic value}
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{Chi-square p-value}
#'   \item{bic_unconstrained}{BIC of unconstrained model}
#'   \item{bic_constrained}{BIC of constrained model}
#'   \item{bic_selected}{Model selected by BIC}
#'   \item{table}{Two-row model summary table}
#' }
#'
#' @details
#' The constrained model is refitted by maximizing the observed-data
#' log-likelihood jointly over every parameter that remains free under the
#' null hypothesis, including the negative binomial innovation sizes when they
#' are not themselves constrained. Fitting alternates blockwise BOBYQA updates
#' with a joint BOBYQA refinement until the log-likelihood stabilizes.
#'
#' For order 2, \code{"alpha1"} constrains the lag-1 coefficients
#' \eqn{\alpha_{i,i-1}}, \eqn{i = 2, \ldots, n}, to a common value, and
#' \code{"alpha2"} constrains the lag-2 coefficients \eqn{\alpha_{i,i-2}},
#' \eqn{i = 3, \ldots, n}.
#'
#' Degrees of freedom are the number of equality constraints imposed under the
#' null model relative to the unconstrained model. Because the constrained
#' parameters are interior under the null, the usual chi-square reference
#' applies.
#'
#' Missing-data inputs are supported through the same \code{na_action} options
#' available in \code{\link{fit_inad}}. If \code{y} has missing values and
#' \code{na_action} is not supplied via \code{...}, this function defaults to
#' \code{na_action = "marginalize"}.
#'
#' @references
#' Li, C. and Zimmerman, D.L. (2026). Integer-valued antedependence models for
#' longitudinal count data. Manuscript under review.
#'
#' @seealso \code{\link{run_stationarity_tests_inad}},
#'   \code{\link{test_order_inad}}, \code{\link{fit_inad}}
#'
#' @importFrom stats pchisq optim
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
#' out <- test_stationarity_inad(
#'   y,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   constrain = "both",
#'   max_iter = 8
#' )
#' out$p_value
#' @export
test_stationarity_inad <- function(y, order = 1, thinning = "binom", innovation = "pois",
                                  blocks = NULL, constrain = "both",
                                  fit_unconstrained = NULL, verbose = FALSE, ...) {
  data_name <- deparse(substitute(y))
  if (!is.matrix(y)) y <- as.matrix(y)
  y_obs <- y[!is.na(y)]
  if (any(y_obs < 0) || any(y_obs != floor(y_obs))) stop("y must contain non-negative integers")
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (!order %in% c(1, 2)) stop("order must be 1 or 2 for stationarity testing")
  if (order >= n_time) stop("order must be less than n_time")
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  n_blocks <- length(unique(blocks))
  fit_args <- list(...)
  if (anyNA(y) && is.null(fit_args$na_action)) {
    fit_args$na_action <- "marginalize"
  }
  na_action_inference <- if (!is.null(fit_args$na_action)) fit_args$na_action else "fail"
  nb_inno_size_ub <- if (!is.null(fit_args$nb_inno_size_ub)) fit_args$nb_inno_size_ub else 50

  constraint_info <- .parse_constraint(constrain, order, n_time, innovation)

  if (is.null(fit_unconstrained)) {
    if (verbose) cat("Fitting unconstrained model...\n")
    fit_unconstrained <- do.call(
      fit_inad,
      c(
        list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = blocks),
        fit_args
      )
    )
  }
  if (verbose) cat("Fitting constrained model...\n")
  fit_constrained <- .fit_constrained(
    y, order, thinning, innovation, blocks, constraint_info, fit_unconstrained,
    na_action = na_action_inference, nb_inno_size_ub = nb_inno_size_ub,
    verbose = verbose
  )

  logL_uncon <- fit_unconstrained$log_l; logL_con <- fit_constrained$log_l
  lrt_stat <- 2 * (logL_uncon - logL_con)
  # A constrained fit can exceed an unconstrained one only by rounding, unless
  # the unconstrained fit is not at its maximum.
  if (lrt_stat < -1e-6 * max(1, abs(logL_uncon))) {
    warning("The constrained fit has a larger log-likelihood; LRT set to 0")
  }
  lrt_stat <- max(lrt_stat, 0)
  df <- constraint_info$df
  p_value <- pchisq(lrt_stat, df = df, lower.tail = FALSE)

  n_params_uncon <- .count_params_inad_fit(fit_unconstrained)
  n_params_con <- n_params_uncon - df
  bic_uncon <- -2 * logL_uncon + n_params_uncon * log(n_subjects)
  bic_con <- -2 * logL_con + n_params_con * log(n_subjects)
  bic_selected <- if (bic_uncon <= bic_con) "unconstrained" else "constrained"

  table <- data.frame(model = c("Unconstrained", paste("Constrained:", constraint_info$description)),
                      logLik = c(logL_uncon, logL_con),
                      n_params = c(n_params_uncon, n_params_con),
                      BIC = c(bic_uncon, bic_con), stringsAsFactors = FALSE)

  result <- list(
    statistic  = stats::setNames(lrt_stat, "chi-squared"),
    parameter  = stats::setNames(df, "df"),
    p.value    = p_value,
    method     = paste0("LRT for INAD Stationarity (", constraint_info$description, ")"),
    data.name  = data_name,
    inference  = "lrt",
    fit_unconstrained = fit_unconstrained, fit_constrained = fit_constrained,
    constraint = constraint_info$description, lrt_stat = lrt_stat, df = df, p_value = p_value,
    bic_unconstrained = bic_uncon, bic_constrained = bic_con, bic_selected = bic_selected,
    estimates = list(unconstrained = list(alpha = fit_unconstrained$alpha, theta = fit_unconstrained$theta,
                                          tau = fit_unconstrained$tau,
                                          nb_inno_size = fit_unconstrained$nb_inno_size),
                     constrained = list(alpha = fit_constrained$alpha, theta = fit_constrained$theta,
                                        tau = fit_constrained$tau,
                                        nb_inno_size = fit_constrained$nb_inno_size)),
    table = table, settings = list(order = order, thinning = thinning, innovation = innovation,
                                   n_subjects = n_subjects, n_time = n_time, n_blocks = n_blocks, constrain = constrain,
                                   na_action = na_action_inference))
  class(result) <- c("test_stationarity_inad", "htest")
  result
}

#' @export
print.test_stationarity_inad <- function(x, digits = 4, ...) {
  stat <- if (!is.null(x$statistic)) x$statistic else x$lrt_stat
  cat("\nLikelihood Ratio Test for INAD Stationarity\n============================================\n\n")
  cat("Constraint: H0:", x$constraint, "\n\n")
  print(x$table, row.names = FALSE, digits = digits)
  cat("\nLRT:", round(stat, digits), " df:", x$df, " p-value:", format.pval(x$p_value, digits = digits))
  cat("\nBIC selects:", x$bic_selected, "\n")
  invisible(x)
}

#' Parse a time-invariance constraint specification
#'
#' @keywords internal
.parse_constraint <- function(constrain, order, n_time, innovation = "pois") {
  if (is.list(constrain) && !is.null(constrain$df)) return(constrain)
  tokens <- tolower(as.character(constrain))
  tokens[tokens == "nb_inno_size"] <- "kappa"
  is_nb <- identical(innovation, "nbinom")
  valid <- c("alpha", "alpha1", "alpha2", "theta", "kappa", "innovation", "both", "all")
  bad <- setdiff(tokens, valid)
  if (length(bad) > 0L) {
    stop("Invalid constraint: ", paste(bad, collapse = ", "),
         ". Use any of: ", paste(valid, collapse = ", "), call. = FALSE)
  }

  alpha1_const <- any(tokens %in% c("alpha", "alpha1", "both", "all"))
  alpha2_const <- order == 2 && any(tokens %in% c("alpha", "alpha2", "both", "all"))
  theta_const <- any(tokens %in% c("theta", "innovation", "both", "all"))
  kappa_const <- is_nb && any(tokens %in% c("kappa", "innovation", "all"))

  if (order == 1 && "alpha2" %in% tokens) stop("'alpha2' is only available for order 2")
  if (!is_nb && "kappa" %in% tokens) stop("'kappa' requires innovation = 'nbinom'")

  # Lag-1 coefficients alpha_{i,i-1}, i = 2..n; lag-2 coefficients alpha_{i,i-2}, i = 3..n.
  df <- (if (alpha1_const) n_time - 2 else 0) +
    (if (alpha2_const) n_time - 3 else 0) +
    (if (theta_const) n_time - 1 else 0) +
    (if (kappa_const) n_time - 1 else 0)
  if (df <= 0) stop("The requested constraint imposes no restriction")

  lab_alpha1 <- if (order == 1) "alpha constant" else "alpha1 constant"
  desc <- paste(c(if (alpha1_const) lab_alpha1,
                  if (alpha2_const) "alpha2 constant",
                  if (theta_const) "theta constant",
                  if (kappa_const) "kappa constant"), collapse = ", ")

  list(order = order, alpha_const = alpha1_const, alpha1_const = alpha1_const,
       alpha2_const = alpha2_const, theta_const = theta_const, kappa_const = kappa_const,
       df = df, description = desc)
}

#' Fit an INAD model under time-invariance constraints
#'
#' Maximizes the observed-data log-likelihood jointly over all parameters left
#' free by the constraint. With complete data it first takes damped Newton
#' steps with the exact score and information from pooled starting values
#' (\code{settings$path = "newton"}); if the first-order conditions do not hold
#' at the result, it continues with blockwise BOBYQA updates and a joint BOBYQA
#' refinement, repeated until the log-likelihood stabilizes
#' (\code{settings$path = "search"}).
#'
#' @return An \code{inad_fit} object whose \code{settings$constraint_df} records
#'   the number of equality constraints, so that \code{n_params}, \code{aic},
#'   and \code{bic} count only the free parameters.
#' @keywords internal
.fit_constrained <- function(y, order, thinning, innovation, blocks, cinfo, init_from = NULL,
                             na_action = "fail", nb_inno_size_ub = 50,
                             max_rounds = 25, tol = 1e-8, verbose = FALSE, ...) {
  if (!is.matrix(y)) y <- as.matrix(y)
  n <- nrow(y); N <- ncol(y)
  cinfo <- .parse_constraint(cinfo, order, N, innovation)
  block_info <- .normalize_blocks(blocks, n)
  blocks_id <- block_info$blocks_id
  B <- block_info$n_blocks
  is_nb <- innovation == "nbinom"
  eps <- 1e-10
  # Family-specific limits for the thinning coefficients (alpha <= 1 only for
  # binomial thinning); a coefficient common to all occasions gets its own limit.
  ub_a1 <- if (cinfo$alpha1_const) .inad_alpha_upper(y, 1, thinning, common = TRUE)
           else .inad_alpha_upper(y, 1, thinning)[2:N]
  ub_a2 <- if (order == 2) {
    if (cinfo$alpha2_const) .inad_alpha_upper(y, 2, thinning, common = TRUE)
    else .inad_alpha_upper(y, 2, thinning)[3:N]
  } else NULL

  # Starting values from the unconstrained fit when available.
  theta_default <- colMeans(y, na.rm = TRUE)
  theta_default[!is.finite(theta_default) | theta_default <= 0] <- 1
  if (innovation == "bell") theta_default <- pmax(.bell_theta_from_mean(theta_default), 1e-3)
  if (!is.null(init_from$alpha)) {
    au <- .unpack_alpha(order, init_from$alpha, N)
    a1_0 <- au$a1; a2_0 <- au$a2
  } else {
    a1_0 <- c(0, rep(0.3, N - 1))
    a2_0 <- if (order == 2) c(0, 0, rep(0.1, N - 2)) else NULL
  }
  theta0 <- if (!is.null(init_from$theta)) as.numeric(init_from$theta) else theta_default
  if (length(theta0) == 1L) theta0 <- rep(theta0, N)
  kappa0 <- NULL
  if (is_nb) {
    kappa0 <- if (!is.null(init_from$nb_inno_size)) as.numeric(init_from$nb_inno_size) else rep(1, N)
    if (length(kappa0) == 1L) kappa0 <- rep(kappa0, N)
    kappa0 <- pmin(pmax(kappa0, 1e-6), nb_inno_size_ub)
  }
  tau0 <- rep(0, B)
  if (B > 1L && !is.null(init_from$tau) && length(init_from$tau) >= B) {
    tau0 <- as.numeric(init_from$tau)[seq_len(B)]
  }
  tau0[1] <- 0

  # Fast path with complete data: pooled starting values, then Newton steps
  # with the exact score and information. It is accepted when the first-order
  # (Karush-Kuhn-Tucker) conditions hold; otherwise the derivative-free
  # algorithm below continues from its result.
  if (!anyNA(y)) {
    fast <- .fit_constrained_newton(y, order, thinning, innovation,
                                    blocks_id = if (B > 1L) blocks_id else NULL, B = B,
                                    cinfo = cinfo, a1 = a1_0, a2 = a2_0, theta = theta0,
                                    kappa = kappa0, tau = tau0,
                                    nb_inno_size_ub = nb_inno_size_ub,
                                    ub_a1 = ub_a1, ub_a2 = ub_a2)
    if (!is.null(fast)) {
      if (isTRUE(fast$kkt_ok)) {
        f <- fast$fit
        f$convergence <- 0L
        f$settings <- list(order = order, thinning = thinning, innovation = innovation,
                           n_subjects = n, n_time = N, method = "constrained", path = "newton",
                           constraint = cinfo$description, constraint_df = cinfo$df,
                           constraint_info = cinfo, nb_inno_size_ub = nb_inno_size_ub)
        if (B == 1L) f$tau <- init_from$tau
        return(.finalize_inad_fit(f, y = y, na_action = na_action,
                                  blocks_id = if (B > 1L) blocks_id else NULL,
                                  block_levels = if (B > 1L) block_info$block_levels else NULL))
      }
      au <- .unpack_alpha(order, fast$fit$alpha, N)
      a1_0 <- au$a1; a2_0 <- au$a2
      theta0 <- fast$fit$theta
      if (is_nb) kappa0 <- fast$fit$nb_inno_size
      if (B > 1L) tau0 <- fast$fit$tau
    }
  }

  # Free-parameter layout: one group per parameter type.
  groups <- list(); par0 <- numeric(0); lower <- numeric(0); upper <- numeric(0)
  add_group <- function(name, init, lb, ub) {
    idx <- length(par0) + seq_along(init)
    groups[[name]] <<- idx
    par0 <<- c(par0, init)
    lower <<- c(lower, rep_len(lb, length(init)))
    upper <<- c(upper, rep_len(ub, length(init)))
  }
  clamp_a <- function(a, ub) pmin(pmax(a, pmin(1e-6, ub / 2)), ub * (1 - 1e-6))
  if (cinfo$alpha1_const) add_group("a1", clamp_a(mean(a1_0[2:N]), ub_a1), 0, ub_a1)
  else add_group("a1", clamp_a(a1_0[2:N], ub_a1), 0, ub_a1)
  if (order == 2) {
    if (cinfo$alpha2_const) add_group("a2", clamp_a(mean(a2_0[3:N]), ub_a2), 0, ub_a2)
    else add_group("a2", clamp_a(a2_0[3:N], ub_a2), 0, ub_a2)
  }
  if (cinfo$theta_const) add_group("theta", max(mean(theta0), 1e-6), eps, 1e6)
  else add_group("theta", pmax(theta0, 1e-6), eps, 1e6)
  if (is_nb) {
    if (cinfo$kappa_const) add_group("kappa", exp(mean(log(kappa0))), eps, nb_inno_size_ub)
    else add_group("kappa", kappa0, eps, nb_inno_size_ub)
  }
  if (B > 1L) add_group("tau", tau0[-1], -1e3, 1e3)

  unpack <- function(par) {
    a1 <- numeric(N)
    a1[2:N] <- par[groups$a1]  # a scalar recycles when alpha1 is constant
    alpha <- a1
    if (order == 2) {
      a2 <- numeric(N)
      a2[3:N] <- par[groups$a2]
      alpha <- cbind(a1, a2)
    }
    theta <- if (cinfo$theta_const) rep(par[groups$theta], N) else par[groups$theta]
    kappa <- NULL
    if (is_nb) kappa <- if (cinfo$kappa_const) rep(par[groups$kappa], N) else par[groups$kappa]
    tau <- if (B > 1L) c(0, par[groups$tau]) else 0
    list(alpha = alpha, theta = theta, kappa = kappa, tau = tau)
  }

  feasible <- function(theta, tau) {
    mu <- .inad_effective_innovation_mean(theta, 0, innovation)
    all(is.finite(mu)) && all(outer(tau, mu, "+") > 0)
  }

  objective <- function(par) {
    p <- unpack(par)
    if (!feasible(p$theta, p$tau)) return(1e20)
    ll <- tryCatch(
      logL_inad(y = y, order = order, thinning = thinning, innovation = innovation,
                alpha = p$alpha, theta = p$theta, nb_inno_size = p$kappa,
                blocks = if (B > 1L) blocks_id else NULL, tau = p$tau,
                na_action = na_action),
      error = function(e) -Inf
    )
    if (!is.finite(ll)) return(1e20)
    -ll
  }

  bobyqa <- function(x0, fn, lb, ub, xtol, maxeval) {
    opt <- tryCatch(
      nloptr::nloptr(x0 = x0, eval_f = fn, lb = lb, ub = ub,
                     opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = xtol, maxeval = maxeval)),
      error = function(e) NULL
    )
    if (is.null(opt)) return(list(par = x0, value = fn(x0)))
    list(par = opt$solution, value = opt$objective)
  }

  par <- pmin(pmax(par0, lower), upper)
  value <- objective(par)
  if (value >= 1e20) stop("Constrained fit: starting values are infeasible")
  converged <- FALSE
  for (round in seq_len(max_rounds)) {
    value_round <- value
    # Blockwise sweeps.
    for (sweep in seq_len(50)) {
      value_sweep <- value
      for (g in names(groups)) {
        idx <- groups[[g]]
        fn_g <- function(v) { p <- par; p[idx] <- v; objective(p) }
        if (length(idx) == 1L) {
          lo <- lower[idx]; hi <- upper[idx]
          if (g == "tau") {
            mu <- .inad_effective_innovation_mean(unpack(par)$theta, 0, innovation)
            lo <- max(lo, -min(mu) + 1e-8); hi <- min(hi, max(abs(par[idx]), 1) * 50)
          }
          if (g == "theta" || g == "kappa") hi <- min(hi, max(par[idx], 1) * 100)
          opt <- stats::optimize(fn_g, interval = c(lo, hi), tol = 1e-10)
          cand <- list(par = opt$minimum, value = opt$objective)
        } else {
          cand <- bobyqa(par[idx], fn_g, lower[idx], upper[idx], 1e-8, 3000)
        }
        if (cand$value < value) { par[idx] <- cand$par; value <- cand$value }
      }
      if (value_sweep - value < tol) break
    }
    # Joint refinement over all free parameters.
    cand <- bobyqa(par, objective, lower, upper, 1e-10, 20000)
    if (cand$value < value) { par <- cand$par; value <- cand$value }
    if (verbose) message(sprintf("constrained fit round %d: logL = %.8f", round, -value))
    if (value_round - value < tol) { converged <- TRUE; break }
  }

  est <- unpack(par)
  fit <- list(
    alpha = est$alpha, theta = est$theta,
    tau = if (B > 1L) est$tau else init_from$tau,
    nb_inno_size = est$kappa, log_l = -value,
    convergence = if (converged) 0L else 1L,
    settings = list(order = order, thinning = thinning, innovation = innovation,
                    n_subjects = n, n_time = N, method = "constrained", path = "search",
                    constraint = cinfo$description, constraint_df = cinfo$df,
                    constraint_info = cinfo, nb_inno_size_ub = nb_inno_size_ub)
  )
  # Finish with Newton steps using the exact (Louis) score and information.
  if (!anyNA(y)) {
    pol <- .inad_newton_polish(y, fit, blocks = if (B > 1L) blocks_id else NULL,
                               nb_inno_size_ub = nb_inno_size_ub)
    if (is.finite(pol$log_l) && pol$log_l >= fit$log_l) {
      fit$alpha <- pol$fit$alpha; fit$theta <- pol$fit$theta
      fit$nb_inno_size <- pol$fit$nb_inno_size
      if (B > 1L) fit$tau <- pol$fit$tau
      fit$log_l <- pol$log_l
    }
  }
  .finalize_inad_fit(fit, y = y, na_action = na_action,
                     blocks_id = if (B > 1L) blocks_id else NULL,
                     block_levels = if (B > 1L) block_info$block_levels else NULL)
}

#' Constrained fit by Newton steps from pooled starting values
#'
#' Replaces each constrained group of the unconstrained estimates by a pooled
#' value (thinning coefficients weighted by their predecessor totals, with the
#' innovation means shifted so that the fitted means are unchanged), then
#' takes Newton steps with the exact score and information.
#'
#' @return \code{NULL} if the starting values are infeasible, otherwise a list
#'   with the \code{fit} and \code{kkt_ok}.
#' @keywords internal
.fit_constrained_newton <- function(y, order, thinning, innovation, blocks_id, B, cinfo,
                                    a1, a2, theta, kappa, tau, nb_inno_size_ub,
                                    ub_a1, ub_a2) {
  N <- ncol(y)
  S <- colSums(y)
  ybar <- colMeans(y)
  mu <- .inad_effective_innovation_mean(theta, 0, innovation)
  pool <- function(a, lag, ub) {
    idx <- (lag + 1L):N
    ac <- sum(a[idx] * S[idx - lag]) / max(sum(S[idx - lag]), 1e-12)
    ac <- min(max(ac, 1e-6), ub * (1 - 1e-6))
    mu[idx] <<- mu[idx] + (a[idx] - ac) * ybar[idx - lag]
    a[idx] <- ac
    a
  }
  if (cinfo$alpha1_const) a1 <- pool(a1, 1L, ub_a1)
  if (order == 2 && cinfo$alpha2_const) a2 <- pool(a2, 2L, ub_a2)
  mu <- pmax(mu, 0.05 * pmax(ybar, 1e-3))
  if (cinfo$theta_const) mu[] <- mean(mu)
  theta <- if (innovation == "bell") .bell_theta_from_mean(mu) else mu
  if (!is.null(kappa) && cinfo$kappa_const) kappa[] <- exp(mean(log(kappa)))
  fit <- list(
    alpha = if (order == 2) cbind(a1, a2) else a1, theta = theta, nb_inno_size = kappa,
    tau = if (B > 1L) tau else NULL,
    settings = list(order = order, thinning = thinning, innovation = innovation,
                    n_subjects = nrow(y), n_time = N, constraint_info = cinfo,
                    nb_inno_size_ub = nb_inno_size_ub)
  )
  ll0 <- tryCatch(logL_inad(y = y, order = order, thinning = thinning, innovation = innovation,
                            alpha = fit$alpha, theta = theta, nb_inno_size = kappa,
                            blocks = blocks_id, tau = if (B > 1L) tau else 0),
                  error = function(e) -Inf)
  if (!is.finite(ll0)) return(NULL)
  fit$log_l <- ll0
  pol <- tryCatch(.inad_newton_polish(y, fit, blocks = blocks_id, nb_inno_size_ub = nb_inno_size_ub,
                                      max_iter = 100),
                  error = function(e) NULL)
  if (is.null(pol) || !is.finite(pol$log_l)) return(NULL)
  fit <- pol$fit
  fit$log_l <- pol$log_l
  k <- tryCatch(.inad_kkt_core(y, fit, blocks_id), error = function(e) NULL)
  list(fit = fit, kkt_ok = !is.null(k) && isTRUE(k$kkt_ok))
}

#' Run All Stationarity-Related Tests for INAD
#'
#' Runs the standard battery of time-invariance tests against a common
#' unconstrained fit. For order 1 the constraints are \code{"alpha"},
#' \code{"theta"}, and \code{"both"}; with negative binomial innovations they
#' are \code{"alpha"}, \code{"innovation"} (\code{theta} and \code{kappa}
#' jointly), and \code{"all"}. For order 2 they are \code{"alpha1"},
#' \code{"alpha2"}, \code{"alpha"}, \code{"theta"} (or \code{"innovation"}),
#' and \code{"all"}.
#'
#' @param y Integer matrix.
#' @param order Model order (1 or 2).
#' @param thinning Thinning operator.
#' @param innovation Innovation distribution.
#' @param blocks Optional block assignments.
#' @param verbose Logical.
#' @param ... Additional arguments.
#' @return A list with class "stationarity_tests_inad".
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
#' out <- run_stationarity_tests_inad(
#'   y,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   verbose = FALSE,
#'   max_iter = 8
#' )
#' out$summary
#' @export
run_stationarity_tests_inad <- function(y, order = 1, thinning = "binom", innovation = "pois",
                                        blocks = NULL, verbose = FALSE, ...) {
  if (!is.matrix(y)) y <- as.matrix(y)
  n_subjects <- nrow(y)
  fit_args <- list(...)
  if (anyNA(y) && is.null(fit_args$na_action)) {
    fit_args$na_action <- "marginalize"
  }
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  if (verbose) cat("Fitting unconstrained model...\n")
  fit_uncon <- do.call(
    fit_inad,
    c(
      list(y = y, order = order, thinning = thinning, innovation = innovation, blocks = blocks),
      fit_args
    )
  )
  is_nb <- identical(innovation, "nbinom")
  tests <- if (order == 1) {
    if (is_nb) c("alpha", "innovation", "all") else c("alpha", "theta", "both")
  } else {
    c("alpha1", "alpha2", "alpha", if (is_nb) "innovation" else "theta", "all")
  }
  results <- list()
  for (test in tests) {
    if (verbose) cat("Testing:", test, "constant...\n")
    results[[test]] <- do.call(
      test_stationarity_inad,
      c(list(y = y, order = order, thinning = thinning, innovation = innovation,
             blocks = blocks, constrain = test, fit_unconstrained = fit_uncon,
             verbose = FALSE), fit_args)
    )
  }
  summary_table <- data.frame(constraint = tests,
                              statistic = sapply(results, function(r) r$statistic),
                              df = sapply(results, function(r) r$df),
                              LRT = sapply(results, function(r) r$lrt_stat),
                              p_value = sapply(results, function(r) r$p_value),
                              BIC_unconstrained = sapply(results, function(r) r$bic_unconstrained),
                              BIC_constrained = sapply(results, function(r) r$bic_constrained),
                              BIC_selected = sapply(results, function(r) r$bic_selected),
                              stringsAsFactors = FALSE, row.names = NULL)
  output <- list(method = "lrt", fit_unconstrained = fit_uncon, tests = results, summary = summary_table,
                 settings = list(order = order, thinning = thinning, innovation = innovation,
                                 n_subjects = n_subjects, n_time = ncol(y)))
  class(output) <- "stationarity_tests_inad"
  output
}

#' @export
print.stationarity_tests_inad <- function(x, digits = 4, ...) {
  cat("\nStationarity Tests for INAD Model\n\n")
  tbl <- x$summary
  tbl$p_value <- format.pval(tbl$p_value, digits = 4)
  print(tbl, row.names = FALSE, digits = digits)
  invisible(x)
}
