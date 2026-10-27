#' EM Algorithm for INAD Model Estimation
#'
#' Fits INAD models using the Expectation-Maximization algorithm, treating the
#' thinned counts as missing data. This is an alternative to direct likelihood
#' optimization by \code{\link{fit_inad}}.
#'
#' For Gaussian and CAT EM entry points, see \code{\link{em_gau}} and
#' \code{\link{em_cat}}. For CAT specifically, \code{fit_cat()} supports
#' \code{na_action = "em"} for orders 0/1 and \code{na_action = "marginalize"}
#' for order 2 missing-data fits.
#'
#' @details
#' The complete data for subject \eqn{s} at time \eqn{i} are the thinned
#' counts \eqn{Z_{si,j} = \alpha_{i,i-j} \odot Y_{s,i-j}} and the innovation
#' \eqn{w_{si} = y_{si} - \sum_j Z_{si,j}}. Given the observed data, the
#' latent vector at time \eqn{i} has finite support and a conditional
#' distribution that depends only on \eqn{y_{si}} and its \eqn{p}
#' predecessors. The E-step computes this conditional distribution exactly, so
#' the expected complete-data log-likelihood is evaluated without
#' approximation for every thinning-innovation combination.
#'
#' For binomial, Poisson, and negative binomial thinning, the
#' parameter-dependent part of the thinning log-likelihood is linear in the
#' thinned count, and the M-step for the thinning coefficient has the closed
#' form \eqn{\hat\alpha_{i,i-j} = \sum_s E(Z_{si,j} \mid y_s) / \sum_s
#' y_{s,i-j}}. For Poisson and Bell innovations the innovation log-likelihood is
#' also linear in \eqn{w_{si}}. For negative binomial innovations it contains
#' \eqn{\log\Gamma(w_{si} + \kappa_i)}, which is not linear in \eqn{w_{si}};
#' its conditional expectation is computed from the conditional distribution of
#' \eqn{w_{si}}, and the size parameters \eqn{\kappa_i} are updated by
#' maximizing the resulting expected log-likelihood. Without block effects the
#' innovation updates are exact maximizations. With block effects the
#' innovation means, sizes, and block effects are updated by conditional
#' maximization (an ECM algorithm), which also never decreases the
#' observed-data log-likelihood.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param order Model order (1 or 2). Order 0 does not require EM.
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional vector of length n_subjects for block effects.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance for the absolute change in the
#'   observed-data log-likelihood between iterations.
#' @param alpha_init Optional initial values for alpha parameters.
#' @param theta_init Optional initial values for theta parameters.
#' @param tau_init Optional initial values for tau parameters.
#' @param nb_inno_size Optional starting value(s), scalar or length n_time,
#'   for the negative binomial innovation size parameters, which are estimated.
#'   Default 1.
#' @param nb_inno_size_ub Upper bound for the negative binomial innovation size
#'   parameters. Default 50, as in \code{\link{fit_inad}}.
#' @param safeguard Logical; if TRUE, warn when the log-likelihood decreases by
#'   more than numerical tolerance.
#' @param newton_polish Logical; if TRUE, finish with Newton-Raphson steps that
#'   use the exact score and observed information from Louis' identity (see
#'   \code{\link{ci_inad}}). EM converges linearly, so this removes the small
#'   remaining distance to the maximum. Default FALSE (pure EM).
#' @param verbose Logical; if TRUE, print iteration progress.
#'
#' @return A list with class "inad_fit" containing estimated parameters,
#'   the maximized log-likelihood, \code{n_params}, \code{aic}, \code{bic},
#'   the number of iterations, and the log-likelihood trace.
#'
#' @references
#' Dempster, A.P., Laird, N.M. and Rubin, D.B. (1977). Maximum likelihood from
#' incomplete data via the EM algorithm. \emph{Journal of the Royal Statistical
#' Society, Series B}, 39, 1-38.
#'
#' Meng, X.-L. and Rubin, D.B. (1993). Maximum likelihood estimation via the
#' ECM algorithm: a general framework. \emph{Biometrika}, 80, 267-278.
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(
#'   n_subjects = 50,
#'   n_time = 5,
#'   order = 1,
#'   thinning = "binom",
#'   innovation = "pois",
#'   alpha = 0.25,
#'   theta = 2
#' )
#' fit <- em_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 20, tol = 1e-6)
#' fit$log_l
#' @seealso \code{\link{em_gau}}, \code{\link{em_cat}}, \code{\link{fit_inad}},
#'   \code{\link{fit_cat}}
#' @export
em_inad <- function(y, order = 1, thinning = "binom", innovation = "pois",
                    blocks = NULL, max_iter = 2000, tol = 1e-9,
                    alpha_init = NULL, theta_init = NULL, tau_init = NULL,
                    nb_inno_size = NULL, nb_inno_size_ub = 50,
                    safeguard = TRUE, newton_polish = FALSE, verbose = FALSE) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (anyNA(y)) stop("em_inad requires complete data; use fit_inad(na_action = ...) for missing values")
  if (any(y < 0) || any(y != floor(y))) stop("y must contain non-negative integers")
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (!order %in% c(1, 2)) stop("EM algorithm supports order 1 or 2.")
  if (order >= n_time) stop("order must be less than n_time")
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  block_info <- .normalize_blocks(blocks, n_subjects)
  blocks_id <- block_info$blocks_id
  n_blocks <- block_info$n_blocks

  init <- .em_init(y, order, blocks_id, alpha_init, theta_init, tau_init, innovation = innovation)
  au <- .unpack_alpha(order, init$alpha, n_time)
  par <- list(a1 = au$a1, a2 = au$a2, theta = init$theta, tau = init$tau,
              kappa = if (innovation == "nbinom") {
                k0 <- if (is.null(nb_inno_size)) 1 else as.numeric(nb_inno_size)
                pmin(pmax(rep_len(k0, n_time), 1e-6), nb_inno_size_ub)
              } else NULL)

  es <- .em_inad_estep(y, order, thinning, innovation, par, blocks_id, n_blocks)
  loglik_old <- es$loglik
  if (!is.finite(loglik_old)) stop("Initial values give a non-finite log-likelihood")
  loglik_trace <- loglik_old; converged <- FALSE; iter <- 0

  if (verbose) cat(sprintf("EM iter 0: logL = %.8f\n", loglik_old))

  for (iter in seq_len(max_iter)) {
    par_new <- .em_inad_mstep(y, order, thinning, innovation, es, par, blocks_id, n_blocks,
                              nb_inno_size_ub)
    es_new <- .em_inad_estep(y, order, thinning, innovation, par_new, blocks_id, n_blocks)
    loglik_new <- es_new$loglik
    if (!is.finite(loglik_new)) {
      warning("EM produced a non-finite log-likelihood; returning the last finite iterate")
      break
    }
    if (safeguard && loglik_new < loglik_old - 1e-6 * (1 + abs(loglik_old))) {
      warning(sprintf("EM log-likelihood decreased at iteration %d (%.8f -> %.8f)",
                      iter, loglik_old, loglik_new))
    }
    par <- par_new; es <- es_new
    loglik_trace <- c(loglik_trace, loglik_new)
    if (verbose) cat(sprintf("EM iter %d: logL = %.8f\n", iter, loglik_new))
    if (abs(loglik_new - loglik_old) < tol) { converged <- TRUE; loglik_old <- loglik_new; break }
    loglik_old <- loglik_new
  }

  alpha <- if (order == 1) par$a1 else cbind(alpha1 = par$a1, alpha2 = par$a2)
  result <- list(
    alpha = alpha, theta = par$theta,
    tau = if (n_blocks > 1L) par$tau else NULL,
    nb_inno_size = par$kappa,
    log_l = loglik_old, iterations = iter, converged = converged,
    convergence = if (converged) 0L else 1L,
    loglik_trace = loglik_trace,
    settings = list(order = order, thinning = thinning, innovation = innovation,
                    n_subjects = n_subjects, n_time = n_time, method = "em",
                    nb_inno_size_ub = nb_inno_size_ub)
  )
  if (newton_polish) {
    pol <- .inad_newton_polish(y, result, blocks = if (n_blocks > 1L) blocks_id else NULL,
                               nb_inno_size_ub = nb_inno_size_ub)
    if (is.finite(pol$log_l) && pol$log_l >= result$log_l) {
      result$alpha <- pol$fit$alpha; result$theta <- pol$fit$theta
      result$nb_inno_size <- pol$fit$nb_inno_size
      if (n_blocks > 1L) result$tau <- pol$fit$tau
      result$log_l_em <- result$log_l
      result$log_l <- pol$log_l
      result$settings$method <- "em+newton"
    }
  }
  .finalize_inad_fit(result, y = y, na_action = "fail",
                     blocks_id = if (n_blocks > 1L) blocks_id else NULL,
                     block_levels = if (n_blocks > 1L) block_info$block_levels else NULL)
}

#' @keywords internal
.em_init <- function(y, order, blocks, alpha_init, theta_init, tau_init, innovation = "pois") {
  n_time <- ncol(y); n_blocks <- length(unique(blocks))
  theta <- if (is.null(theta_init)) {
    mu0 <- pmax(colMeans(y), 1e-6)
    if (innovation == "bell") {
      th0 <- .bell_theta_from_mean(mu0)
      th0[!is.finite(th0)] <- 1
      pmax(th0, 1e-6)
    } else {
      mu0
    }
  } else pmax(rep_len(as.numeric(theta_init), n_time), 1e-6)
  if (is.null(alpha_init)) {
    if (order == 1) {
      alpha <- numeric(n_time); alpha[1] <- 0
      for (i in 2:n_time) alpha[i] <- pmin(pmax(mean(y[,i]) / (2 * max(mean(y[,i-1]), 1e-8)), 1e-6), 0.99)
    } else {
      alpha1 <- alpha2 <- numeric(n_time)
      alpha1[2] <- pmin(pmax(mean(y[,2]) / (2 * max(mean(y[,1]), 1e-8)), 0.1), 0.9)
      for (i in 3:n_time) {
        alpha1[i] <- pmin(pmax(mean(y[,i]) / (3 * max(mean(y[,i-1]), 1e-8)), 0.1), 0.9)
        alpha2[i] <- pmin(pmax(mean(y[,i]) / (3 * max(mean(y[,i-2]), 1e-8)), 0.1), 0.9)
      }
      alpha <- cbind(alpha1, alpha2)
    }
  } else alpha <- alpha_init
  tau <- if (is.null(tau_init)) rep(0, n_blocks) else {
    tau_init <- rep_len(as.numeric(tau_init), n_blocks); tau_init[1] <- 0; tau_init
  }
  list(alpha = alpha, theta = theta, tau = tau)
}

#' Innovation pmf parameter by block and time
#'
#' Returns a \code{n_blocks} by \code{n_time} matrix of the parameter used by
#' \code{.innov_vec()}: the mean for Poisson and negative binomial innovations,
#' and the Bell parameter for Bell innovations.
#'
#' @keywords internal
.em_inad_lam <- function(theta, tau, n_blocks, innovation) {
  out <- matrix(NA_real_, n_blocks, length(theta))
  for (b in seq_len(n_blocks)) {
    out[b, ] <- .inad_effective_innovation_param(theta, tau[b], innovation)
  }
  out
}

#' EM E-step for INAD models
#'
#' For every subject and time, computes the exact conditional distribution of
#' the latent thinned counts given the observed data, and accumulates (i) the
#' conditional means of the thinned counts, (ii) the conditional pmf of the
#' innovation, aggregated by block and time, and (iii) the observed-data
#' log-likelihood, which is the sum of the log normalizing constants.
#'
#' @keywords internal
.em_inad_estep <- function(y, order, thinning, innovation, par, blocks, n_blocks) {
  n <- nrow(y); N <- ncol(y); ymax <- max(y)
  lam <- .em_inad_lam(par$theta, par$tau, n_blocks, innovation)
  if (any(!is.finite(lam)) || any(lam <= 0)) return(list(loglik = -Inf))
  kappa <- par$kappa
  EZ1 <- EZ2 <- matrix(0, n, N)
  Wtab <- array(0, dim = c(n_blocks, ymax + 1L, N))
  loglik <- 0
  for (i in seq_len(N)) {
    sz <- if (innovation == "nbinom") kappa[i] else NA_real_
    for (s in seq_len(n)) {
      b <- blocks[s]; yi <- y[s, i]; vals <- 0:yi
      lfw <- log(.innov_vec(vals, lam[b, i], innovation, sz))  # index w + 1
      if (i == 1) {
        loglik <- loglik + lfw[yi + 1L]
        Wtab[b, yi + 1L, 1L] <- Wtab[b, yi + 1L, 1L] + 1
        next
      }
      if (order == 1 || i == 2) {
        lp <- log(.thin_vec(vals, y[s, i - 1], par$a1[i], thinning)) + rev(lfw)
        lZ <- .log_sum_exp(lp)
        if (!is.finite(lZ)) return(list(loglik = -Inf))
        post <- exp(lp - lZ)
        loglik <- loglik + lZ
        EZ1[s, i] <- sum(vals * post)
        idx <- yi - vals + 1L
        Wtab[b, idx, i] <- Wtab[b, idx, i] + post
      } else {
        lg1 <- log(.thin_vec(vals, y[s, i - 1], par$a1[i], thinning))
        lg2 <- log(.thin_vec(vals, y[s, i - 2], par$a2[i], thinning))
        K <- outer(vals, vals, "+")
        ok <- K <= yi
        lp <- outer(lg1, lg2, "+")
        lp[!ok] <- -Inf
        lp[ok] <- lp[ok] + lfw[yi - K[ok] + 1L]
        lZ <- .log_sum_exp(lp)
        if (!is.finite(lZ)) return(list(loglik = -Inf))
        post <- exp(lp - lZ)
        loglik <- loglik + lZ
        EZ1[s, i] <- sum(rowSums(post) * vals)
        EZ2[s, i] <- sum(colSums(post) * vals)
        pk <- as.numeric(rowsum(post[ok], K[ok], reorder = TRUE))  # P(Z1 + Z2 = k), k = 0..yi
        idx <- yi - vals + 1L
        Wtab[b, idx, i] <- Wtab[b, idx, i] + pk
      }
    }
  }
  list(loglik = loglik, EZ1 = EZ1, EZ2 = EZ2, Wtab = Wtab)
}

#' EM M-step for INAD models
#'
#' Thinning coefficients use the closed-form update shared by all three
#' operators. Innovation parameters maximize the expected complete-data
#' innovation log-likelihood, evaluated from the conditional innovation pmf.
#'
#' @keywords internal
.em_inad_mstep <- function(y, order, thinning, innovation, es, par, blocks, n_blocks,
                           nb_inno_size_ub = 50, n_cm = 3) {
  N <- ncol(y)
  a_lo <- 1e-12
  a_hi <- if (thinning == "binom") 1 - 1e-12 else Inf
  a1 <- par$a1; a2 <- par$a2
  # When all predecessor counts are zero the coefficient is not identifiable
  # (the update is 0/0) and is set to zero.
  for (i in 2:N) {
    den <- sum(y[, i - 1])
    a1[i] <- if (den > 0) min(max(sum(es$EZ1[, i]) / den, a_lo), a_hi) else 0
  }
  if (order == 2) {
    for (i in 3:N) {
      den <- sum(y[, i - 2])
      a2[i] <- if (den > 0) min(max(sum(es$EZ2[, i]) / den, a_lo), a_hi) else 0
    }
  }

  W <- es$Wtab; wv <- 0:(dim(W)[2] - 1L)
  theta <- par$theta; tau <- par$tau; kappa <- par$kappa
  is_nb <- innovation == "nbinom"
  n_b <- tabulate(blocks, n_blocks)
  ymax <- max(y)

  # Expected complete-data innovation log-likelihood at time i for block b.
  q_ib <- function(i, b, th, kp, tb) {
    wts <- W[b, , i]
    nz <- wts > 0
    if (!any(nz)) return(0)
    lam <- .inad_effective_innovation_param(th, tb, innovation)
    if (!is.finite(lam) || lam <= 0) return(-Inf)
    p <- .innov_vec(wv[nz], lam, innovation, kp)
    if (any(!(p > 0))) return(-Inf)
    sum(wts[nz] * log(p))
  }
  q_i <- function(i, th, kp, tau_vec) {
    sum(vapply(seq_len(n_blocks), function(b) q_ib(i, b, th, kp, tau_vec[b]), numeric(1)))
  }
  maximize_1d <- function(f, lo, hi) {
    opt <- stats::optimize(function(x) -f(x), interval = c(lo, hi), tol = 1e-12)
    list(par = opt$minimum, value = -opt$objective)
  }
  # Lower bound on theta (on its own scale) so every block's innovation mean is positive.
  theta_lo <- function(tau_vec) {
    m_lo <- max(1e-10, -min(tau_vec) + 1e-10)
    if (innovation == "bell") .bell_theta_from_mean(m_lo) else m_lo
  }
  theta_hi <- if (innovation == "bell") .bell_theta_from_mean(10 * ymax + 10) else 10 * ymax + 10
  lk_lo <- log(1e-6); lk_hi <- log(nb_inno_size_ub)

  if (n_blocks == 1L) {
    for (i in seq_len(N)) {
      mbar <- max(sum(wv * W[1, , i]) / n_b[1], 1e-10)
      theta[i] <- if (innovation == "bell") .bell_theta_from_mean(mbar) else mbar
      if (is_nb) {
        # The mean update is exact for any size, so the joint maximizer is
        # (mbar, argmax over size).
        kappa[i] <- exp(maximize_1d(function(lk) q_i(i, theta[i], exp(lk), tau), lk_lo, lk_hi)$par)
      }
    }
  } else {
    for (cm in seq_len(n_cm)) {
      for (i in seq_len(N)) {
        lo <- theta_lo(tau)
        cur <- q_i(i, theta[i], if (is_nb) kappa[i] else NA_real_, tau)
        for (inner in seq_len(if (is_nb) 3L else 1L)) {
          cand <- maximize_1d(function(th) q_i(i, th, if (is_nb) kappa[i] else NA_real_, tau),
                              lo, max(theta_hi, 2 * theta[i]))
          if (cand$value >= cur) { theta[i] <- cand$par; cur <- cand$value }
          if (is_nb) {
            cand <- maximize_1d(function(lk) q_i(i, theta[i], exp(lk), tau), lk_lo, lk_hi)
            if (cand$value >= cur) { kappa[i] <- exp(cand$par); cur <- cand$value }
          }
        }
      }
      mu_theta <- .inad_effective_innovation_mean(theta, 0, innovation)
      for (b in 2:n_blocks) {
        q_b <- function(tb) sum(vapply(seq_len(N), function(i) {
          q_ib(i, b, theta[i], if (is_nb) kappa[i] else NA_real_, tb)
        }, numeric(1)))
        cur <- q_b(tau[b])
        cand <- maximize_1d(q_b, -min(mu_theta) + 1e-10, max(10 * ymax + 10, 2 * abs(tau[b])))
        if (cand$value >= cur) tau[b] <- cand$par
      }
    }
  }
  list(a1 = a1, a2 = a2, theta = theta, tau = tau, kappa = kappa)
}

#' Summary Method for INAD Model Fits
#'
#' Summary method for inad_fit objects
#'
#' @param object An \code{inad_fit} object.
#' @param ... Unused.
#' @return \code{object}, invisibly.
#' @export
summary.inad_fit <- function(object, ...) {
    cat("INAD Antedependence Model -- Summary\n")
    cat("====================================\n\n")
    cat(sprintf("  Order            : %d\n",  object$settings$order))
    cat(sprintf("  Thinning         : %s\n",  object$settings$thinning))
    cat(sprintf("  Innovation       : %s\n",  object$settings$innovation))
    cat(sprintf("  Subjects (n)     : %d\n",  object$settings$n_subjects))
    cat(sprintf("  Time points (T)  : %d\n",  object$settings$n_time))
    if (!is.null(object$settings$constraint))
        cat(sprintf("  Constraint       : %s\n",  object$settings$constraint))
    if (!is.null(object$n_params))
        cat(sprintf("  Parameters       : %d\n",  object$n_params))
    if (!is.null(object$n_missing))
        cat(sprintf("  Missing values   : %d (%.1f%%)\n",
                    object$n_missing, object$pct_missing))
    cat("\nFit statistics:\n")
    cat(sprintf("  Log-likelihood   : %.4f\n", object$log_l))
    if (!is.null(object$aic)) cat(sprintf("  AIC              : %.4f\n", object$aic))
    if (!is.null(object$bic)) cat(sprintf("  BIC              : %.4f\n", object$bic))
    if (!is.null(object$theta)) {
        cat("\nEstimated innovation mean (theta):\n")
        print(round(as.numeric(object$theta), 4))
    }
    if (!is.null(object$alpha)) {
        cat("\nEstimated thinning parameter (alpha):\n")
        print(round(as.numeric(object$alpha), 4))
    }
    invisible(object)
}

#' Print Method for INAD Model Fits
#'
#' @param x An object of class \code{inad_fit}.
#' @param digits Number of digits to print.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.inad_fit <- function(x, digits = 4, ...) {
  cat("\nINAD Model Fit")
  if (!is.null(x$settings$method) && startsWith(x$settings$method, "em")) cat(" (EM Algorithm)")
  cat("\n", paste(rep("=", 40), collapse = ""), "\n\n")
  cat("Order:", x$settings$order, " Thinning:", x$settings$thinning, " Innovation:", x$settings$innovation, "\n")
  if (!is.null(x$settings$constraint)) cat("Constraint:", x$settings$constraint, "\n")
  cat("Log-likelihood:", round(x$log_l, digits), "\n")
  if (!is.null(x$aic)) cat("AIC:", round(x$aic, digits), "\n")
  if (!is.null(x$bic)) cat("BIC:", round(x$bic, digits), "\n")
  if (!is.null(x$n_params)) cat("Parameters:", x$n_params, "\n")
  if (!is.null(x$n_missing)) cat("Missing values:", x$n_missing, "\n")
  if (!is.null(x$convergence)) cat("Convergence code:", x$convergence, "\n")
  if (!is.null(x$iterations)) cat("Iterations:", x$iterations, " Converged:", x$converged, "\n")
  invisible(x)
}
