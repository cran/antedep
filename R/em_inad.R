#' EM algorithm for INAD model estimation
#'
#' Fits INAD models using the Expectation-Maximization algorithm.
#' This is an alternative to direct likelihood optimization.
#'
#' For Gaussian and CAT EM entry points, see \code{\link{em_gau}} and
#' \code{\link{em_cat}}. For CAT specifically, \code{fit_cat()} supports
#' \code{na_action = "em"} for orders 0/1 and \code{na_action = "marginalize"}
#' for order 2 missing-data fits.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#' @param order Model order (1 or 2). Order 0 does not require EM.
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional integer vector of length n_subjects for block effects.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance for log-likelihood change.
#' @param alpha_init Optional initial values for alpha parameters.
#' @param theta_init Optional initial values for theta parameters.
#' @param tau_init Optional initial values for tau parameters.
#' @param nb_inno_size Size parameter for negative binomial innovation (if used).
#' @param safeguard Logical; if TRUE, use step-halving when likelihood decreases.
#' @param verbose Logical; if TRUE, print iteration progress.
#'
#' @return A list with class "inad_fit" containing estimated parameters.
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
                    blocks = NULL, max_iter = 200, tol = 1e-7,
                    alpha_init = NULL, theta_init = NULL, tau_init = NULL,
                    nb_inno_size = NULL, safeguard = TRUE, verbose = FALSE) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (any(y < 0) || any(y != floor(y))) stop("y must contain non-negative integers")
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (!order %in% c(1, 2)) stop("EM algorithm supports order 1 or 2.")
  if (order >= n_time) stop("order must be less than n_time")
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  block_info <- .normalize_blocks(blocks, n_subjects)
  blocks <- block_info$blocks_id
  n_blocks <- block_info$n_blocks
  
  if (!(thinning == "binom" && innovation == "pois")) {
    warning("EM is optimized for binom-pois. Other combinations use approximations.")
  }
  
  init <- .em_init(y, order, blocks, alpha_init, theta_init, tau_init, innovation = innovation)
  alpha <- init$alpha; theta <- init$theta; tau <- init$tau
  
  loglik_old <- .em_ll(y, order, thinning, innovation, alpha, theta, tau, blocks, nb_inno_size)
  loglik_trace <- loglik_old; converged <- FALSE; iter <- 0; loglik_new <- loglik_old
  
  if (verbose) cat(sprintf("EM iter 0: logL = %.6f\n", loglik_old))
  
  for (iter in 1:max_iter) {
    E_vals <- .em_e(y, order, thinning, innovation, alpha, theta, tau, blocks, nb_inno_size)
    updates <- .em_m(y, order, thinning, innovation, E_vals, blocks, nb_inno_size)
    alpha_new <- updates$alpha; theta_new <- updates$theta; tau_new <- updates$tau
    loglik_new <- .em_ll(y, order, thinning, innovation, alpha_new, theta_new, tau_new, blocks, nb_inno_size)
    
    if (safeguard && loglik_new + 1e-10 < loglik_old) {
      lam <- 0.5
      while (lam > 1/1024) {
        alpha_try <- (1 - lam) * alpha + lam * alpha_new
        theta_try <- (1 - lam) * theta + lam * theta_new
        tau_try <- (1 - lam) * tau + lam * tau_new
        loglik_try <- .em_ll(y, order, thinning, innovation, alpha_try, theta_try, tau_try, blocks, nb_inno_size)
        if (loglik_try + 1e-10 >= loglik_old) {
          alpha_new <- alpha_try; theta_new <- theta_try; tau_new <- tau_try; loglik_new <- loglik_try; break
        }
        lam <- lam / 2
      }
    }
    
    loglik_trace <- c(loglik_trace, loglik_new)
    if (verbose) cat(sprintf("EM iter %d: logL = %.6f\n", iter, loglik_new))
    
    if (abs(loglik_new - loglik_old) / (1 + abs(loglik_old)) < tol) {
      converged <- TRUE; alpha <- alpha_new; theta <- theta_new; tau <- tau_new; break
    }
    alpha <- alpha_new; theta <- theta_new; tau <- tau_new; loglik_old <- loglik_new
  }
  
  result <- list(
    alpha = alpha, theta = theta, tau = tau,
    log_l = loglik_new, iterations = iter, converged = converged, loglik_trace = loglik_trace,
    settings = list(order = order, thinning = thinning, innovation = innovation,
                    blocks = if (n_blocks > 1L) blocks else NULL,
                    block_levels = block_info$block_levels,
                    n_subjects = n_subjects, n_time = n_time, n_blocks = n_blocks, method = "em")
  )
  class(result) <- "inad_fit"
  result
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
  } else pmax(theta_init, 1e-6)
  if (is.null(alpha_init)) {
    if (order == 1) {
      alpha <- numeric(n_time); alpha[1] <- 0
      for (i in 2:n_time) alpha[i] <- pmin(pmax(mean(y[,i]) / (2 * max(mean(y[,i-1]), 1e-8)), 1e-6), 0.99)
    } else {
      alpha1 <- alpha2 <- numeric(n_time); alpha1[1:2] <- alpha2[1:2] <- 0
      for (i in 3:n_time) {
        alpha1[i] <- pmin(pmax(mean(y[,i]) / (3 * max(mean(y[,i-1]), 1e-8)), 0.1), 0.9)
        alpha2[i] <- pmin(pmax(mean(y[,i]) / (3 * max(mean(y[,i-2]), 1e-8)), 0.1), 0.9)
      }
      alpha <- cbind(alpha1, alpha2)
    }
  } else alpha <- alpha_init
  tau <- if (is.null(tau_init)) rep(0, n_blocks) else { tau_init[1] <- 0; tau_init }
  list(alpha = alpha, theta = theta, tau = tau)
}

#' @keywords internal
.em_e <- function(y, order, thinning, innovation, alpha, theta, tau, blocks, nb_inno_size) {
  n <- nrow(y); N <- ncol(y)
  alpha_list <- .unpack_alpha(order, alpha, N); a1 <- alpha_list$a1; a2 <- alpha_list$a2
  if (order == 1) {
    E_Z <- E_W <- matrix(0, n, N); E_Z[,1] <- 0; E_W[,1] <- y[,1]
    for (i in 2:N) for (s in 1:n) {
      lam <- .inad_effective_innovation_param(theta[i], tau[blocks[s]], innovation)
      y_i <- y[s,i]; y_i1 <- y[s,i-1]
      if (y_i1 == 0) { E_Z[s,i] <- 0; E_W[s,i] <- y_i } else {
        z_vals <- 0:min(y_i1, y_i); w_vals <- y_i - z_vals
        log_pz <- log(pmax(.thin_vec(z_vals, y_i1, a1[i], thinning), 1e-300))
        log_pw <- log(pmax(.innov_vec(w_vals, lam, innovation, nb_inno_size), 1e-300))
        post <- exp(log_pz + log_pw - .log_sum_exp(log_pz + log_pw))
        E_Z[s,i] <- sum(z_vals * post); E_W[s,i] <- y_i - E_Z[s,i]
      }
    }
    list(E_Z = E_Z, E_W = E_W)
  } else {
    E_Z1 <- E_Z2 <- E_W <- matrix(0, n, N); E_W[,1] <- y[,1]
    for (s in 1:n) {
      lam <- .inad_effective_innovation_param(theta[2], tau[blocks[s]], innovation)
      y_2 <- y[s,2]; y_1 <- y[s,1]
      if (y_1 > 0) {
        z_vals <- 0:min(y_1, y_2); w_vals <- y_2 - z_vals
        log_pz <- log(pmax(.thin_vec(z_vals, y_1, a1[2], thinning), 1e-300))
        log_pw <- log(pmax(.innov_vec(w_vals, lam, innovation, nb_inno_size), 1e-300))
        post <- exp(log_pz + log_pw - .log_sum_exp(log_pz + log_pw))
        E_Z1[s,2] <- sum(z_vals * post)
      }
      E_W[s,2] <- y[s,2] - E_Z1[s,2]
    }
    for (i in 3:N) for (s in 1:n) {
      lam <- .inad_effective_innovation_param(theta[i], tau[blocks[s]], innovation)
      y_i <- y[s,i]; y_i1 <- y[s,i-1]; y_i2 <- y[s,i-2]
      if (y_i1 == 0 && y_i2 == 0) { E_W[s,i] <- y_i; next }
      grid <- expand.grid(z1 = 0:min(y_i1, y_i), z2 = 0:min(y_i2, y_i))
      grid <- grid[grid$z1 + grid$z2 <= y_i, ]
      if (nrow(grid) == 0) { E_W[s,i] <- y_i; next }
      log_pz1 <- log(pmax(.thin_vec(grid$z1, y_i1, a1[i], thinning), 1e-300))
      log_pz2 <- log(pmax(.thin_vec(grid$z2, y_i2, a2[i], thinning), 1e-300))
      log_pw <- log(pmax(.innov_vec(y_i - grid$z1 - grid$z2, lam, innovation, nb_inno_size), 1e-300))
      post <- exp(log_pz1 + log_pz2 + log_pw - .log_sum_exp(log_pz1 + log_pz2 + log_pw))
      E_Z1[s,i] <- sum(grid$z1 * post); E_Z2[s,i] <- sum(grid$z2 * post)
      E_W[s,i] <- y_i - E_Z1[s,i] - E_Z2[s,i]
    }
    list(E_Z1 = E_Z1, E_Z2 = E_Z2, E_W = E_W)
  }
}

#' @keywords internal
.em_m <- function(y, order, thinning, innovation, E_vals, blocks, nb_inno_size) {
  n <- nrow(y); N <- ncol(y); B <- length(unique(blocks))
  if (order == 1) {
    E_Z <- E_vals$E_Z; E_W <- E_vals$E_W
    alpha <- numeric(N); alpha[1] <- 0
    for (i in 2:N) alpha[i] <- pmax(pmin(sum(E_Z[,i]) / max(sum(y[,i-1]), 1), 0.9999), 1e-12)
    theta <- numeric(N); tau <- rep(0, B)
    for (k in 1:10) {
      if (innovation == "bell") {
        for (i in 1:N) {
          mu_i <- pmax(mean(E_W[, i] - tau[blocks]), 1e-12)
          th_i <- .bell_theta_from_mean(mu_i)
          theta[i] <- if (is.finite(th_i)) pmax(th_i, 1e-12) else 1e-12
        }
        mu_theta <- .bell_mean_from_theta(theta)
        if (B > 1) {
          for (b in 2:B) {
            idx <- which(blocks == b)
            if (length(idx) > 0) {
              tau[b] <- mean(colMeans(E_W[idx, , drop = FALSE]) - mu_theta)
            }
          }
        }
      } else {
        for (i in 1:N) theta[i] <- pmax(mean(E_W[,i] - tau[blocks]), 1e-12)
        if (B > 1) for (b in 2:B) { idx <- which(blocks == b); if (length(idx) > 0) tau[b] <- mean(rowMeans(E_W[idx,,drop=FALSE]) - mean(theta)) }
      }
    }
    list(alpha = alpha, theta = theta, tau = tau)
  } else {
    E_Z1 <- E_vals$E_Z1; E_Z2 <- E_vals$E_Z2; E_W <- E_vals$E_W
    alpha1 <- alpha2 <- numeric(N); alpha1[1:2] <- alpha2[1:2] <- 0
    alpha1[2] <- pmax(pmin(sum(E_Z1[,2]) / max(sum(y[,1]), 1), 0.9999), 1e-12)
    for (i in 3:N) {
      alpha1[i] <- pmax(pmin(sum(E_Z1[,i]) / max(sum(y[,i-1]), 1), 0.9999), 1e-12)
      alpha2[i] <- pmax(pmin(sum(E_Z2[,i]) / max(sum(y[,i-2]), 1), 0.9999), 1e-12)
    }
    theta <- numeric(N); tau <- rep(0, B)
    for (k in 1:10) {
      if (innovation == "bell") {
        for (i in 1:N) {
          mu_i <- pmax(mean(E_W[, i] - tau[blocks]), 1e-12)
          th_i <- .bell_theta_from_mean(mu_i)
          theta[i] <- if (is.finite(th_i)) pmax(th_i, 1e-12) else 1e-12
        }
        mu_theta <- .bell_mean_from_theta(theta)
        if (B > 1) {
          for (b in 2:B) {
            idx <- which(blocks == b)
            if (length(idx) > 0) {
              tau[b] <- mean(colMeans(E_W[idx, , drop = FALSE]) - mu_theta)
            }
          }
        }
      } else {
        for (i in 1:N) theta[i] <- pmax(mean(E_W[,i] - tau[blocks]), 1e-12)
        if (B > 1) for (b in 2:B) { idx <- which(blocks == b); if (length(idx) > 0) tau[b] <- mean(rowMeans(E_W[idx,,drop=FALSE]) - mean(theta)) }
      }
    }
    list(alpha = cbind(alpha1, alpha2), theta = theta, tau = tau)
  }
}

#' @keywords internal
.em_ll <- function(y, order, thinning, innovation, alpha, theta, tau, blocks, nb_inno_size) {
  tryCatch(logL_inad(y = y, order = order, thinning = thinning, innovation = innovation,
                     alpha = alpha, theta = theta, blocks = blocks, tau = tau, nb_inno_size = nb_inno_size),
           error = function(e) -Inf)
}

#' Print method for INAD model fits
#'
#' @param x An object of class \code{inad_fit}.
#' @param digits Number of digits to print.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @export
print.inad_fit <- function(x, digits = 4, ...) {
  cat("\nINAD Model Fit")
  if (!is.null(x$settings$method) && x$settings$method == "em") cat(" (EM Algorithm)")
  cat("\n", paste(rep("=", 40), collapse = ""), "\n\n")
  cat("Order:", x$settings$order, " Thinning:", x$settings$thinning, " Innovation:", x$settings$innovation, "\n")
  cat("Log-likelihood:", round(x$log_l, digits), "\n")
  if (!is.null(x$aic)) cat("AIC:", round(x$aic, digits), "\n")
  if (!is.null(x$bic)) cat("BIC:", round(x$bic, digits), "\n")
  if (!is.null(x$n_params)) cat("Parameters:", x$n_params, "\n")
  if (!is.null(x$n_missing)) cat("Missing values:", x$n_missing, "\n")
  if (!is.null(x$convergence)) cat("Convergence code:", x$convergence, "\n")
  if (!is.null(x$iterations)) cat("Iterations:", x$iterations, " Converged:", x$converged, "\n")
  invisible(x)
}
