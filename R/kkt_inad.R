#' First-Order Optimality Check for a Fitted INAD Model
#'
#' Evaluates the observed-data score of a fitted INAD model at its estimates,
#' computed exactly by Louis' identity (Fisher's identity for the score), and
#' checks the first-order (Karush-Kuhn-Tucker) conditions for a maximum over
#' the parameter space.
#'
#' The conditions are a zero score for every interior estimate, a nonpositive
#' one-sided score for a thinning coefficient or an innovation parameter at
#' zero (an innovation mean on the boundary), a nonnegative one-sided score
#' for a binomial thinning coefficient at one, and a nonnegative score for a
#' negative binomial innovation size at its numerical cap
#' (\code{nb_inno_size_ub}). At \eqn{\alpha = 0}, and at binomial
#' \eqn{\alpha = 1}, the thinned count is degenerate and the score is the exact
#' one-sided derivative of the observed-data log-likelihood.
#'
#' Interior stationarity is checked with the observed information \eqn{I} of
#' the interior parameters, which must be positive definite. The Newton
#' decrement \eqn{s^\top I^{-1} s}, where \eqn{s} is the score, is about twice
#' the log-likelihood increase that a further Newton step would give, and must
#' not exceed \code{tol}; each interior parameter must also satisfy
#' \eqn{s_j^2 / I_{jj} \le} \code{tol}, about twice the increase available by
#' moving it alone (this is implied by the decrement condition). A boundary
#' estimate whose score points into the parameter space passes only if the
#' largest log-likelihood increase available by moving it alone into the
#' interior, found by a line search on the observed-data log-likelihood, is at
#' most \code{tol / 2}. Thinning coefficients whose predecessor counts are all
#' zero are not identifiable; they are reported with status
#' \code{"not identifiable"} and are not checked.
#'
#' Satisfying these conditions shows that the fit is a stationary point that
#' meets the boundary conditions, to the stated tolerance. It does not show
#' that the fit is the global maximum.
#'
#' @param fit An \code{inad_fit} object from \code{\link{fit_inad}} or
#'   \code{\link{em_inad}}, fitted to complete data.
#' @param y Count matrix. Defaults to the data stored in \code{fit}.
#' @param blocks Optional block labels. Defaults to those stored in \code{fit}.
#' @param tol Tolerance, on the scale of twice a log-likelihood increase, for
#'   the Newton decrement, the single-parameter checks, and the boundary
#'   line searches.
#'
#' @return An object of class \code{"kkt_inad"}: a list with \code{table} (one
#'   row per free parameter, with columns \code{param}, \code{est},
#'   \code{score}, \code{status}, \code{check}, the quantity compared with
#'   \code{tol}, and \code{ok}), \code{info_pd}, \code{min_info_eigen} (the
#'   smallest eigenvalue of the interior observed information),
#'   \code{newton_decrement}, \code{max_abs_interior_score}, \code{n_lower},
#'   \code{n_upper}, \code{n_size_cap}, and \code{kkt_ok}.
#'
#' @seealso \code{\link{fit_inad}}, \code{\link{em_inad}}, \code{\link{ci_inad}}
#'
#' @examples
#' set.seed(1)
#' y <- simulate_inad(n_subjects = 60, n_time = 5, order = 1, thinning = "binom",
#'                    innovation = "pois", alpha = 0.3, theta = 2)
#' fit <- em_inad(y, order = 1, thinning = "binom", innovation = "pois",
#'                newton_polish = TRUE)
#' kkt_inad(fit)
#' @export
kkt_inad <- function(fit, y = NULL, blocks = NULL, tol = 1e-6) {
    if (!inherits(fit, "inad_fit")) stop("fit must be an 'inad_fit' object")
    if (is.null(y)) y <- fit$.y
    if (is.null(y)) stop("y must be supplied when it is not stored in fit")
    if (!is.matrix(y)) y <- as.matrix(y)
    if (anyNA(y)) stop("kkt_inad requires complete data")
    if (is.null(blocks)) blocks <- fit$.blocks
    if (fit$settings$order < 1) stop("kkt_inad requires order 1 or 2")
    .inad_kkt_core(y, fit, blocks, tol)
}

#' First-order conditions for a fit-like list (see kkt_inad)
#' @keywords internal
.inad_kkt_core <- function(y, fit, blocks, tol = 1e-6) {
    N <- ncol(y)
    thinning <- fit$settings$thinning
    kappa_ub <- if (!is.null(fit$settings$nb_inno_size_ub)) fit$settings$nb_inno_size_ub else 50

    li <- .louis_full_info(y, fit, blocks)
    tab <- li$fp$tab
    score <- as.numeric(li$score)
    est <- tab$est
    is_alpha <- tab$group %in% c("a1", "a2")

    unid <- rep(FALSE, nrow(tab))
    for (lag in 1:2) {
        g <- paste0("a", lag)
        k <- which(tab$group == g & !is.na(tab$index))
        if (length(k)) unid[k] <- .inad_alpha_unidentified(y, lag)[tab$index[k]]
    }
    lower <- ((is_alpha & !unid) | tab$group == "theta") & est <= 1e-8
    upper <- is_alpha & thinning == "binom" & est >= 1 - 1e-8
    cap <- tab$group == "kappa" & est >= kappa_ub * (1 - 1e-6)
    interior <- !(lower | upper | cap | unid)

    status <- rep("interior", nrow(tab))
    status[lower] <- "lower bound"
    status[upper] <- "upper bound"
    status[cap] <- "size cap"
    status[unid] <- "not identifiable"

    ok <- rep(TRUE, nrow(tab))
    check <- rep(0, nrow(tab))

    # Interior parameters: positive definite information, Newton decrement,
    # and the single-parameter quantities s_j^2 / I_jj.
    I_in <- li$info[interior, interior, drop = FALSE]
    I_in <- (I_in + t(I_in)) / 2
    s_in <- score[interior]
    info_pd <- TRUE
    min_eig <- NA_real_
    dec <- 0
    if (any(interior)) {
        min_eig <- min(eigen(I_in, symmetric = TRUE, only.values = TRUE)$values)
        R <- if (is.finite(min_eig) && min_eig > 0) tryCatch(chol(I_in), error = function(e) NULL) else NULL
        info_pd <- !is.null(R)
        dec <- if (info_pd) sum(backsolve(R, s_in, transpose = TRUE)^2) else NA_real_
        d_in <- diag(I_in)
        check[interior] <- ifelse(d_in > 0, s_in^2 / d_in, Inf)
        ok[interior] <- check[interior] <= tol
    }

    # Boundary parameters whose score points into the parameter space: twice
    # the log-likelihood increase from a line search along that parameter.
    wrong <- (lower & score > 0) | ((upper | cap) & score < 0)
    if (any(wrong)) {
        loglik <- .inad_loglik_fun(y, fit, blocks)
        ll0 <- loglik(fit)
        ymax <- max(y)
        for (k in which(wrong)) {
            lim <- if (upper[k]) 0
                   else if (cap[k]) est[k] * 1e-3
                   else if (is_alpha[k]) (if (thinning == "binom") 1 - 1e-10 else est[k] + 1e3)
                   else est[k] + 10 * ymax + 10
            ikk <- li$info[k, k]
            t0 <- if (is.finite(ikk) && ikk > 0) abs(score[k]) / ikk else 1e-4 * max(1, abs(est[k]))
            cs <- .inad_coord_search(loglik, fit, li$fp, k, lim, t0, N, ll0 = ll0)
            check[k] <- 2 * cs$gain
            ok[k] <- check[k] <= tol
        }
    }

    out <- list(
        table = data.frame(param = tab$param, est = est, score = score, status = status,
                           check = check, ok = ok, stringsAsFactors = FALSE),
        info_pd = info_pd,
        min_info_eigen = min_eig,
        newton_decrement = dec,
        max_abs_interior_score = if (any(interior)) max(abs(s_in)) else 0,
        n_lower = sum(lower), n_upper = sum(upper), n_size_cap = sum(cap),
        n_unidentified = sum(unid),
        kkt_ok = info_pd && is.finite(dec) && dec <= tol && all(ok),
        tol = tol
    )
    class(out) <- "kkt_inad"
    out
}

#' @export
print.kkt_inad <- function(x, digits = 4, ...) {
    cat("\nFirst-order optimality check for an INAD fit\n\n")
    if (isTRUE(x$info_pd)) {
        cat("Interior parameters: Newton decrement", format(x$newton_decrement, digits = digits),
            ", max |score|", format(x$max_abs_interior_score, digits = digits), "\n")
    } else {
        cat("Interior parameters: observed information NOT positive definite",
            "(smallest eigenvalue", format(x$min_info_eigen, digits = digits), ")\n")
    }
    cat("At lower bound:", x$n_lower, " at binomial upper bound:", x$n_upper,
        " at size cap:", x$n_size_cap, " not identifiable:", x$n_unidentified, "\n")
    bad <- x$table[!x$table$ok, , drop = FALSE]
    if (nrow(bad)) {
        cat("\nParameters failing the check (column 'check' is compared with the tolerance):\n")
        print(bad, row.names = FALSE, digits = digits)
    }
    cat("\nFirst-order conditions", if (isTRUE(x$kkt_ok)) "satisfied" else "NOT satisfied",
        "at tolerance", format(x$tol), "\n")
    invisible(x)
}
