#' Predict from an INAD fit
#'
#' Generates plug-in predictive means or Monte Carlo sample paths from a fitted
#' integer-valued antedependence model. Version 1 supports order-1 INAD fits and
#' uses fitted maximum-likelihood parameters without integrating over parameter
#' uncertainty.
#'
#' @param object An \code{inad_fit} object returned by \code{\link{fit_inad}}.
#' @param newdata Optional nonnegative integer matrix containing observed
#'   histories for the subjects to forecast. If \code{NULL}, the training data
#'   stored on \code{object} are used through time \code{n_time - h}, so the
#'   method predicts the last \code{h} training occasions.
#' @param blocks Optional block/group vector for the rows of \code{newdata}.
#'   Required when the fitted model has nonzero block effects and
#'   \code{newdata} is supplied. If \code{newdata = NULL}, fitted training
#'   blocks are used.
#' @param h Positive integer forecast horizon.
#' @param type Prediction type. \code{"mean"} returns the Monte Carlo
#'   predictive mean. \code{"sample"} returns raw simulated forecast paths.
#' @param n_sims Positive integer number of Monte Carlo paths.
#' @param seed Optional random seed for reproducible Monte Carlo forecasts.
#' @param ... Unused.
#'
#' @return For \code{type = "mean"}, an \code{n_test x h} numeric matrix. For
#'   \code{type = "sample"}, an \code{n_test x h x n_sims} integer array.
#'
#' @details
#' Predictions are conditional on the fitted parameter estimates, following the
#' usual plug-in convention used by many \code{predict()} methods. They do not
#' account for parameter-estimation uncertainty. When \code{newdata} is supplied,
#' the forecast starts at \code{ncol(newdata) + 1}; therefore
#' \code{ncol(newdata) + h} must not exceed the fitted time grid.
#'
#' @examples
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
#' fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois")
#' predict(fit, h = 1, type = "mean", n_sims = 20, seed = 1)
#'
#' @export
predict.inad_fit <- function(object,
                             newdata = NULL,
                             blocks = NULL,
                             h = 1L,
                             type = c("mean", "sample"),
                             n_sims = 1000L,
                             seed = NULL,
                             ...) {
    type <- match.arg(type)
    h <- .inad_predict_positive_integer(h, "h")
    n_sims <- .inad_predict_positive_integer(n_sims, "n_sims")

    settings <- object$settings
    order <- as.integer(settings$order)
    if (length(order) != 1L || is.na(order)) {
        stop("object$settings$order is not available.", call. = FALSE)
    }
    if (order != 1L) {
        stop("predict.inad_fit() currently supports order = 1 fits only.", call. = FALSE)
    }

    n_time <- as.integer(settings$n_time)
    if (length(n_time) != 1L || is.na(n_time)) n_time <- length(object$theta)

    if (is.null(newdata)) {
        if (is.null(object$.y)) {
            stop("newdata is required because the fit does not contain training data.", call. = FALSE)
        }
        if (h >= ncol(object$.y)) {
            stop("h must be smaller than the fitted number of time points when newdata = NULL.", call. = FALSE)
        }
        history <- object$.y[, seq_len(ncol(object$.y) - h), drop = FALSE]
        blocks_pred <- if (is.null(blocks)) object$.blocks else
            .inad_predict_normalize_blocks(blocks, nrow(history), object)
    } else {
        history <- .inad_predict_history_matrix(newdata)
        blocks_pred <- .inad_predict_normalize_blocks(blocks, nrow(history), object)
    }

    samples <- .simulate_inad_forward(
        fit = object,
        history = history,
        blocks = blocks_pred,
        start_time = ncol(history) + 1L,
        h = h,
        n_sims = n_sims,
        seed = seed
    )

    if (type == "sample") return(samples)

    out <- apply(samples, c(1L, 2L), mean)
    dim(out) <- c(dim(samples)[1L], dim(samples)[2L])
    rownames(out) <- rownames(history)
    colnames(out) <- paste0("h", seq_len(h))
    out
}

#' @keywords internal
.simulate_inad_forward <- function(fit,
                                   history,
                                   blocks = NULL,
                                   start_time = ncol(history) + 1L,
                                   h = 1L,
                                   n_sims = 1000L,
                                   seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    history <- .inad_predict_history_matrix(history)

    n_test <- nrow(history)
    n_obs <- ncol(history)
    start_time <- as.integer(start_time)
    h <- .inad_predict_positive_integer(h, "h")
    n_sims <- .inad_predict_positive_integer(n_sims, "n_sims")

    if (length(start_time) != 1L || is.na(start_time) || start_time != n_obs + 1L) {
        stop("start_time must equal ncol(history) + 1.", call. = FALSE)
    }

    settings <- fit$settings
    order <- as.integer(settings$order)
    thinning <- settings$thinning
    innovation <- settings$innovation
    n_time <- as.integer(settings$n_time)
    if (length(n_time) != 1L || is.na(n_time)) n_time <- length(fit$theta)

    if (order != 1L) {
        stop(".simulate_inad_forward() currently supports order = 1 only.", call. = FALSE)
    }
    if (start_time + h - 1L > n_time) {
        stop("Requested horizon exceeds the fitted time grid.", call. = FALSE)
    }

    tau <- fit$tau
    if (is.null(tau)) tau <- 0
    tau <- as.numeric(tau)
    has_tau <- length(tau) > 1L && any(abs(tau[-1L]) > 1e-12)
    if (has_tau && is.null(blocks)) {
        stop("blocks must be supplied when the fit has nonzero block effects.", call. = FALSE)
    }
    if (is.null(blocks)) blocks <- rep(1L, n_test)
    blocks <- as.integer(blocks)
    if (length(blocks) != n_test) {
        stop("blocks must have length nrow(history).", call. = FALSE)
    }
    if (any(blocks < 1L | blocks > length(tau))) {
        stop("blocks contain levels outside fit$tau.", call. = FALSE)
    }

    alpha <- .inad_predict_time_vector(fit$alpha, n_time, "alpha")
    theta <- .inad_predict_time_vector(fit$theta, n_time, "theta")

    nb_inno_size <- NULL
    if (innovation == "nbinom") {
        if (is.null(fit$nb_inno_size)) {
            stop("nb_inno_size is required for negative-binomial innovation.", call. = FALSE)
        }
        nb_inno_size <- .inad_predict_time_vector(fit$nb_inno_size, n_time, "nb_inno_size")
    }

    draw_innovation <- function(t_index) {
        eff_mean <- .inad_effective_innovation_mean(theta[t_index], tau[blocks], innovation)
        if (any(!is.finite(eff_mean)) || any(eff_mean <= 0)) {
            stop("Non-positive effective innovation mean.", call. = FALSE)
        }

        if (innovation == "pois") {
            return(stats::rpois(n_test, lambda = eff_mean))
        }
        if (innovation == "bell") {
            eff_theta <- .bell_theta_from_mean(eff_mean)
            return(as.integer(vapply(eff_theta, function(th) rbell(1L, theta = th), integer(1L))))
        }
        stats::rnbinom(n_test, size = nb_inno_size[t_index], mu = eff_mean)
    }

    draw_thinning <- function(y_prev, alpha_t) {
        if (alpha_t <= 0) return(integer(n_test))
        if (thinning == "binom") {
            return(stats::rbinom(n_test, size = y_prev, prob = alpha_t))
        }
        if (thinning == "pois") {
            return(stats::rpois(n_test, lambda = alpha_t * y_prev))
        }

        out <- integer(n_test)
        idx <- y_prev > 0L
        if (any(idx)) {
            out[idx] <- stats::rnbinom(
                n = sum(idx),
                size = y_prev[idx],
                prob = 1 / (1 + alpha_t)
            )
        }
        out
    }

    out <- array(NA_integer_, dim = c(n_test, h, n_sims))
    for (ss in seq_len(n_sims)) {
        y_prev <- history[, n_obs]
        for (step in seq_len(h)) {
            t_index <- start_time + step - 1L
            y_next <- draw_thinning(y_prev, alpha[t_index]) + draw_innovation(t_index)
            out[, step, ss] <- as.integer(y_next)
            y_prev <- y_next
        }
    }
    out
}

#' @keywords internal
.inad_predict_history_matrix <- function(x) {
    if (is.null(x)) stop("history cannot be NULL.", call. = FALSE)
    if (is.vector(x) && !is.list(x)) x <- matrix(x, nrow = 1L)
    if (!is.matrix(x)) x <- as.matrix(x)
    if (nrow(x) < 1L || ncol(x) < 1L) {
        stop("newdata/history must have at least one row and one column.", call. = FALSE)
    }
    if (any(is.na(x))) {
        stop("Prediction histories cannot contain missing values.", call. = FALSE)
    }
    if (any(!is.finite(x) | x < 0 | x != floor(x))) {
        stop("Prediction histories must be nonnegative integer values.", call. = FALSE)
    }
    storage.mode(x) <- "integer"
    x
}

#' @keywords internal
.inad_predict_positive_integer <- function(x, name) {
    x <- as.integer(x)
    if (length(x) != 1L || is.na(x) || x < 1L) {
        stop(name, " must be a positive integer.", call. = FALSE)
    }
    x
}

#' @keywords internal
.inad_predict_time_vector <- function(x, n_time, name) {
    out <- as.numeric(x)
    if (length(out) == 1L) out <- rep(out, n_time)
    if (length(out) != n_time) {
        stop(name, " length does not match the fitted time grid.", call. = FALSE)
    }
    out
}

#' @keywords internal
.inad_predict_normalize_blocks <- function(blocks, n, fit) {
    tau <- fit$tau
    if (is.null(tau)) tau <- 0
    tau <- as.numeric(tau)
    has_tau <- length(tau) > 1L && any(abs(tau[-1L]) > 1e-12)
    if (is.null(blocks)) {
        if (has_tau) {
            stop("blocks must be supplied when the fit has nonzero block effects.", call. = FALSE)
        }
        return(NULL)
    }
    if (length(blocks) != n) {
        stop("blocks must have length nrow(newdata).", call. = FALSE)
    }

    levels <- fit$settings$block_levels
    if (!is.null(levels) && !is.numeric(blocks) && !is.integer(blocks)) {
        out <- match(as.character(blocks), as.character(levels))
        if (any(is.na(out))) stop("blocks contain levels not present in the fitted model.", call. = FALSE)
        return(as.integer(out))
    }

    out <- as.integer(blocks)
    if (any(is.na(out))) stop("blocks must be coercible to integer block ids.", call. = FALSE)
    out
}
