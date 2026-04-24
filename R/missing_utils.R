#' Missing data utilities
#'
#' Internal helper functions used by antedependence model fitting when the input
#' contains missing values (NA). These functions are not exported.
#'
#' @importFrom stats complete.cases
#' @keywords internal
#' @noRd
NULL

.validate_missing <- function(y) {
    if (!is.matrix(y)) y <- as.matrix(y)

    n_subjects <- nrow(y)
    n_time <- ncol(y)

    if (!any(is.na(y))) {
        return(list(
            has_missing = FALSE,
            n_complete = n_subjects,
            n_intermittent = 0,
            patterns = rep("complete", n_subjects),
            pct_missing = 0
        ))
    }

    patterns <- apply(y, 1, function(row) {
        obs <- which(!is.na(row))

        if (length(obs) == 0) return("all_missing")
        if (length(obs) == n_time) return("complete")

        if (length(obs) >= 2 && all(diff(obs) == 1)) {
            if (obs[1] == 1) return("dropout")
            if (obs[length(obs)] == n_time) return("dropin")
            return("monotone_middle")
        }

        "intermittent"
    })

    list(
        has_missing = TRUE,
        n_complete = sum(patterns == "complete"),
        n_intermittent = sum(patterns == "intermittent"),
        patterns = patterns,
        pct_missing = mean(is.na(y)) * 100
    )
}

.get_truncation_bound <- function(y, subject_idx, time_idx,
                                  buffer = 1.2, min_bound = 10) {
    if (!is.matrix(y)) y <- as.matrix(y)

    subject_obs <- y[subject_idx, ]
    subj_vals <- subject_obs[!is.na(subject_obs)]
    subject_max <- if (length(subj_vals) > 0) max(subj_vals) else 0

    time_obs <- y[-subject_idx, time_idx]
    time_vals <- time_obs[!is.na(time_obs)]
    time_max <- if (length(time_vals) > 0) max(time_vals) else 0

    raw_bound <- buffer * max(subject_max, time_max)

    if (!is.finite(raw_bound) || raw_bound < min_bound) {
        raw_bound <- min_bound
    }

    floor(raw_bound)
}

.safe_log <- function(x) {
    out <- rep(NA_real_, length(x))

    is_zero <- (x == 0)
    is_neg <- (x < 0)
    is_pos <- (x > 0)

    out[is_zero] <- 0
    out[is_neg] <- -Inf
    out[is_pos] <- log(x[is_pos])

    out
}

.warn_intermittent_missing <- function(y, missing_info, max_intermittent = 3) {
    if (!is.matrix(y)) y <- as.matrix(y)

    if (is.null(missing_info$n_intermittent) || missing_info$n_intermittent <= 0) {
        return(invisible(NULL))
    }

    miss_per_subject <- rowSums(is.na(y))
    max_per_subject <- if (length(miss_per_subject) > 0) max(miss_per_subject) else 0

    if (is.finite(max_per_subject) && max_per_subject > max_intermittent) {
        warning(
            paste0(
                "Subject(s) with more than ", max_intermittent,
                " missing values. Computation may be slow. ",
                "Consider na_action = 'complete'."
            ),
            call. = FALSE
        )
    }

    invisible(NULL)
}

.extract_complete_cases <- function(y, blocks = NULL, warn = TRUE) {
    if (!is.matrix(y)) y <- as.matrix(y)

    complete_rows <- stats::complete.cases(y)
    n_complete <- sum(complete_rows)
    n_total <- nrow(y)

    if (n_complete == 0) {
        stop("No complete cases available for initialization.", call. = FALSE)
    }

    if (isTRUE(warn) && n_complete < n_total) {
        pct_complete <- round(100 * n_complete / n_total)
        warning(
            paste0(
                "Only ", n_complete, "/", n_total,
                " subjects have complete data (", pct_complete, "%).\n",
                "Consider using na_action = 'em' for less information loss."
            ),
            call. = FALSE
        )
    }

    y_complete <- y[complete_rows, , drop = FALSE]
    blocks_complete <- if (!is.null(blocks)) blocks[complete_rows] else NULL

    list(
        y = y_complete,
        blocks = blocks_complete,
        n_complete = n_complete,
        n_removed = n_total - n_complete
    )
}

.normalize_blocks <- function(blocks, n_subjects) {
    if (is.null(blocks)) {
        return(list(
            blocks_id = rep(1L, n_subjects),
            block_levels = "1",
            n_blocks = 1L
        ))
    }

    if (length(blocks) != n_subjects) {
        stop("blocks must have length equal to number of subjects (", n_subjects, ")", call. = FALSE)
    }
    if (any(is.na(blocks))) {
        stop("blocks must not contain missing values", call. = FALSE)
    }

    block_fac <- factor(blocks)
    list(
        blocks_id = as.integer(block_fac),
        block_levels = levels(block_fac),
        n_blocks = nlevels(block_fac)
    )
}
