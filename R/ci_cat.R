# ci_cat.R - Confidence intervals for categorical antedependence model parameters

#' Confidence intervals for fitted categorical AD models
#'
#' Computes Wald-based confidence intervals for the transition probability
#' parameters of a fitted categorical antedependence model.
#'
#' @param fit A fitted model object of class \code{"cat_fit"} from \code{fit_cat()}.
#' @param y Optional data matrix. If NULL, \code{fit$cell_counts} is used
#'   (observed counts for closed-form fits; expected counts for EM fits).
#' @param level Confidence level (default 0.95).
#' @param parameters Which parameters to compute CIs for: "all" (default),
#'   "marginal", or "transition".
#'
#' @return A list of class \code{"cat_ci"} containing:
#'   \item{marginal}{Data frame of CIs for marginal parameters (if requested)}
#'   \item{transition}{List of data frames of CIs for transition parameters (if requested)}
#'   \item{level}{Confidence level used}
#'   \item{settings}{Model settings from fit}
#'
#' @details
#' Confidence intervals are computed using the Wald method based on the
#' asymptotic normality of maximum likelihood estimators.
#'
#' For a probability estimate \eqn{\hat{\pi}} based on count N, the standard error is:
#' \deqn{SE(\hat{\pi}) = \sqrt{\frac{\hat{\pi}(1-\hat{\pi})}{N}}}
#'
#' For conditional probabilities \eqn{\hat{\pi}_{j|i}} based on conditioning count
#' \eqn{N_i}, the standard error is:
#' \deqn{SE(\hat{\pi}_{j|i}) = \sqrt{\frac{\hat{\pi}_{j|i}(1-\hat{\pi}_{j|i})}{N_i}}}
#'
#' The confidence interval is then:
#' \deqn{\hat{\pi} \pm z_{\alpha/2} \times SE(\hat{\pi})}
#'
#' Note: CIs are truncated to the interval from 0 to 1 when they exceed these bounds.
#'
#' Missing-data fits with \code{na_action = "marginalize"} are not currently
#' supported because observed cell counts are not stored for that path.
#'
#' @examples
#' \donttest{
#' # Fit a model
#' set.seed(123)
#' y <- simulate_cat(200, 5, order = 1, n_categories = 2)
#' fit <- fit_cat(y, order = 1)
#'
#' # Compute confidence intervals
#' ci <- ci_cat(fit)
#' print(ci)
#'
#' # Just marginal CIs
#' ci_marg <- ci_cat(fit, parameters = "marginal")
#' }
#'
#' @references
#' Agresti, A. (2013). Categorical Data Analysis (3rd ed.). Wiley.
#'
#' @seealso \code{\link{fit_cat}}
#'
#' @export
ci_cat <- function(fit, y = NULL, level = 0.95, parameters = "all") {
  
  if (!inherits(fit, "cat_fit")) {
    stop("fit must be a cat_fit object from fit_cat()")
  }

  fit_na_action <- NULL
  if (!is.null(fit$settings$na_action_effective)) {
    fit_na_action <- fit$settings$na_action_effective
  } else if (!is.null(fit$settings$na_action)) {
    fit_na_action <- fit$settings$na_action
  }
  if (identical(fit_na_action, "marginalize") || is.null(fit$cell_counts)) {
    stop(
      "ci_cat currently supports complete-data fits only. Missing-data CAT confidence intervals are not implemented yet; refit with na_action = 'complete'.",
      call. = FALSE
    )
  }

  if (!is.null(y) && anyNA(y)) {
    stop(
      "ci_cat currently supports complete y only. Missing-data CAT confidence intervals are not implemented yet.",
      call. = FALSE
    )
  }
  
  if (level <= 0 || level >= 1) {
    stop("level must be between 0 and 1")
  }
  
  parameters <- match.arg(parameters, c("all", "marginal", "transition"))
  
  # Get z critical value
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  
  # Extract settings
  p <- fit$settings$order
  c <- fit$settings$n_categories
  n_time <- fit$settings$n_time
  n_subjects <- fit$settings$n_subjects
  homogeneous <- fit$settings$homogeneous
  n_blocks <- fit$settings$n_blocks
  
  # Determine if we have heterogeneous groups
  if (!homogeneous && n_blocks > 1) {
    # Handle heterogeneous case
    marginal_ci <- vector("list", n_blocks)
    transition_ci <- vector("list", n_blocks)
    
    for (g in seq_len(n_blocks)) {
      marg_g <- fit$marginal[[g]]
      trans_g <- fit$transition[[g]]
      counts_g <- fit$cell_counts[[g]]
      
      # Get sample size for this block
      if (!is.null(fit$settings$blocks)) {
        n_g <- sum(fit$settings$blocks == g)
      } else {
        n_g <- n_subjects / n_blocks  # Approximate
      }
      
      if (parameters %in% c("all", "marginal")) {
        marginal_ci[[g]] <- .compute_marginal_ci(marg_g, counts_g, p, c, n_time, n_g, z, level)
      }
      
      if (parameters %in% c("all", "transition") && p >= 1) {
        transition_ci[[g]] <- .compute_transition_ci(trans_g, counts_g, p, c, n_time, z, level)
      }
    }
    names(marginal_ci) <- paste0("block_", seq_len(n_blocks))
    names(transition_ci) <- paste0("block_", seq_len(n_blocks))
    
  } else {
    # Homogeneous case
    marginal_ci <- NULL
    transition_ci <- NULL
    
    if (parameters %in% c("all", "marginal")) {
      marginal_ci <- .compute_marginal_ci(fit$marginal, fit$cell_counts, 
                                           p, c, n_time, n_subjects, z, level)
    }
    
    if (parameters %in% c("all", "transition") && p >= 1) {
      transition_ci <- .compute_transition_ci(fit$transition, fit$cell_counts,
                                               p, c, n_time, z, level)
    }
  }
  
  out <- list(
    marginal = marginal_ci,
    transition = transition_ci,
    level = level,
    settings = fit$settings
  )
  
  class(out) <- "cat_ci"
  out
}


#' Compute CIs for marginal parameters
#'
#' @param marginal Marginal parameter list from fit
#' @param cell_counts Cell counts from fit
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param n_subjects Number of subjects
#' @param z Z critical value
#' @param level Confidence level
#'
#' @return Data frame with CIs
#'
#' @keywords internal
.compute_marginal_ci <- function(marginal, cell_counts, p, c, n_time, n_subjects, z, level) {
  
  results <- list()
  
  if (p == 0) {
    # Independence model: marginal at each time point
    for (k in seq_len(n_time)) {
      probs <- marginal[[k]]
      
      for (cat in seq_len(c)) {
        prob <- probs[cat]
        se <- sqrt(prob * (1 - prob) / n_subjects)
        lower <- max(0, prob - z * se)
        upper <- min(1, prob + z * se)
        
        results[[length(results) + 1]] <- data.frame(
          parameter = paste0("pi_t", k, "_cat", cat),
          time = k,
          category = cat,
          estimate = prob,
          se = se,
          lower = lower,
          upper = upper,
          level = level,
          stringsAsFactors = FALSE
        )
      }
    }
  } else {
    # AD(p): marginal at t1
    probs_t1 <- marginal[["t1"]]
    
    for (cat in seq_len(c)) {
      prob <- probs_t1[cat]
      se <- sqrt(prob * (1 - prob) / n_subjects)
      lower <- max(0, prob - z * se)
      upper <- min(1, prob + z * se)
      
      results[[length(results) + 1]] <- data.frame(
        parameter = paste0("pi_t1_cat", cat),
        time = 1,
        category = cat,
        estimate = prob,
        se = se,
        lower = lower,
        upper = upper,
        level = level,
        stringsAsFactors = FALSE
      )
    }
    
    # For order >= 2, also have conditionals in marginal
    if (p >= 2 && n_time >= 2) {
      for (k in 2:min(p, n_time)) {
        trans_name <- paste0("t", k, "_given_1to", k - 1)
        if (!is.null(marginal[[trans_name]])) {
          # These are conditional probabilities - need conditioning counts
          # For now, skip detailed CIs for these (complex structure)
        }
      }
    }
  }
  
  if (length(results) > 0) {
    do.call(rbind, results)
  } else {
    NULL
  }
}


#' Compute CIs for transition parameters
#'
#' @param transition Transition parameter list from fit
#' @param cell_counts Cell counts from fit
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param z Z critical value
#' @param level Confidence level
#'
#' @return List of data frames with CIs for each time point
#'
#' @keywords internal
.compute_transition_ci <- function(transition, cell_counts, p, c, n_time, z, level) {
  
  results_list <- list()
  
  for (k in (p + 1):n_time) {
    trans_name <- paste0("t", k)
    trans_k <- transition[[trans_name]]
    
    # Get corresponding counts
    count_name <- paste0("t", k - p, "_to_t", k)
    counts_k <- cell_counts[[count_name]]
    
    if (is.null(trans_k) || is.null(counts_k)) next
    
    results <- list()
    
    if (p == 1) {
      # Order 1: trans_k is c x c matrix
      # trans_k[i, j] = P(Y_k = j | Y_{k-1} = i)
      # Conditioning count is sum over j of counts_k[i, j]
      
      for (i in seq_len(c)) {
        n_cond <- sum(counts_k[i, ])  # Conditioning count
        
        for (j in seq_len(c)) {
          prob <- trans_k[i, j]
          
          if (n_cond > 0) {
            se <- sqrt(prob * (1 - prob) / n_cond)
          } else {
            se <- NA
          }
          
          lower <- if (is.na(se)) NA else max(0, prob - z * se)
          upper <- if (is.na(se)) NA else min(1, prob + z * se)
          
          results[[length(results) + 1]] <- data.frame(
            parameter = paste0("pi_t", k, "_", j, "|", i),
            time = k,
            from = i,
            to = j,
            estimate = prob,
            se = se,
            lower = lower,
            upper = upper,
            n_cond = n_cond,
            level = level,
            stringsAsFactors = FALSE
          )
        }
      }
      
    } else if (p == 2) {
      # Order 2: trans_k is c x c x c array
      # trans_k[i, j, k] = P(Y_k = k | Y_{k-2} = i, Y_{k-1} = j)
      
      for (i in seq_len(c)) {
        for (j in seq_len(c)) {
          n_cond <- sum(counts_k[i, j, ])  # Conditioning count
          
          for (cat in seq_len(c)) {
            prob <- trans_k[i, j, cat]
            
            if (n_cond > 0) {
              se <- sqrt(prob * (1 - prob) / n_cond)
            } else {
              se <- NA
            }
            
            lower <- if (is.na(se)) NA else max(0, prob - z * se)
            upper <- if (is.na(se)) NA else min(1, prob + z * se)
            
            results[[length(results) + 1]] <- data.frame(
              parameter = paste0("pi_t", k, "_", cat, "|", i, ",", j),
              time = k,
              from1 = i,
              from2 = j,
              to = cat,
              estimate = prob,
              se = se,
              lower = lower,
              upper = upper,
              n_cond = n_cond,
              level = level,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(results) > 0) {
      results_list[[trans_name]] <- do.call(rbind, results)
    }
  }
  
  results_list
}


#' Print method for cat_ci objects
#'
#' @param x A cat_ci object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns \code{x}.
#' @export
print.cat_ci <- function(x, ...) {
  cat("Confidence Intervals for Categorical AD Model\n")
  cat("==============================================\n\n")
  
  cat("Confidence level:", x$level * 100, "%\n")
  cat("Order:", x$settings$order, "\n")
  cat("Categories:", x$settings$n_categories, "\n")
  cat("Time points:", x$settings$n_time, "\n\n")
  
  if (!is.null(x$marginal)) {
    cat("Marginal probabilities:\n")
    if (is.data.frame(x$marginal)) {
      print(x$marginal[, c("parameter", "estimate", "se", "lower", "upper")], 
            row.names = FALSE, digits = 4)
    } else {
      # Heterogeneous case
      for (g in names(x$marginal)) {
        cat("\n", g, ":\n", sep = "")
        if (!is.null(x$marginal[[g]])) {
          print(x$marginal[[g]][, c("parameter", "estimate", "se", "lower", "upper")],
                row.names = FALSE, digits = 4)
        }
      }
    }
    cat("\n")
  }
  
  if (!is.null(x$transition)) {
    cat("Transition probabilities:\n")
    if (is.list(x$transition) && !is.null(names(x$transition))) {
      # Check if it's per-time or per-block
      if (grepl("^t[0-9]+$", names(x$transition)[1])) {
        # Per-time
        for (t_name in names(x$transition)) {
          cat("\n", t_name, ":\n", sep = "")
          df <- x$transition[[t_name]]
          cols_to_show <- intersect(c("parameter", "estimate", "se", "lower", "upper", "n_cond"),
                                    names(df))
          print(df[, cols_to_show], row.names = FALSE, digits = 4)
        }
      } else {
        # Per-block (heterogeneous)
        for (g in names(x$transition)) {
          cat("\n", g, ":\n", sep = "")
          for (t_name in names(x$transition[[g]])) {
            cat("  ", t_name, ":\n", sep = "")
            df <- x$transition[[g]][[t_name]]
            cols_to_show <- intersect(c("parameter", "estimate", "se", "lower", "upper"),
                                      names(df))
            print(df[, cols_to_show], row.names = FALSE, digits = 4)
          }
        }
      }
    }
  }
  
  invisible(x)
}


#' Summary method for cat_ci objects
#'
#' @param object A cat_ci object
#' @param ... Additional arguments (ignored)
#'
#' @return A data frame summarizing all CIs
#'
#' @export
summary.cat_ci <- function(object, ...) {
  
  all_cis <- list()
  
  # Collect marginal CIs
  if (!is.null(object$marginal)) {
    if (is.data.frame(object$marginal)) {
      object$marginal$type <- "marginal"
      all_cis[[length(all_cis) + 1]] <- object$marginal
    } else {
      for (g in names(object$marginal)) {
        if (!is.null(object$marginal[[g]])) {
          df <- object$marginal[[g]]
          df$type <- "marginal"
          df$block <- g
          all_cis[[length(all_cis) + 1]] <- df
        }
      }
    }
  }
  
  # Collect transition CIs
  if (!is.null(object$transition)) {
    if (is.list(object$transition) && !is.null(names(object$transition))) {
      if (grepl("^t[0-9]+$", names(object$transition)[1])) {
        # Per-time
        for (t_name in names(object$transition)) {
          df <- object$transition[[t_name]]
          df$type <- "transition"
          all_cis[[length(all_cis) + 1]] <- df
        }
      } else {
        # Per-block
        for (g in names(object$transition)) {
          for (t_name in names(object$transition[[g]])) {
            df <- object$transition[[g]][[t_name]]
            df$type <- "transition"
            df$block <- g
            all_cis[[length(all_cis) + 1]] <- df
          }
        }
      }
    }
  }
  
  if (length(all_cis) > 0) {
    # Combine, keeping only common columns
    result <- do.call(rbind, lapply(all_cis, function(df) {
      df[, intersect(names(df), c("parameter", "type", "estimate", "se", "lower", "upper", "level"))]
    }))
    return(result)
  }
  
  NULL
}
