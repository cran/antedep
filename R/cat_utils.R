# cat_utils.R - Internal utilities for categorical antedependence models

#' Validate categorical data matrix
#'
#' @param y Data matrix
#' @param n_categories Number of categories (NULL to infer)
#' @param allow_na Logical; if TRUE, allow missing values.
#'
#' @return List with validated y and n_categories
#'
#' @keywords internal
.validate_y_cat <- function(y, n_categories = NULL, allow_na = FALSE) {
  # Convert to matrix if needed

  if (!is.matrix(y)) y <- as.matrix(y)
  
  # Check for numeric/integer
  if (!is.numeric(y)) stop("y must be numeric")
  
  # Check for integer values
  if (any(y != floor(y), na.rm = TRUE)) {
    stop("y must be integer-valued (category codes)")
  }
  
  # Check for positive values
  if (any(y < 1, na.rm = TRUE)) {
    stop("y must have category codes >= 1")
  }
  
  # Check for NA
  if (!allow_na && any(is.na(y))) {
    stop("Missing data not supported in this call; y must be complete")
  }
  
  # Infer number of categories
  c_inferred <- max(y, na.rm = TRUE)
  if (!is.finite(c_inferred)) {
    stop("Unable to infer categories: all values are NA")
  }
  
  if (!is.null(n_categories)) {
    if (n_categories < c_inferred) {
      stop("n_categories (", n_categories, ") is smaller than max observed category (", 
           c_inferred, ")")
    }
    c_inferred <- n_categories
  }
  
  # Convert to integer storage
  storage.mode(y) <- "integer"
  
  list(y = y, n_categories = c_inferred)
}


#' Validate blocks parameter
#'
#' @param blocks Block membership vector
#' @param n_subjects Number of subjects
#'
#' @return Validated blocks vector (integer, 1-indexed)
#'
#' @keywords internal
.validate_blocks_cat <- function(blocks, n_subjects) {
  if (is.null(blocks)) {
    return(rep(1L, n_subjects))
  }
  
  if (length(blocks) != n_subjects) {
    stop("blocks must have length equal to number of subjects (", n_subjects, ")")
  }
  
  blocks <- as.integer(blocks)
  
  if (any(blocks < 1)) {
    stop("blocks must have values >= 1")
  }
  
  # Ensure contiguous from 1
  unique_blocks <- sort(unique(blocks))
  if (!identical(unique_blocks, seq_along(unique_blocks))) {
    # Remap to contiguous
    blocks <- match(blocks, unique_blocks)
  }
  
  blocks
}


#' Count free parameters for AD(p) categorical model
#'
#' @param order AD order p
#' @param n_categories Number of categories c
#' @param n_time Number of time points n
#' @param n_blocks Number of blocks/groups
#' @param homogeneous Whether parameters are shared across blocks
#'
#' @return Number of free parameters
#'
#' @keywords internal
.count_params_cat <- function(order, n_categories, n_time, n_blocks = 1, 
                               homogeneous = TRUE) {
  c <- n_categories
  n <- n_time
  p <- order
  
  if (p >= n) {
    stop("order must be less than n_time")
  }
  
  # Number of free parameters per population
  # For AD(p): (c-1) * sum_{k=1}^{n} c^{min(k-1, p)}
  # = (c-1) * [c^0 + c^1 + ... + c^{p-1} + (n-p) * c^p]  for p >= 1
  # = (c-1) * n  for p = 0
  
  if (p == 0) {
    n_params_per_pop <- (c - 1) * n
  } else {
    # Sum of c^0 + c^1 + ... + c^{p-1} = (c^p - 1) / (c - 1)
    sum_geometric <- (c^p - 1) / (c - 1)
    n_params_per_pop <- (c - 1) * (sum_geometric + (n - p) * c^p)
  }
  
  # Multiply by number of populations if heterogeneous
  if (homogeneous || n_blocks == 1) {
    return(as.integer(n_params_per_pop))
  } else {
    return(as.integer(n_params_per_pop * n_blocks))
  }
}


#' Get all category combinations of given length
#'
#' @param n_categories Number of categories c
#' @param length Number of positions
#'
#' @return Matrix where each row is a combination (c^length rows, length columns)
#'
#' @keywords internal
.get_combinations_cat <- function(n_categories, length) {
  if (length == 0) {
    return(matrix(nrow = 1, ncol = 0))
  }
  
  if (length == 1) {
    return(matrix(1:n_categories, ncol = 1))
  }
  
  # Use expand.grid and convert to matrix
  args <- rep(list(1:n_categories), length)
  combos <- expand.grid(args)
  as.matrix(combos)
}


#' Count cells for contingency table
#'
#' Counts occurrences of each combination of categories at specified time indices.
#'
#' @param y Data matrix (n_subjects x n_time)
#' @param time_indices Integer vector of time indices to count
#' @param n_categories Number of categories
#' @param subject_mask Logical vector indicating which subjects to include (NULL = all)
#'
#' @return Array of counts with dimensions c x c x ... (length = length(time_indices))
#'
#' @keywords internal
.count_cells_cat <- function(y, time_indices, n_categories, subject_mask = NULL) {
  n_subjects <- nrow(y)
  n_positions <- length(time_indices)
  c <- n_categories
  
  # Apply subject mask
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  
  n_sub <- nrow(y_sub)
  
  if (n_positions == 0) {
    return(array(n_sub, dim = 1))
  }
  
  # Extract relevant columns
  y_cols <- y_sub[, time_indices, drop = FALSE]
  
  # Create count array
  dims <- rep(c, n_positions)
  counts <- array(0L, dim = dims)
  
  # Count each combination
  for (i in seq_len(n_sub)) {
    # Get the index into the array
    idx <- as.list(y_cols[i, ])
    counts[matrix(unlist(idx), nrow = 1)] <- counts[matrix(unlist(idx), nrow = 1)] + 1L
  }
  
  counts
}


#' Count cells efficiently using table()
#'
#' Alternative implementation using R's table function.
#'
#' @param y Data matrix (n_subjects x n_time)
#' @param time_indices Integer vector of time indices
#' @param n_categories Number of categories
#' @param subject_mask Logical vector for subject selection
#'
#' @return Array of counts
#'
#' @keywords internal
.count_cells_table_cat <- function(y, time_indices, n_categories, subject_mask = NULL) {
  c <- n_categories
  n_positions <- length(time_indices)
  
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  
  if (n_positions == 0) {
    return(array(nrow(y_sub), dim = 1))
  }
  
  # Extract columns and convert to factors
  y_cols <- y_sub[, time_indices, drop = FALSE]
  
  # Create list of factors with consistent levels
  factors <- lapply(seq_len(n_positions), function(j) {
    factor(y_cols[, j], levels = 1:c)
  })
  
  # Use table to count
  counts <- table(factors)
  
  # Convert to array (table already returns array-like object)
  array(as.integer(counts), dim = rep(c, n_positions))
}


#' Convert counts to probabilities (with safe division)
#'
#' @param counts Array of counts
#' @param margin Which margins to condition on (sum over). NULL = no conditioning (joint).
#'
#' @return Array of probabilities
#'
#' @keywords internal
.counts_to_probs_cat <- function(counts, margin = NULL) {
  if (is.null(margin)) {
    # Joint probability: divide by total
    total <- sum(counts)
    if (total == 0) {
      return(counts * 0)  # All zeros
    }
    return(counts / total)
  }
  
  # Conditional probability: divide by marginal over specified dimensions
  # margin specifies which dimensions to keep (condition on)
  marginal <- apply(counts, margin, sum)
  
  # Divide each slice by its marginal
  # This is tricky with arrays - use sweep or manual
  probs <- counts
  dims <- dim(counts)
  n_dims <- length(dims)
  other_dims <- setdiff(seq_len(n_dims), margin)
  
  if (length(other_dims) == 1) {
    # Simple case: condition on all but one dimension
    # The last dimension (other_dims) is what we're computing P(Y_last | rest)
    for (idx in seq_len(prod(dims[margin]))) {
      # Convert linear index to subscripts for margin dimensions
      subs <- arrayInd(idx, dims[margin])
      
      # Build index for full array
      full_idx <- vector("list", n_dims)
      full_idx[margin] <- as.list(subs)
      full_idx[other_dims] <- list(TRUE)
      
      # Get denominator
      denom <- do.call(`[`, c(list(marginal), as.list(subs)))
      
      if (denom == 0) {
        # Set to 0 (or could use uniform 1/c)
        do.call(`[<-`, c(list(probs), full_idx, list(0)))
      } else {
        slice <- do.call(`[`, c(list(counts), full_idx))
        do.call(`[<-`, c(list(probs), full_idx, list(slice / denom)))
      }
    }
  } else {
    # More complex case - use sweep
    # Actually, let's just loop over all conditioning combinations
    cond_combos <- .get_combinations_cat(dims[margin[1]], length(margin))
    # This gets complicated - let's use a simpler approach
    
    # For transition probabilities, we typically have:
    # counts[y_{k-p}, ..., y_{k-1}, y_k] and want P(y_k | y_{k-p}, ..., y_{k-1})
    # So margin = 1:(n_dims-1) and we divide by sum over last dimension
    
    probs <- counts / (marginal + (marginal == 0))  # Avoid division by zero
    probs[marginal == 0] <- 0
  }
  
  probs
}


#' Compute transition probabilities from counts (simpler version)
#'
#' Given counts array with last dimension being Y_k and earlier dimensions
#' being conditioning variables, compute P(Y_k | conditioning).
#'
#' @param counts Array of counts, last dimension is response
#'
#' @return Array of transition probabilities (same dimensions)
#'
#' @keywords internal
.counts_to_transition_probs <- function(counts) {
  dims <- dim(counts)
  n_dims <- length(dims)
  
  if (n_dims == 1) {
    # Just marginal probabilities
    total <- sum(counts)
    if (total == 0) return(rep(0, dims[1]))
    return(counts / total)
  }
  
  # Sum over last dimension to get denominator
  margin_dims <- seq_len(n_dims - 1)
  denom <- apply(counts, margin_dims, sum)
  
  # Divide - need to replicate denom to match counts dimensions
  # Use sweep for this
  probs <- sweep(counts, margin_dims, denom, FUN = function(x, d) {
    ifelse(d == 0, 0, x / d)
  })
  
  probs
}


#' Safe log function (0 * log(0) = 0)
#'
#' @param x Numeric vector or array
#'
#' @return Log of x, with 0 for x = 0
#'
#' @keywords internal
.safe_log <- function(x) {
  result <- log(x)
  result[x == 0] <- 0
  result[x < 0] <- -Inf
  result
}


#' Compute log-likelihood contribution from counts and probabilities
#'
#' Computes sum of N * log(pi) over all cells.
#'
#' @param counts Array of counts
#' @param probs Array of probabilities (same dimensions as counts)
#'
#' @return Scalar log-likelihood contribution
#'
#' @keywords internal
.loglik_contribution <- function(counts, probs) {
  # N * log(pi), with 0 * log(0) = 0
  sum(counts * .safe_log(probs))
}


#' Extract history for a subject at time k
#'
#' @param y_row Single row of data (one subject)
#' @param k Current time index
#' @param order AD order p
#'
#' @return Integer vector of length min(k-1, order) with category history
#'
#' @keywords internal
.get_history <- function(y_row, k, order) {
  if (k == 1) return(integer(0))
  
  history_length <- min(k - 1, order)
  start_idx <- k - history_length
  
  y_row[start_idx:(k - 1)]
}


#' Convert history to array index
#'
#' @param history Integer vector of category values
#' @param n_categories Number of categories
#'
#' @return Integer index (1-based) into flattened array
#'
#' @keywords internal
.history_to_index <- function(history, n_categories) {
  if (length(history) == 0) return(1L)
  
  # history[1] is earliest, history[end] is most recent
  # Treat as base-c number
  idx <- 1L
  multiplier <- 1L
  for (i in seq_along(history)) {
    idx <- idx + (history[i] - 1L) * multiplier
    multiplier <- multiplier * n_categories
  }
  idx
}


#' Print method for cat_fit objects
#'
#' @param x A cat_fit object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns \code{x}.
#' @export
print.cat_fit <- function(x, ...) {
  cat("Categorical Antedependence Model Fit\n")
  cat("====================================\n\n")
  
  cat("Order:", x$settings$order, "\n")
  cat("Categories:", x$settings$n_categories, "\n")
  cat("Time points:", x$settings$n_time, "\n")
  cat("Subjects:", x$settings$n_subjects, "\n")
  
  if (!is.null(x$settings$blocks) && !x$settings$homogeneous) {
    cat("Groups:", length(unique(x$settings$blocks)), "(heterogeneous)\n")
  } else if (!is.null(x$settings$blocks)) {
    cat("Groups:", length(unique(x$settings$blocks)), "(homogeneous)\n")
  }
  
  cat("\nLog-likelihood:", round(x$log_l, 4), "\n")
  cat("AIC:", round(x$aic, 4), "\n")
  cat("BIC:", round(x$bic, 4), "\n")
  cat("Parameters:", x$n_params, "\n")
  if (!is.null(x$n_missing)) {
    cat("Missing values:", x$n_missing, "\n")
  }
  if (!is.null(x$convergence)) {
    cat("Convergence code:", x$convergence, "\n")
  }
  
  invisible(x)
}
