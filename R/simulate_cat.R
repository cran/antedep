# simulate_cat.R - Data simulation for categorical antedependence models

#' Simulate categorical antedependence series
#'
#' Generate simulated longitudinal categorical data from an AD(p) model with
#' specified transition probabilities.
#'
#' @param n_subjects Number of subjects to simulate.
#' @param n_time Number of time points.
#' @param order Antedependence order p. Must be 0, 1, or 2. Default is 1.
#' @param n_categories Number of categories c. Default is 2 (binary).
#' @param marginal List of marginal/joint probabilities for initial time points.
#'   If NULL, uniform probabilities are used. See Details for structure.
#' @param transition List of transition probability arrays for time points
#'   k = p+1 to n. If NULL, uniform transitions are used. See Details.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership. Used with homogeneous = FALSE.
#' @param homogeneous Logical. If TRUE (default), same parameters for all subjects.
#'   If FALSE, marginal and transition should be lists indexed by block.
#' @param seed Optional random seed for reproducibility.
#'
#' @return Integer matrix with n_subjects rows and n_time columns, where each
#'   entry is a category code from 1 to c.
#'
#' @details
#' Data are simulated sequentially:
#' \enumerate{
#'   \item For k = 1: Draw Y(1) from marginal distribution
#'   \item For k = 2 to p: Draw Y(k) conditional on Y(1), ..., Y(k-1)
#'   \item For k = p+1 to n: Draw Y(k) conditional on Y(k-p), ..., Y(k-1)
#' }
#'
#' Parameter structure for marginal:
#' \itemize{
#'   \item Order 0: List with elements t1, t2, ..., tn, each a vector of length c
#'     summing to 1
#'   \item Order 1: List with element t1 (vector of length c)
#'   \item Order 2: List with t1 (vector), t2_given_1to1 (c x c matrix where rows
#'     represent conditioning values and columns represent outcomes)
#' }
#'
#' Parameter structure for transition:
#' \itemize{
#'   \item Order 0: Not used (NULL)
#'   \item Order 1: List with elements t2, t3, ..., tn, each c x c matrix where
#'     rows are previous values and columns are current values (rows sum to 1)
#'   \item Order 2: List with elements t3, t4, ..., tn, each c x c x c array where
#'     first two indices are conditioning values and third is outcome
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @examples
#' y <- simulate_cat(n_subjects = 30, n_time = 5, order = 1, n_categories = 3, seed = 1)
#' dim(y)
#'
#' @export
simulate_cat <- function(n_subjects, n_time, order = 1, n_categories = 2,
                         marginal = NULL, transition = NULL,
                         blocks = NULL, homogeneous = TRUE, seed = NULL) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  n_subjects <- as.integer(n_subjects)
  n_time <- as.integer(n_time)
  p <- as.integer(order)
  c <- as.integer(n_categories)
  
  if (n_subjects < 1) stop("n_subjects must be at least 1")
  if (n_time < 1) stop("n_time must be at least 1")
  if (p < 0) stop("order must be non-negative")
  if (p > 2) stop("order > 2 not currently supported")
  if (p >= n_time) stop("order must be less than n_time")
  if (c < 2) stop("n_categories must be at least 2")
  
  # Process blocks
  if (is.null(blocks)) {
    blocks <- rep(1L, n_subjects)
  } else {
    blocks <- .validate_blocks_cat(blocks, n_subjects)
  }
  n_blocks <- max(blocks)
  
  # Set up default parameters if not provided
  if (is.null(marginal)) {
    marginal <- .default_marginal_cat(p, c, n_time, n_blocks, homogeneous)
  }
  if (is.null(transition) && p >= 1) {
    transition <- .default_transition_cat(p, c, n_time, n_blocks, homogeneous)
  }
  
  # Allocate output matrix
  y <- matrix(0L, nrow = n_subjects, ncol = n_time)
  
  # Simulate each subject
  for (s in seq_len(n_subjects)) {
    # Get parameters for this subject's group
    if (homogeneous || n_blocks == 1) {
      marg_s <- marginal
      trans_s <- transition
    } else {
      g <- blocks[s]
      marg_s <- marginal[[g]]
      trans_s <- transition[[g]]
    }
    
    y[s, ] <- .simulate_subject_cat(n_time, p, c, marg_s, trans_s)
  }
  
  y
}


#' Simulate one subject's trajectory
#'
#' @param n_time Number of time points
#' @param p Order
#' @param c Number of categories
#' @param marginal Marginal parameters
#' @param transition Transition parameters
#'
#' @return Integer vector of length n_time
#'
#' @keywords internal
.simulate_subject_cat <- function(n_time, p, c, marginal, transition) {
  y <- integer(n_time)
  
  if (p == 0) {
    # Independence: draw from marginal at each time
    for (k in seq_len(n_time)) {
      probs <- marginal[[k]]
      y[k] <- sample.int(c, size = 1, prob = probs)
    }
    return(y)
  }
  
  # AD(p) with p >= 1
  
  # Draw Y_1 from marginal
  y[1] <- sample.int(c, size = 1, prob = marginal[["t1"]])
  
  # Draw Y_2 to Y_p (if applicable)
  if (p >= 2 && n_time >= 2) {
    for (k in 2:min(p, n_time)) {
      trans_name <- paste0("t", k, "_given_1to", k - 1)
      trans_k <- marginal[[trans_name]]
      
      # Get the row of transition matrix corresponding to history
      # history is y[1:(k-1)]
      idx <- y[1:(k-1)]
      
      if (k == 2) {
        probs <- trans_k[idx[1], ]
      } else {
        # Need to index into multidimensional array
        probs <- trans_k[matrix(c(idx, NA), nrow = 1)[, -k]]
        # Actually, extract the conditional distribution
        probs <- .extract_conditional_probs(trans_k, idx)
      }
      
      y[k] <- sample.int(c, size = 1, prob = probs)
    }
  }
  
  # Draw Y_(p+1) to Y_n
  if (n_time > p) {
    for (k in (p + 1):n_time) {
      trans_name <- paste0("t", k)
      trans_k <- transition[[trans_name]]
      
      # History is y[(k-p):(k-1)]
      history <- y[(k - p):(k - 1)]
      
      # Extract conditional probabilities
      probs <- .extract_conditional_probs(trans_k, history)
      
      y[k] <- sample.int(c, size = 1, prob = probs)
    }
  }
  
  y
}


#' Extract conditional probabilities from transition array
#'
#' @param trans_array Transition array with last dimension being outcome
#' @param history Vector of conditioning values
#'
#' @return Vector of probabilities for current time point
#'
#' @keywords internal
.extract_conditional_probs <- function(trans_array, history) {
  dims <- dim(trans_array)
  n_dims <- length(dims)
  
  if (n_dims == 1) {
    # Marginal distribution (no conditioning)
    return(trans_array)
  }
  
  if (length(history) == 1 && n_dims == 2) {
    # Order 1: matrix
    return(trans_array[history[1], ])
  }
  
  if (length(history) == 2 && n_dims == 3) {
    # Order 2: 3D array
    return(trans_array[history[1], history[2], ])
  }
  
  # General case: build index
  idx <- as.list(history)
  idx[[length(idx) + 1]] <- TRUE  # Select all values of outcome
  
  probs <- do.call(`[`, c(list(trans_array), idx))
  as.numeric(probs)
}


#' Create default marginal parameters (uniform)
#'
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param n_blocks Number of groups
#' @param homogeneous Whether to share parameters
#'
#' @return List of marginal parameters
#'
#' @keywords internal
.default_marginal_cat <- function(p, c, n_time, n_blocks, homogeneous) {
  # Create uniform marginal for one population
  create_one <- function() {
    marg <- list()
    
    if (p == 0) {
      # Independence: uniform marginal at each time point
      for (k in seq_len(n_time)) {
        marg[[paste0("t", k)]] <- rep(1/c, c)
      }
    } else {
      # AD(p): need marginal for t1 and conditionals for t2 to t_p
      marg[["t1"]] <- rep(1/c, c)
      
      if (p >= 2 && n_time >= 2) {
        for (k in 2:min(p, n_time)) {
          # Uniform conditional: each row sums to 1
          dims <- rep(c, k)
          trans_k <- array(1/c, dim = dims)
          marg[[paste0("t", k, "_given_1to", k - 1)]] <- trans_k
        }
      }
    }
    marg
  }
  
  if (homogeneous || n_blocks == 1) {
    return(create_one())
  } else {
    result <- vector("list", n_blocks)
    for (g in seq_len(n_blocks)) {
      result[[g]] <- create_one()
    }
    return(result)
  }
}


#' Create default transition parameters (uniform)
#'
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param n_blocks Number of groups
#' @param homogeneous Whether to share parameters
#'
#' @return List of transition parameters
#'
#' @keywords internal
.default_transition_cat <- function(p, c, n_time, n_blocks, homogeneous) {
  if (p == 0) return(NULL)
  
  # Create uniform transitions for one population
  create_one <- function() {
    trans <- list()
    
    for (k in (p + 1):n_time) {
      # Uniform transition: array of dimension c^(p+1), each "row" sums to 1
      dims <- rep(c, p + 1)
      trans_k <- array(1/c, dim = dims)
      trans[[paste0("t", k)]] <- trans_k
    }
    trans
  }
  
  if (homogeneous || n_blocks == 1) {
    return(create_one())
  } else {
    result <- vector("list", n_blocks)
    for (g in seq_len(n_blocks)) {
      result[[g]] <- create_one()
    }
    return(result)
  }
}


#' Validate transition probability parameters
#'
#' Checks that transition probabilities are valid (non-negative, sum to 1).
#'
#' @param marginal Marginal parameters
#' @param transition Transition parameters
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param tol Tolerance for checking sums
#'
#' @return TRUE if valid, otherwise throws error
#'
#' @keywords internal
.validate_params_cat <- function(marginal, transition, p, c, n_time, tol = 1e-8) {
  # Check marginal
  if (p == 0) {
    for (k in seq_len(n_time)) {
      probs <- marginal[[k]]
      if (any(probs < 0)) stop("Negative probabilities in marginal[[", k, "]]")
      if (abs(sum(probs) - 1) > tol) stop("Probabilities don't sum to 1 in marginal[[", k, "]]")
    }
  } else {
    # Check t1
    if (any(marginal[["t1"]] < 0)) stop("Negative probabilities in marginal$t1")
    if (abs(sum(marginal[["t1"]]) - 1) > tol) stop("Probabilities don't sum to 1 in marginal$t1")
    
    # Check conditionals for t2 to tp
    if (p >= 2) {
      for (k in 2:min(p, n_time)) {
        trans_name <- paste0("t", k, "_given_1to", k - 1)
        trans_k <- marginal[[trans_name]]
        if (any(trans_k < 0)) stop("Negative probabilities in marginal$", trans_name)
      }
    }
  }
  
  # Check transitions
  if (p >= 1 && n_time > p) {
    for (k in (p + 1):n_time) {
      trans_name <- paste0("t", k)
      trans_k <- transition[[trans_name]]
      if (any(trans_k < 0)) stop("Negative probabilities in transition$", trans_name)
    }
  }
  
  TRUE
}
