# Internal helpers for CAT hypothesis-test statistics

.array_get <- function(x, idx) {
  if (length(idx) == 0L) {
    return(x)
  }
  do.call(`[`, c(list(x), as.list(idx), list(drop = TRUE)))
}

.cat_pk <- function(k, p) {
  min(k - 1L, p)
}

.cat_prob_under_order <- function(k, p, history, y_cat, marginal, transition) {
  if (k == 1L) {
    return(as.numeric(marginal[["t1"]][y_cat]))
  }

  if (p == 0L) {
    return(as.numeric(marginal[[paste0("t", k)]][y_cat]))
  }

  if (p == 1L) {
    mat_k <- transition[[paste0("t", k)]]
    return(as.numeric(mat_k[history[length(history)], y_cat]))
  }

  if (p == 2L) {
    if (k == 2L) {
      mat_2 <- marginal[["t2_given_1to1"]]
      return(as.numeric(mat_2[history[length(history)], y_cat]))
    }
    arr_k <- transition[[paste0("t", k)]]
    i1 <- history[length(history) - 1L]
    i2 <- history[length(history)]
    return(as.numeric(arr_k[i1, i2, y_cat]))
  }

  stop("Only orders 0, 1, and 2 are supported")
}

.cat_split_populations <- function(y, blocks, fit) {
  y <- as.matrix(y)
  homogeneous <- isTRUE(fit$settings$homogeneous)

  if (homogeneous) {
    return(list(list(
      y = y,
      marginal = fit$marginal,
      transition = fit$transition
    )))
  }

  blocks_eff <- blocks
  if (is.null(blocks_eff)) {
    blocks_eff <- fit$settings$blocks
  }
  if (is.null(blocks_eff)) {
    stop("blocks are required for heterogeneous CAT score statistics")
  }

  block_levels <- fit$settings$block_levels
  if (is.null(block_levels)) {
    block_levels <- as.character(unique(blocks_eff))
  }

  pops <- vector("list", length(block_levels))
  for (i in seq_along(block_levels)) {
    mask_i <- as.character(blocks_eff) == block_levels[[i]]
    pops[[i]] <- list(
      y = y[mask_i, , drop = FALSE],
      marginal = fit$marginal[[i]],
      transition = fit$transition[[i]]
    )
  }
  pops
}

.cat_score_order_single <- function(y, p, q, c, marginal, transition) {
  y <- as.matrix(y)
  n_time <- ncol(y)
  stat <- 0

  # First summand: j = p+1, ..., q
  if (q >= (p + 1L)) {
    for (j in seq.int(p + 1L, min(q, n_time))) {
      if (j == 1L) {
        counts_j <- .count_cells_table_cat(y, 1L, c, subject_mask = NULL)
        N <- nrow(y)
        idx <- which(counts_j > 0)
        for (cat_j in idx) {
          obs <- counts_j[cat_j]
          prob <- .cat_prob_under_order(
            k = 1L, p = p, history = integer(0), y_cat = cat_j,
            marginal = marginal, transition = transition
          )
          expct <- N * prob
          if (expct > 0) {
            stat <- stat + (obs - expct)^2 / expct
          }
        }
      } else {
        counts_j <- .count_cells_table_cat(y, seq_len(j), c, subject_mask = NULL)
        denom_j <- apply(counts_j, seq_len(j - 1L), sum)
        idx <- which(counts_j > 0, arr.ind = TRUE)
        for (r in seq_len(nrow(idx))) {
          id <- idx[r, ]
          history <- id[seq_len(j - 1L)]
          y_cat <- id[j]
          obs <- .array_get(counts_j, id)
          denom <- .array_get(denom_j, history)
          prob <- .cat_prob_under_order(
            k = j, p = p, history = history, y_cat = y_cat,
            marginal = marginal, transition = transition
          )
          expct <- denom * prob
          if (expct > 0) {
            stat <- stat + (obs - expct)^2 / expct
          }
        }
      }
    }
  }

  # Second summand: k = q+1, ..., n
  if (n_time >= (q + 1L)) {
    for (k in seq.int(q + 1L, n_time)) {
      counts_k <- .count_cells_table_cat(y, (k - q):k, c, subject_mask = NULL)
      denom_k <- apply(counts_k, seq_len(q), sum)
      idx <- which(counts_k > 0, arr.ind = TRUE)
      for (r in seq_len(nrow(idx))) {
        id <- idx[r, ]
        history <- id[seq_len(q)]
        y_cat <- id[q + 1L]
        obs <- .array_get(counts_k, id)
        denom <- .array_get(denom_k, history)
        prob <- .cat_prob_under_order(
          k = k, p = p, history = history, y_cat = y_cat,
          marginal = marginal, transition = transition
        )
        expct <- denom * prob
        if (expct > 0) {
          stat <- stat + (obs - expct)^2 / expct
        }
      }
    }
  }

  stat
}

.cat_score_homogeneity <- function(y, blocks, p, c, fit_null) {
  y <- as.matrix(y)
  blocks <- as.integer(blocks)
  n_time <- ncol(y)
  stat <- 0

  marginal_pool <- fit_null$marginal
  transition_pool <- fit_null$transition
  groups <- sort(unique(blocks))

  for (g in groups) {
    y_g <- y[blocks == g, , drop = FALSE]
    N_g <- nrow(y_g)

    for (k in seq_len(n_time)) {
      p_k <- .cat_pk(k, p)

      if (p_k == 0L) {
        counts_k <- .count_cells_table_cat(y_g, k, c, subject_mask = NULL)
        idx <- which(counts_k > 0)
        for (y_cat in idx) {
          obs <- counts_k[y_cat]
          prob <- .cat_prob_under_order(
            k = k, p = p, history = integer(0), y_cat = y_cat,
            marginal = marginal_pool, transition = transition_pool
          )
          expct <- N_g * prob
          if (expct > 0) {
            stat <- stat + (obs - expct)^2 / expct
          }
        }
      } else {
        counts_k <- .count_cells_table_cat(y_g, (k - p_k):k, c, subject_mask = NULL)
        denom_k <- apply(counts_k, seq_len(p_k), sum)
        idx <- which(counts_k > 0, arr.ind = TRUE)
        for (r in seq_len(nrow(idx))) {
          id <- idx[r, ]
          history <- id[seq_len(p_k)]
          y_cat <- id[p_k + 1L]
          obs <- .array_get(counts_k, id)
          denom <- .array_get(denom_k, history)
          prob <- .cat_prob_under_order(
            k = k, p = p, history = history, y_cat = y_cat,
            marginal = marginal_pool, transition = transition_pool
          )
          expct <- denom * prob
          if (expct > 0) {
            stat <- stat + (obs - expct)^2 / expct
          }
        }
      }
    }
  }

  stat
}

.cat_score_timeinvariance_single <- function(y, p, c, pooled_transition) {
  y <- as.matrix(y)
  n_time <- ncol(y)
  stat <- 0

  for (k in seq.int(p + 1L, n_time)) {
    counts_k <- .count_cells_table_cat(y, (k - p):k, c, subject_mask = NULL)
    denom_k <- apply(counts_k, seq_len(p), sum)
    idx <- which(counts_k > 0, arr.ind = TRUE)
    for (r in seq_len(nrow(idx))) {
      id <- idx[r, ]
      history <- id[seq_len(p)]
      y_cat <- id[p + 1L]
      obs <- .array_get(counts_k, id)
      denom <- .array_get(denom_k, history)
      prob <- .array_get(pooled_transition, c(history, y_cat))
      expct <- denom * prob
      if (expct > 0) {
        stat <- stat + (obs - expct)^2 / expct
      }
    }
  }

  stat
}

.cat_score_stationarity_single <- function(y, p, c, stationary_marginal, stationary_transition) {
  y <- as.matrix(y)
  n_time <- ncol(y)
  N <- nrow(y)
  stat <- 0

  if (p == 0L) {
    for (k in seq_len(n_time)) {
      counts_k <- .count_cells_table_cat(y, k, c, subject_mask = NULL)
      idx <- which(counts_k > 0)
      for (y_cat in idx) {
        obs <- counts_k[y_cat]
        prob <- stationary_marginal[y_cat]
        expct <- N * prob
        if (expct > 0) {
          stat <- stat + (obs - expct)^2 / expct
        }
      }
    }
    return(stat)
  }

  if (p %in% c(1L, 2L)) {
    # Under the current stationary null implementation, the first p time points
    # are evaluated using the stationary marginal, and k > p uses a pooled
    # pth-order transition.
    for (k in seq_len(min(p, n_time))) {
      counts_k <- .count_cells_table_cat(y, k, c, subject_mask = NULL)
      idx <- which(counts_k > 0)
      for (y_cat in idx) {
        obs <- counts_k[y_cat]
        prob <- stationary_marginal[y_cat]
        expct <- N * prob
        if (expct > 0) {
          stat <- stat + (obs - expct)^2 / expct
        }
      }
    }

    if (n_time >= (p + 1L)) {
      for (k in seq.int(p + 1L, n_time)) {
        counts_k <- .count_cells_table_cat(y, (k - p):k, c, subject_mask = NULL)
        denom_k <- apply(counts_k, seq_len(p), sum)
        idx <- which(counts_k > 0, arr.ind = TRUE)
        for (r in seq_len(nrow(idx))) {
          id <- idx[r, ]
          history <- id[seq_len(p)]
          y_cat <- id[p + 1L]
          obs <- .array_get(counts_k, id)
          denom <- .array_get(denom_k, history)
          prob <- .array_get(stationary_transition, c(history, y_cat))
          expct <- denom * prob
          if (expct > 0) {
            stat <- stat + (obs - expct)^2 / expct
          }
        }
      }
    }

    return(stat)
  }

  stop(
    "score stationarity test for categorical AD currently supports orders 0, 1, and 2 only",
    call. = FALSE
  )
}

.cat_f_lawley <- function(N, pi) {
  N * pi * log(N * pi) + (1 - pi) / 2 + (1 - pi^2) / (12 * N * pi)
}

.cat_precompute_order_probs <- function(marginal, transition, q, c, n_time) {
  probs <- list(
    marginal = vector("list", n_time),
    pair = vector("list", n_time),
    triple = vector("list", n_time)
  )

  if (q == 0L) {
    for (k in seq_len(n_time)) {
      probs$marginal[[k]] <- as.numeric(marginal[[paste0("t", k)]])
    }
    return(probs)
  }

  probs$marginal[[1L]] <- as.numeric(marginal[["t1"]])

  if (q == 1L) {
    for (k in 2:n_time) {
      t_k <- transition[[paste0("t", k)]]
      pair_k <- diag(probs$marginal[[k - 1L]]) %*% t_k
      probs$pair[[k]] <- pair_k
      probs$marginal[[k]] <- colSums(pair_k)
    }
    return(probs)
  }

  # q == 2
  t2 <- as.matrix(marginal[["t2_given_1to1"]])
  pair_2 <- diag(probs$marginal[[1L]]) %*% t2
  probs$pair[[2L]] <- pair_2
  probs$marginal[[2L]] <- colSums(pair_2)

  if (n_time >= 3L) {
    for (k in 3:n_time) {
      arr_k <- transition[[paste0("t", k)]]
      triple_k <- array(0, dim = c(c, c, c))
      for (a in seq_len(c)) {
        for (b in seq_len(c)) {
          triple_k[a, b, ] <- probs$pair[[k - 1L]][a, b] * arr_k[a, b, ]
        }
      }
      probs$triple[[k]] <- triple_k

      pair_k <- apply(triple_k, c(2, 3), sum)
      probs$pair[[k]] <- pair_k
      probs$marginal[[k]] <- colSums(pair_k)
    }
  }

  probs
}

.cat_hist_prob_from_precomputed <- function(pre, k, ord, history) {
  if (ord == 0L) {
    return(1)
  }
  if (ord == 1L) {
    return(pre$marginal[[k - 1L]][history[1L]])
  }
  if (ord == 2L) {
    return(pre$pair[[k - 1L]][history[1L], history[2L]])
  }
  stop("Only orders 0, 1, and 2 are supported")
}

.cat_event_prob_from_precomputed <- function(pre, k, ord, history, y_cat) {
  if (ord == 0L) {
    return(pre$marginal[[k]][y_cat])
  }
  if (ord == 1L) {
    return(pre$pair[[k]][history[1L], y_cat])
  }
  if (ord == 2L) {
    return(pre$triple[[k]][history[1L], history[2L], y_cat])
  }
  stop("Only orders 0, 1, and 2 are supported")
}

.cat_mlrt_order_e_single <- function(N, p, q, c, n_time, marginal_q, transition_q) {
  if (!isTRUE(q > p)) {
    stop("q must be greater than p for modified order test")
  }

  base <- list(order = p, marginal = vector("list", n_time), pair = vector("list", n_time))

  if (p == 0L) {
    for (k in seq_len(n_time)) {
      base$marginal[[k]] <- as.numeric(marginal_q[[paste0("t", k)]])
    }
  } else if (p == 1L) {
    base$marginal[[1L]] <- as.numeric(marginal_q[["t1"]])
    base$transition <- transition_q
    for (k in 2:n_time) {
      t_k <- as.matrix(transition_q[[paste0("t", k)]])
      pair_k <- diag(base$marginal[[k - 1L]]) %*% t_k
      base$pair[[k]] <- pair_k
      base$marginal[[k]] <- colSums(pair_k)
    }
  } else if (p == 2L) {
    pre2 <- .cat_precompute_order_probs(
      marginal = marginal_q,
      transition = transition_q,
      q = 2L,
      c = c,
      n_time = n_time
    )
    base$marginal <- pre2$marginal
    base$pair <- pre2$pair
    base$triple <- pre2$triple
  } else {
    stop("Only orders 0, 1, and 2 are supported")
  }

  hist_prob_base <- function(k, ord, history) {
    if (ord == 0L) return(1)

    if (ord == 1L) {
      if (p == 0L) {
        return(base$marginal[[k - 1L]][history[1L]])
      }
      return(base$marginal[[k - 1L]][history[1L]])
    }

    if (ord == 2L) {
      if (p == 0L) {
        return(base$marginal[[k - 2L]][history[1L]] * base$marginal[[k - 1L]][history[2L]])
      }
      return(base$pair[[k - 1L]][history[1L], history[2L]])
    }

    stop("Only orders 0, 1, and 2 are supported")
  }

  event_prob_base <- function(k, ord, history, y_cat) {
    if (ord == 0L) {
      return(base$marginal[[k]][y_cat])
    }

    if (ord == 1L) {
      if (p == 0L) {
        return(base$marginal[[k - 1L]][history[1L]] * base$marginal[[k]][y_cat])
      }
      return(base$pair[[k]][history[1L], y_cat])
    }

    if (ord == 2L) {
      if (p == 0L) {
        return(
          base$marginal[[k - 2L]][history[1L]] *
            base$marginal[[k - 1L]][history[2L]] *
            base$marginal[[k]][y_cat]
        )
      }
      if (p == 1L) {
        t_k <- as.matrix(base$transition[[paste0("t", k)]])
        return(base$pair[[k - 1L]][history[1L], history[2L]] * t_k[history[2L], y_cat])
      }
      return(base$triple[[k]][history[1L], history[2L], y_cat])
    }

    stop("Only orders 0, 1, and 2 are supported")
  }

  e_val <- 0
  for (k in 2:n_time) {
    qk <- .cat_pk(k, q)
    pk <- .cat_pk(k, p)
    if (qk <= pk) {
      next
    }

    # First summation in e(p,q): q-order terms, evaluated under AD(p) fit
    term_q <- 0
    if (qk == 1L) {
      for (h1 in seq_len(c)) {
        pi_h <- hist_prob_base(k, qk, c(h1))
        if (pi_h <= 0) next
        for (yc in seq_len(c)) {
          pi_e <- event_prob_base(k, qk, c(h1), yc)
          if (pi_e <= 0) next
          pi_cond <- pi_e / pi_h
          term_q <- term_q + (.cat_f_lawley(N, pi_e) - pi_cond * .cat_f_lawley(N, pi_h))
        }
      }
    } else if (qk == 2L) {
      for (h1 in seq_len(c)) {
        for (h2 in seq_len(c)) {
          pi_h <- hist_prob_base(k, qk, c(h1, h2))
          if (pi_h <= 0) next
          for (yc in seq_len(c)) {
            pi_e <- event_prob_base(k, qk, c(h1, h2), yc)
            if (pi_e <= 0) next
            pi_cond <- pi_e / pi_h
            term_q <- term_q + (.cat_f_lawley(N, pi_e) - pi_cond * .cat_f_lawley(N, pi_h))
          }
        }
      }
    }

    # Second summation in e(p,q): p-order terms
    term_p <- 0
    if (pk == 0L) {
      for (yc in seq_len(c)) {
        pi_e <- event_prob_base(k, 0L, integer(0), yc)
        if (pi_e <= 0) next
        term_p <- term_p + (.cat_f_lawley(N, pi_e) - pi_e * .cat_f_lawley(N, 1))
      }
    } else if (pk == 1L) {
      for (h1 in seq_len(c)) {
        pi_h <- hist_prob_base(k, 1L, c(h1))
        if (pi_h <= 0) next
        for (yc in seq_len(c)) {
          pi_e <- event_prob_base(k, 1L, c(h1), yc)
          if (pi_e <= 0) next
          pi_cond <- pi_e / pi_h
          term_p <- term_p + (.cat_f_lawley(N, pi_e) - pi_cond * .cat_f_lawley(N, pi_h))
        }
      }
    } else {
      for (h1 in seq_len(c)) {
        for (h2 in seq_len(c)) {
          pi_h <- hist_prob_base(k, 2L, c(h1, h2))
          if (pi_h <= 0) next
          for (yc in seq_len(c)) {
            pi_e <- event_prob_base(k, 2L, c(h1, h2), yc)
            if (pi_e <= 0) next
            pi_cond <- pi_e / pi_h
            term_p <- term_p + (.cat_f_lawley(N, pi_e) - pi_cond * .cat_f_lawley(N, pi_h))
          }
        }
      }
    }

    e_val <- e_val + 2 * (term_q - term_p)
  }

  e_val
}

.cat_mlrt_order_e <- function(fit_null, fit_alt) {
  p <- fit_null$settings$order
  q <- fit_alt$settings$order
  c <- fit_null$settings$n_categories
  n_time <- fit_null$settings$n_time

  if (!isTRUE(fit_null$settings$homogeneous)) {
    blocks <- fit_null$settings$blocks
    if (is.null(blocks)) {
      stop("Cannot compute CAT modified LRT for heterogeneous fit without block assignments")
    }
    levels <- fit_null$settings$block_levels
    if (is.null(levels)) {
      levels <- as.character(unique(blocks))
    }

    e_total <- 0
    for (i in seq_along(levels)) {
      N_i <- sum(as.character(blocks) == levels[[i]])
      e_total <- e_total + .cat_mlrt_order_e_single(
        N = N_i,
        p = p,
        q = q,
        c = c,
        n_time = n_time,
        marginal_q = fit_null$marginal[[i]],
        transition_q = fit_null$transition[[i]]
      )
    }
    return(e_total)
  }

  .cat_mlrt_order_e_single(
    N = fit_alt$settings$n_subjects,
    p = p,
    q = q,
    c = c,
    n_time = n_time,
    marginal_q = fit_null$marginal,
    transition_q = fit_null$transition
  )
}

.cat_prob_vector_under_order <- function(k, p, history, c, marginal, transition) {
  probs <- numeric(c)
  for (y_cat in seq_len(c)) {
    probs[y_cat] <- .cat_prob_under_order(
      k = k,
      p = p,
      history = history,
      y_cat = y_cat,
      marginal = marginal,
      transition = transition
    )
  }
  probs
}

.cat_wald_row_multinom <- function(pi_alt, pi_null, N_hist) {
  c <- length(pi_alt)
  if (c <= 1L || N_hist <= 0) {
    return(0)
  }

  p_alt_sub <- as.numeric(pi_alt[seq_len(c - 1L)])
  p_null_sub <- as.numeric(pi_null[seq_len(c - 1L)])
  d <- p_alt_sub - p_null_sub

  if (all(abs(d) <= .Machine$double.eps^0.5)) {
    return(0)
  }

  p_alt_sub <- as.numeric(p_alt_sub)
  V <- diag(p_alt_sub, nrow = length(p_alt_sub), ncol = length(p_alt_sub)) -
    tcrossprod(p_alt_sub)
  ridge <- 1e-10
  V <- V + diag(ridge, nrow = nrow(V))

  w <- tryCatch(
    as.vector(qr.solve(V, d, tol = 1e-10)),
    error = function(e) NULL
  )
  if (is.null(w) || any(!is.finite(w))) {
    return(0)
  }

  as.numeric(N_hist * crossprod(d, w))
}

.cat_wald_order_single <- function(y, p, q, c, marginal_null, transition_null,
                                   marginal_alt, transition_alt) {
  y <- as.matrix(y)
  n_time <- ncol(y)
  stat <- 0

  for (k in 2:n_time) {
    pk <- .cat_pk(k, p)
    qk <- .cat_pk(k, q)
    if (qk <= pk) {
      next
    }

    counts_k <- .count_cells_table_cat(y, (k - qk):k, c, subject_mask = NULL)
    denom_k <- apply(counts_k, seq_len(qk), sum)
    idx_hist <- which(denom_k > 0, arr.ind = TRUE)
    if (is.null(dim(idx_hist))) {
      idx_hist <- matrix(idx_hist, nrow = 1)
    }

    for (r in seq_len(nrow(idx_hist))) {
      history_q <- as.integer(idx_hist[r, seq_len(qk)])
      N_hist <- .array_get(denom_k, history_q)
      if (!is.finite(N_hist) || N_hist <= 0) {
        next
      }

      if (pk == 0L) {
        history_p <- integer(0)
      } else {
        history_p <- history_q[(qk - pk + 1L):qk]
      }

      pi_alt <- .cat_prob_vector_under_order(
        k = k,
        p = q,
        history = history_q,
        c = c,
        marginal = marginal_alt,
        transition = transition_alt
      )
      pi_null <- .cat_prob_vector_under_order(
        k = k,
        p = p,
        history = history_p,
        c = c,
        marginal = marginal_null,
        transition = transition_null
      )

      stat <- stat + .cat_wald_row_multinom(pi_alt, pi_null, N_hist)
    }
  }

  stat
}

.cat_resolve_blocks <- function(blocks, fit) {
  if (!is.null(blocks)) {
    return(as.integer(blocks))
  }
  fit_blocks <- fit$settings$blocks
  if (is.null(fit_blocks)) {
    return(rep.int(1L, fit$settings$n_subjects))
  }
  as.integer(fit_blocks)
}

.cat_default_mlrt_nsim <- function() {
  nsim <- getOption("antedep.cat_mlrt_nsim", 120L)
  nsim <- as.integer(nsim)[1L]
  if (!is.finite(nsim) || nsim < 10L) {
    nsim <- 120L
  }
  nsim
}

.cat_mlrt_expected_lrt_boot <- function(simulate_fn, lrt_raw_fn,
                                        nsim = .cat_default_mlrt_nsim()) {
  vals <- {
    out <- numeric(nsim)
    for (b in seq_len(nsim)) {
      y_b <- simulate_fn()
      out[b] <- lrt_raw_fn(y_b)
    }
    out
  }

  keep <- is.finite(vals) & vals >= 0
  if (!any(keep)) {
    return(NA_real_)
  }
  mean(vals[keep])
}

.cat_timeinvariant_sim_params <- function(fit_null) {
  p <- fit_null$settings$order
  n_time <- fit_null$settings$n_time

  make_one <- function(marginal_one, transition_one) {
    trans <- vector("list", max(0L, n_time - p))
    if (n_time > p) {
      names(trans) <- paste0("t", (p + 1L):n_time)
      for (k in (p + 1L):n_time) {
        trans[[paste0("t", k)]] <- transition_one[["pooled"]]
      }
    } else {
      trans <- NULL
    }
    list(marginal = marginal_one, transition = trans)
  }

  if (isTRUE(fit_null$settings$homogeneous)) {
    return(make_one(fit_null$marginal, fit_null$transition))
  }

  G <- fit_null$settings$n_blocks
  marginal <- vector("list", G)
  transition <- vector("list", G)
  for (g in seq_len(G)) {
    one <- make_one(fit_null$marginal[[g]], fit_null$transition[[g]])
    marginal[[g]] <- one$marginal
    transition[[g]] <- one$transition
  }
  list(marginal = marginal, transition = transition)
}

.cat_stationary_sim_params <- function(fit_null) {
  p <- fit_null$settings$order
  n_time <- fit_null$settings$n_time

  if (p > 1L) {
    stop("modified stationarity test currently supports CAT order 0 or 1 only")
  }

  make_one <- function(marginal_one, transition_one) {
    marg <- list()
    trans <- NULL

    if (p == 0L) {
      pi <- as.numeric(marginal_one[["stationary"]])
      for (k in seq_len(n_time)) {
        marg[[paste0("t", k)]] <- pi
      }
      return(list(marginal = marg, transition = NULL))
    }

    pi <- as.numeric(marginal_one[["stationary"]])
    Q <- transition_one[["stationary"]]
    marg[["t1"]] <- pi
    trans <- vector("list", max(0L, n_time - 1L))
    if (n_time >= 2L) {
      names(trans) <- paste0("t", 2:n_time)
      for (k in 2:n_time) {
        trans[[paste0("t", k)]] <- Q
      }
    } else {
      trans <- NULL
    }
    list(marginal = marg, transition = trans)
  }

  if (isTRUE(fit_null$settings$homogeneous)) {
    return(make_one(fit_null$marginal, fit_null$transition))
  }

  G <- fit_null$settings$n_blocks
  marginal <- vector("list", G)
  transition <- vector("list", G)
  for (g in seq_len(G)) {
    one <- make_one(fit_null$marginal[[g]], fit_null$transition[[g]])
    marginal[[g]] <- one$marginal
    transition[[g]] <- one$transition
  }
  list(marginal = marginal, transition = transition)
}
