test_that("fit_inad works for bolus_inad without fixed effect", {
    skip_on_cran()
    skip_if_not(exists("bolus_inad", where = asNamespace("antedep"), inherits = FALSE) ||
                    exists("bolus_inad", where = parent.env(environment()), inherits = TRUE))

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "bell")

    expect_type(fit, "list")
    expect_true(is.null(fit$tau))
    expect_true(is.matrix(y))
    expect_equal(length(fit$alpha), ncol(y))
    expect_equal(length(fit$theta), ncol(y))

    expect_true(is.finite(fit$log_l))
    expect_true(is.finite(fit$aic))
    expect_true(is.finite(fit$bic))
    expect_true(is.integer(fit$n_params) || is.numeric(fit$n_params))
    expect_equal(fit$n_obs, length(y))
    expect_equal(fit$n_missing, 0)
    expect_true(is.integer(fit$convergence))

    expect_equal(fit$alpha[1], 0, tolerance = 1e-12)
    expect_true(all(is.finite(fit$alpha)))
    expect_true(all(fit$alpha >= -1e-12))
    expect_true(all(fit$alpha < 1))

    expect_true(all(is.finite(fit$theta)))
    expect_true(all(fit$theta > 0))
})

test_that("fit_inad returns log_l equal to sum(loglik_i) in no FE case", {
    skip_on_cran()
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "bell")

    expect_true(!is.null(fit$loglik_i))
    expect_equal(fit$log_l, sum(fit$loglik_i), tolerance = 1e-8)
})

test_that("fit_inad works for order 0 without fixed effect", {
    skip_on_cran()
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 0, thinning = "binom", innovation = "bell")

    expect_true(is.null(fit$alpha) || length(fit$alpha) == 0)
    expect_equal(length(fit$theta), ncol(y))
    expect_true(is.finite(fit$log_l))
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(fit$theta > 0))
})

test_that("fit_inad handles order 2 structural zeros without fixed effect", {
    skip_on_cran()
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 2, thinning = "binom", innovation = "bell")

    expect_true(is.matrix(fit$alpha))
    expect_equal(dim(fit$alpha), c(ncol(y), 2))

    expect_equal(fit$alpha[1, 1], 0, tolerance = 1e-12)
    expect_equal(fit$alpha[1, 2], 0, tolerance = 1e-12)
    expect_equal(fit$alpha[2, 2], 0, tolerance = 1e-12)

    expect_true(all(is.finite(fit$alpha)))
    expect_true(all(fit$alpha >= -1e-12))
    expect_true(all(fit$alpha < 1))

    expect_true(is.finite(fit$log_l))
})

test_that("fit_inad works with fixed effect and tau normalization", {
    skip_on_cran()
    skip_if_not_installed("nloptr")

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y
    blocks <- bolus_inad$blocks

    fit1 <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "bell",
        blocks = blocks,
        init_tau = 0.4,
        max_iter = 5
    )

    expect_true(is.finite(fit1$log_l))
    expect_true(!is.null(fit1$tau))
    expect_equal(fit1$tau[1], 0, tolerance = 1e-12)
    expect_equal(sort(fit1$settings$block_levels), c("1", "2"))

    B <- max(as.integer(blocks))
    expect_equal(length(fit1$tau), B)

    fit2 <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "bell",
        blocks = blocks,
        init_tau = c(999, rep(0.2, B - 1)),
        max_iter = 5
    )

    expect_equal(fit2$tau[1], 0, tolerance = 1e-12)
})

test_that("fit_inad preserves original block labels in settings", {
    skip_on_cran()
    skip_if_not_installed("nloptr")

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y
    blocks_labeled <- ifelse(bolus_inad$blocks == 1, 5, 2)

    fit <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "bell",
        blocks = blocks_labeled,
        max_iter = 5
    )

    expect_equal(sort(fit$settings$block_levels), c("2", "5"))
    expect_equal(length(fit$tau), 2)
})

test_that("fit_inad works for innovation nbinom and returns nb_inno_size", {
    skip_on_cran()
    skip_if_not_installed("nloptr")

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "nbinom",
        init_nb_inno_size = 1
    )

    expect_true(!is.null(fit$nb_inno_size))
    expect_equal(length(fit$nb_inno_size), ncol(y))
    expect_true(all(is.finite(fit$nb_inno_size)))
    expect_true(all(fit$nb_inno_size > 0))
    expect_true(is.finite(fit$log_l))
})

test_that("fit_inad respects nb_inno_size_ub upper bound", {
    skip_on_cran()
    set.seed(42)
    y <- matrix(rpois(60, lambda = 3), nrow = 20, ncol = 3)

    ub <- 2
    fit <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "nbinom",
        nb_inno_size_ub = ub,
        init_nb_inno_size = 1,
        max_iter = 10
    )

    expect_true(all(fit$nb_inno_size <= ub + 1e-8))
    expect_true(is.finite(fit$log_l))
})

test_that("fit_inad validates nb_inno_size_ub", {
    skip_on_cran()
    y <- matrix(rpois(40, 2), nrow = 10, ncol = 4)

    expect_error(
        fit_inad(y, order = 1, thinning = "binom", innovation = "nbinom",
                 nb_inno_size_ub = -1),
        "nb_inno_size_ub must be a positive finite scalar"
    )
    expect_error(
        fit_inad(y, order = 1, thinning = "binom", innovation = "nbinom",
                 nb_inno_size_ub = Inf),
        "nb_inno_size_ub must be a positive finite scalar"
    )
})

test_that("fit_inad with na_action='fail' errors on missing data", {
    skip_on_cran()
    y <- matrix(c(1L, 2L, NA, 3L, 2L, 1L), nrow = 2, byrow = TRUE)

    expect_error(
        fit_inad(
            y,
            order = 1,
            thinning = "binom",
            innovation = "pois",
            na_action = "fail"
        ),
        "contains NA"
    )
})

test_that("fit_inad with na_action='complete' runs on complete-case subset", {
    skip_on_cran()
    set.seed(101)
    y <- matrix(rpois(120, lambda = 2), nrow = 30, ncol = 4)
    y[1:8, 4] <- NA

    fit <- suppressWarnings(
        fit_inad(
            y,
            order = 1,
            thinning = "binom",
            innovation = "pois",
            na_action = "complete"
        )
    )

    expect_true(is.finite(fit$log_l))
    expect_true(all(is.finite(fit$alpha)))
    expect_true(all(is.finite(fit$theta)))
})

test_that("fit_inad marginalize handles monotone MAR missingness", {
    skip_on_cran()

    set.seed(202)
    y_complete <- simulate_inad(
        n_subjects = 50,
        n_time = 5,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.35,
        theta = 2.0
    )

    # MAR monotone dropout probability depends on observed baseline count.
    y_mis <- y_complete
    p_drop <- stats::plogis(-1 + 0.2 * y_complete[, 1])
    drop <- runif(nrow(y_complete)) < p_drop
    y_mis[drop, 4:5] <- NA

    fit_marg <- fit_inad(
        y_mis,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        na_action = "marginalize",
        max_iter = 12,
        tol = 1e-5
    )
    fit_full <- fit_inad(
        y_complete,
        order = 1,
        thinning = "binom",
        innovation = "pois"
    )

    expect_true(is.finite(fit_marg$log_l))
    expect_gt(fit_marg$n_missing, 0)
    expect_true(is.finite(fit_marg$aic))
    expect_true(is.finite(fit_marg$bic))
    expect_true(is.integer(fit_marg$convergence))
    expect_true(all(is.finite(fit_marg$alpha)))
    expect_true(all(is.finite(fit_marg$theta)))
    expect_equal(fit_marg$theta, fit_full$theta, tolerance = 1.2)
})

test_that("fit_inad marginalize handles intermittent MAR missingness", {
    skip_on_cran()

    set.seed(303)
    y_complete <- simulate_inad(
        n_subjects = 40,
        n_time = 6,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.30,
        theta = 1.8
    )

    # Intermittent MAR: missingness depends on observed baseline and time.
    y_mis <- y_complete
    for (t in 2:ncol(y_mis)) {
        p_mis_t <- stats::plogis(-2 + 0.18 * y_complete[, 1] + 0.15 * t)
        miss_t <- runif(nrow(y_mis)) < p_mis_t
        y_mis[miss_t, t] <- NA
    }

    fit_marg <- fit_inad(
        y_mis,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        na_action = "marginalize",
        max_iter = 12,
        tol = 1e-5
    )

    expect_true(is.finite(fit_marg$log_l))
    expect_gt(fit_marg$n_missing, 0)
    expect_true(is.finite(fit_marg$aic))
    expect_true(is.finite(fit_marg$bic))
    expect_true(is.integer(fit_marg$convergence))
    expect_true(all(is.finite(fit_marg$alpha)))
    expect_true(all(is.finite(fit_marg$theta)))
    expect_true(!is.null(fit_marg$settings$na_action))
    expect_equal(fit_marg$settings$na_action, "marginalize")
})
