test_that("no fixed effect: logL_inad equals sum of logL_inad_i", {
    skip_on_cran()
    set.seed(1)
    y <- matrix(rpois(60, 2), nrow = 6)
    N <- ncol(y)

    alpha <- rep(0.3, N)
    theta <- rep(2.0, N)

    ll_full <- logL_inad(
        y = y,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = alpha,
        theta = theta,
        blocks = NULL
    )

    ll_sum <- sum(vapply(
        seq_len(N),
        function(i) logL_inad_i(
            y = y,
            i = i,
            order = 1,
            thinning = "binom",
            innovation = "pois",
            alpha = alpha,
            theta = theta
        ),
        numeric(1)
    ))

    expect_equal(ll_full, ll_sum, tolerance = 1e-12)
})

test_that("order 0 reduces to innovation only and matches time sum", {
    skip_on_cran()
    set.seed(2)
    y <- matrix(rpois(30, 3), nrow = 5)
    theta <- rep(1.7, ncol(y))

    ll0 <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.2,
        theta = theta,
        blocks = NULL
    )

    ll0_sum <- sum(vapply(
        seq_len(ncol(y)),
        function(i) logL_inad_i(
            y = y,
            i = i,
            order = 0,
            thinning = "binom",
            innovation = "pois",
            alpha = 0.2,
            theta = theta
        ),
        numeric(1)
    ))

    expect_equal(ll0, ll0_sum, tolerance = 1e-12)
})

test_that("fixed effect: tau length 1 expands and tau[1] is forced to 0", {
    skip_on_cran()
    y <- matrix(0L, nrow = 3, ncol = 4)
    blocks <- c(1L, 2L, 2L)
    theta <- rep(1.0, ncol(y))

    ll_a <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.2,
        theta = theta,
        blocks = blocks,
        tau = 0.4
    )

    ll_b <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.2,
        theta = theta,
        blocks = blocks,
        tau = c(999, 0.4)
    )

    expect_equal(ll_a, ll_b, tolerance = 1e-12)
})

test_that("fixed effect: invalid lambda yields -Inf", {
    skip_on_cran()
    y <- matrix(0L, nrow = 2, ncol = 3)
    blocks <- c(1L, 2L)

    ll <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.1,
        theta = rep(0.5, ncol(y)),
        blocks = blocks,
        tau = -1
    )

    expect_equal(ll, -Inf)
})

test_that("bell FE tau is constrained by innovation mean, not raw Bell theta", {
    skip_on_cran()
    y <- matrix(0L, nrow = 2, ncol = 3)
    blocks <- c(1L, 2L)
    theta <- rep(0.5, ncol(y))

    # Valid under mean-shifted Bell FE: 0.5*exp(0.5) - 0.5 > 0
    ll_ok <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "bell",
        alpha = 0.1,
        theta = theta,
        blocks = blocks,
        tau = -0.5
    )
    expect_true(is.finite(ll_ok))

    # Invalid because Bell innovation mean becomes nonpositive.
    ll_bad <- logL_inad(
        y = y,
        order = 0,
        thinning = "binom",
        innovation = "bell",
        alpha = 0.1,
        theta = theta,
        blocks = blocks,
        tau = -1.0
    )
    expect_equal(ll_bad, -Inf)
})

test_that("innovation nbinom requires nb_inno_size", {
    skip_on_cran()
    set.seed(3)
    y <- matrix(rpois(20, 2), nrow = 4)

    expect_error(
        logL_inad(
            y = y,
            order = 0,
            thinning = "binom",
            innovation = "nbinom",
            alpha = 0.2,
            theta = rep(1.0, ncol(y)),
            blocks = NULL
        )
    )
})

test_that("logL_inad with na_action='fail' errors on missing data", {
    skip_on_cran()
    y <- matrix(c(1L, 2L, NA, 2L, 1L, 0L), nrow = 2, byrow = TRUE)
    alpha <- rep(0.2, ncol(y))
    theta <- rep(1.5, ncol(y))

    expect_error(
        logL_inad(
            y = y,
            order = 1,
            thinning = "binom",
            innovation = "pois",
            alpha = alpha,
            theta = theta,
            na_action = "fail"
        ),
        "contains NA"
    )
})

test_that("logL_inad with na_action='complete' matches complete-case subset", {
    skip_on_cran()
    set.seed(11)
    y_full <- matrix(rpois(80, lambda = 2.0), nrow = 20, ncol = 4)
    y_mis <- y_full
    y_mis[1:6, 4] <- NA

    alpha <- rep(0.3, ncol(y_mis))
    theta <- rep(1.8, ncol(y_mis))

    ll_complete <- logL_inad(
        y = y_mis,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = alpha,
        theta = theta,
        na_action = "complete"
    )

    cc_idx <- complete.cases(y_mis)
    ll_cc <- logL_inad(
        y = y_mis[cc_idx, , drop = FALSE],
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = alpha,
        theta = theta,
        na_action = "fail"
    )

    expect_equal(ll_complete, ll_cc, tolerance = 1e-8)
})

test_that("logL_inad with na_action='marginalize' handles monotone and intermittent missingness", {
    skip_on_cran()
    set.seed(12)
    y <- simulate_inad(
        n_subjects = 40,
        n_time = 5,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = 0.35,
        theta = 2.2
    )

    y_mon <- y
    y_mon[1:8, 4:5] <- NA

    y_int <- y
    y_int[sample(length(y_int), size = 25)] <- NA

    alpha <- rep(0.35, ncol(y))
    theta <- rep(2.2, ncol(y))

    ll_mon <- logL_inad(
        y = y_mon,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = alpha,
        theta = theta,
        na_action = "marginalize"
    )
    ll_int <- logL_inad(
        y = y_int,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = alpha,
        theta = theta,
        na_action = "marginalize"
    )

    expect_true(is.finite(ll_mon))
    expect_true(is.finite(ll_int))
})

test_that(".thin_vec nbinom thinning is degenerate at zero lag count", {
    skip_on_cran()
    k_vals <- 0:6
    probs <- .thin_vec(k_vals, yprev = 0, a = 0.4, thinning = "nbinom")
    expect_equal(probs, c(1, rep(0, length(k_vals) - 1)))
})
