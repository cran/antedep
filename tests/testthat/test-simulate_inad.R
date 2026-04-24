test_that("simulate_inad basic structure and defaults", {
    skip_on_cran()
    set.seed(1)
    y <- simulate_inad(
        n_subjects = 5,
        n_time = 7
    )

    expect_true(is.matrix(y))
    expect_equal(dim(y), c(5L, 7L))
    expect_true(is.integer(y))
    expect_true(all(y >= 0))

    y_seed1 <- simulate_inad(n_subjects = 5, n_time = 7, seed = 777)
    y_seed2 <- simulate_inad(n_subjects = 5, n_time = 7, seed = 777)
    expect_equal(y_seed1, y_seed2)
})

test_that("order 0 ignores thinning choice", {
    skip_on_cran()
    set.seed(1)
    y1 <- simulate_inad(
        n_subjects = 10,
        n_time = 6,
        order = 0,
        thinning = "binom",
        innovation = "pois"
    )

    set.seed(1)
    y2 <- simulate_inad(
        n_subjects = 10,
        n_time = 6,
        order = 0,
        thinning = "nbinom",
        innovation = "pois"
    )

    expect_equal(y1, y2)
})

test_that("order 0 matches order 1 with zero alpha", {
    skip_on_cran()
    n_sub <- 8
    n_time <- 6

    set.seed(2)
    y0 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        theta = 5
    )

    zero_alpha <- rep(0, n_time)

    set.seed(2)
    y1 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = zero_alpha,
        theta = 5
    )

    expect_equal(y0, y1)
})

test_that("order 2 runs and produces nonnegative integers", {
    skip_on_cran()
    set.seed(3)
    y <- simulate_inad(
        n_subjects = 4,
        n_time = 8,
        order = 2,
        thinning = "binom",
        innovation = "bell"
    )

    expect_true(is.matrix(y))
    expect_equal(dim(y), c(4L, 8L))
    expect_true(is.integer(y))
    expect_true(all(y >= 0))
})

test_that("negative binomial innovations use size and mu parameterization", {
    skip_on_cran()
    # nb_inno_size is the size (dispersion) parameter (> 0)
    # theta is the mean parameter
    # Mean of NB(size, mu) = mu = theta
    # Variance = mu + mu^2/size = theta * (1 + theta/size)
    
    # With same theta but different size, variance changes
    # Smaller size -> more overdispersion -> higher variance
    set.seed(4)
    y1 <- simulate_inad(
        n_subjects = 100,
        n_time = 5,
        order = 0,  # order 0 so we only test innovation
        innovation = "nbinom",
        theta = 10,
        nb_inno_size = 1  # high overdispersion
    )

    set.seed(5)
    y2 <- simulate_inad(
        n_subjects = 100,
        n_time = 5,
        order = 0,
        innovation = "nbinom",
        theta = 10,
        nb_inno_size = 100  # low overdispersion (closer to Poisson)
    )

    # Both should have mean close to theta = 10
    expect_true(abs(mean(y1) - 10) < 2)
    expect_true(abs(mean(y2) - 10) < 2)
    
    # y1 should have higher variance (smaller size = more overdispersion)
    expect_true(var(as.vector(y1)) > var(as.vector(y2)))
})

test_that("negative binomial innovation mean matches theta", {
    skip_on_cran()
    # Simulation-estimation consistency check
    set.seed(123)
    true_theta <- 8
    true_size <- 2
    
    y <- simulate_inad(
        n_subjects = 500,
        n_time = 1,
        order = 0,
        innovation = "nbinom",
        theta = true_theta,
        nb_inno_size = true_size
    )
    
    # Sample mean should be close to true_theta
    expect_equal(mean(y), true_theta, tolerance = 0.5)
    
    # Theoretical variance = theta * (1 + theta/size) = 8 * (1 + 8/2) = 40
    theo_var <- true_theta * (1 + true_theta / true_size)
    expect_equal(var(as.vector(y)), theo_var, tolerance = 10)
})

test_that("blocks and tau affect innovations as intended", {
    skip_on_cran()
    n_sub <- 200
    n_time <- 3
    blocks <- c(rep(1L, n_sub / 2), rep(2L, n_sub / 2))

    set.seed(10)
    y0 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 0,
        innovation = "pois",
        theta = 10,
        blocks = blocks,
        tau = 3
    )

    m1 <- mean(y0[blocks == 1L, ])
    m2 <- mean(y0[blocks == 2L, ])
    expect_true(m2 > m1)

    set.seed(11)
    y1 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 0,
        innovation = "pois",
        theta = 10,
        blocks = blocks,
        tau = c(100, 3)
    )

    m1b <- mean(y1[blocks == 1L, ])
    m2b <- mean(y1[blocks == 2L, ])
    expect_true(abs(m1b - 10) < 1.5)
    expect_true(m2b > m1b)

    blocks_nonseq <- c(rep(2L, n_sub / 2), rep(5L, n_sub / 2))
    set.seed(12)
    y_nonseq <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 0,
        innovation = "pois",
        theta = 10,
        blocks = blocks_nonseq,
        tau = c(0, 3)
    )
    mn1 <- mean(y_nonseq[blocks_nonseq == 2L, ])
    mn2 <- mean(y_nonseq[blocks_nonseq == 5L, ])
    expect_true(mn2 > mn1)
})

test_that("nbinom thinning handles zero previous counts without introducing NA", {
    skip_on_cran()
    set.seed(99)
    y1 <- simulate_inad(
        n_subjects = 200,
        n_time = 6,
        order = 1,
        thinning = "nbinom",
        innovation = "pois",
        alpha = rep(0.6, 6),
        theta = 1e-12
    )
    expect_true(any(y1[, 1L] == 0L))
    expect_false(anyNA(y1))

    set.seed(100)
    y2 <- simulate_inad(
        n_subjects = 200,
        n_time = 7,
        order = 2,
        thinning = "nbinom",
        innovation = "pois",
        alpha = rbind(c(0, 0.6, rep(0.6, 5)), c(0, 0, rep(0.4, 5))),
        theta = 1e-12
    )
    expect_true(any(y2[, 1L] == 0L))
    expect_false(anyNA(y2))
})

test_that("bell block effects are applied on innovation mean scale", {
    skip_on_cran()
    n_sub <- 4000
    blocks <- c(rep(1L, n_sub / 2), rep(2L, n_sub / 2))

    set.seed(321)
    y <- simulate_inad(
        n_subjects = n_sub,
        n_time = 1,
        order = 0,
        innovation = "bell",
        theta = 0.5,
        blocks = blocks,
        tau = -0.5
    )

    # Base Bell mean at theta=0.5 is 0.5*exp(0.5) ~= 0.824.
    # Group 2 mean should be lower by about 0.5 on the innovation-mean scale.
    m1 <- mean(y[blocks == 1L, 1L])
    m2 <- mean(y[blocks == 2L, 1L])
    expect_gt(m1, m2)
    expect_equal(m1 - m2, 0.5, tolerance = 0.12)
})
