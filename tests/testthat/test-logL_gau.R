testthat::test_that("logL_gau matches manual computation for order 0", {
    set.seed(123)
    n <- 5
    N <- 3
    y <- matrix(rnorm(n * N), nrow = n, ncol = N)
    mu <- c(1, 2, 3)
    sigma <- c(1, 1.5, 2)
    blocks <- c(1, 1, 2, 2, 3)
    tau <- c(0, 0.5, -0.5)

    ll_fun <- logL_gau(
        y = y, order = 0, mu = mu, sigma = sigma,
        blocks = blocks, tau = tau
    )

    ll_manual <- 0
    for (s in seq_len(n)) {
        for (t in seq_len(N)) {
            m <- mu[t] + tau[blocks[s]]
            ll_manual <- ll_manual + dnorm(y[s, t], mean = m, sd = sigma[t], log = TRUE)
        }
    }

    testthat::expect_equal(ll_fun, ll_manual)
})

testthat::test_that("logL_gau matches manual computation for order 1", {
    set.seed(234)
    n <- 5
    N <- 4
    y <- matrix(rnorm(n * N), nrow = n, ncol = N)
    mu <- c(0, 1, 2, 3)
    phi <- c(0, 0.5, 0.6, 0.7)
    sigma <- rep(1, N)
    blocks <- c(1, 1, 2, 2, 2)
    tau <- c(0, 1)

    ll_fun <- logL_gau(
        y = y, order = 1, mu = mu, phi = phi, sigma = sigma,
        blocks = blocks, tau = tau
    )

    ll_manual <- 0
    for (s in seq_len(n)) {
        m1 <- mu[1] + tau[blocks[s]]
        ll_manual <- ll_manual + dnorm(y[s, 1], mean = m1, sd = sigma[1], log = TRUE)
        for (t in 2:N) {
            m_t <- mu[t] + tau[blocks[s]]
            m_tm1 <- mu[t - 1] + tau[blocks[s]]
            cond_mean <- m_t + phi[t] * (y[s, t - 1] - m_tm1)
            ll_manual <- ll_manual + dnorm(y[s, t], mean = cond_mean, sd = sigma[t], log = TRUE)
        }
    }

    testthat::expect_equal(ll_fun, ll_manual)
})

testthat::test_that("logL_gau matches manual computation for order 2", {
    set.seed(345)
    n <- 6
    N <- 5
    y <- matrix(rnorm(n * N), nrow = n, ncol = N)
    mu <- 1:N
    phi <- matrix(0, nrow = 2, ncol = N)
    phi[1, 2:N] <- c(0.3, 0.4, 0.5, 0.6)
    phi[2, 3:N] <- c(0.1, 0.15, 0.2)
    sigma <- rep(1, N)

    ll_fun <- logL_gau(
        y = y, order = 2, mu = mu, phi = phi, sigma = sigma
    )

    ll_manual <- 0
    for (s in seq_len(n)) {
        ll_manual <- ll_manual + dnorm(y[s, 1], mean = mu[1], sd = sigma[1], log = TRUE)
        ll_manual <- ll_manual + dnorm(y[s, 2],
            mean = mu[2] + phi[1, 2] * (y[s, 1] - mu[1]),
            sd = sigma[2], log = TRUE
        )
        for (t in 3:N) {
            cond_mean <- mu[t] +
                phi[1, t] * (y[s, t - 1] - mu[t - 1]) +
                phi[2, t] * (y[s, t - 2] - mu[t - 2])
            ll_manual <- ll_manual + dnorm(y[s, t], mean = cond_mean, sd = sigma[t], log = TRUE)
        }
    }

    testthat::expect_equal(ll_fun, ll_manual)
})

testthat::test_that("logL_gau input validation behavior", {
    y <- matrix(rnorm(20), nrow = 5, ncol = 4)

    testthat::expect_error(logL_gau(y, order = 3))
    testthat::expect_error(logL_gau(y, order = 1, mu = 1:2))
    testthat::expect_error(logL_gau(y, order = 2, phi = matrix(0, nrow = 2, ncol = 2)))
    testthat::expect_error(logL_gau(y, order = 1, blocks = 1))

    testthat::expect_equal(logL_gau(y, order = 0, sigma = -1), -Inf)

    y_bad <- y
    y_bad[1, 1] <- NA
    testthat::expect_error(
        logL_gau(y_bad, order = 0),
        "na_action = 'marginalize' or 'complete'"
    )
    ll_marginalize <- logL_gau(y_bad, order = 0, na_action = "marginalize")
    testthat::expect_true(is.finite(ll_marginalize))
})

testthat::test_that("logL_gau is higher near the generating parameters", {
    set.seed(456)
    n <- 50
    N <- 6
    true_phi <- 0.6
    y <- matrix(0, nrow = n, ncol = N)
    y[, 1] <- rnorm(n)
    for (t in 2:N) {
        y[, t] <- true_phi * y[, t - 1] + rnorm(n)
    }

    ll_true <- logL_gau(y, order = 1, mu = 0, phi = 0.6, sigma = 1)
    ll_wrong <- logL_gau(y, order = 1, mu = 0, phi = 0.1, sigma = 1)

    testthat::expect_true(ll_true > ll_wrong)
})
