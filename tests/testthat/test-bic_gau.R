# File: tests/testthat/test-bic_gau.R

test_that("bic_gau matches manual formula for order 1 without blocks", {
    skip_on_cran()
    skip_if_not(exists("simulate_gau"))
    skip_if_not(exists("fit_gau"))

    set.seed(1)
    y <- simulate_gau(
        n_subjects = 30,
        n_time = 5,
        order = 1,
        mu = seq(0, 0.4, length.out = 5),
        sigma = rep(1, 5),
        phi = c(0, rep(0.4, 4))
    )

    fit <- fit_gau(y, order = 1, estimate_mu = TRUE)
    b1 <- bic_gau(fit, n_subjects = nrow(y))
    b2 <- bic_gau(fit)

    N <- ncol(y)
    k <- N + N + (N - 1)

    expect_equal(b1, -2 * fit$log_l + k * log(nrow(y)))
    expect_equal(b2, b1)
})

test_that("bic_gau counts tau and order 2 antedependence parameters", {
    skip_on_cran()
    N <- 6

    fit <- list(
        log_l = -100,
        settings = list(order = 2, estimate_mu = FALSE),
        sigma = rep(1, N),
        mu = rep(0, N),
        phi = matrix(0, nrow = 2, ncol = N),
        tau = c(0, 0.1, 0.2)
    )

    b <- bic_gau(fit, n_subjects = 50)

    k <- 0
    k <- k + N
    k <- k + (2 * N - 3)
    k <- k + (3 - 1)

    expect_equal(b, -2 * fit$log_l + k * log(50))
})
