test_that("dbell probabilities approximately sum to one", {
    theta <- 1
    z_vals <- 0:.BELL_MAX_Z
    p <- dbell(z_vals, theta)
    s <- sum(p)
    expect_true(s > 0.99)
    expect_true(s < 1.01)
})

test_that("pbell is nondecreasing and bounded between zero and one", {
    theta <- 1
    x_vals <- 0:20
    F_vals <- pbell(x_vals, theta)
    expect_true(all(diff(F_vals) >= -1e-12))
    expect_true(all(F_vals >= 0))
    expect_true(all(F_vals <= 1))
})

test_that("rbell returns nonnegative integer values", {
    set.seed(1)
    n <- 50
    theta <- 1
    x <- rbell(n, theta = theta)
    expect_equal(length(x), n)
    expect_true(is.integer(x))
    expect_true(all(x >= 0))
})

test_that("rbell mean is close to theoretical mean for large n", {
    set.seed(2)
    n <- 5000
    theta <- 1
    x <- rbell(n, theta = theta)
    sample_mean <- mean(x)
    theo_mean <- exp(theta)
    expect_equal(sample_mean, theo_mean, tolerance = 0.2)
})

test_that("qbell inverts pbell approximately", {
    theta <- 1
    probs <- c(0.1, 0.5, 0.9)
    q <- qbell(probs, theta = theta, max_z = 80L)
    F_q <- pbell(q, theta = theta)
    expect_true(all(F_q >= probs))
    expect_true(all(F_q - probs < 0.1))
})
