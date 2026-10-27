test_that("chi-bar-square p-value treats chi^2_0 as a point mass at zero", {
    r <- 10
    stat <- 25.0687
    expected <- sum(choose(r, 1:r) * 2^-r * stats::pchisq(stat, 1:r, lower.tail = FALSE))
    expect_equal(antedep:::.chi_bar_pvalue(stat, r), expected, tolerance = 1e-12)
    expect_equal(antedep:::.chi_bar_pvalue(0, r), 1)
    expect_equal(antedep:::.chi_bar_pvalue(1, 0), 0)
})

test_that("parameter counts agree between fits and model tables", {
    expect_equal(antedep:::.count_params_inad(1, 12, 2, "nbinom", "nbinom"), 11 + 12 + 1 + 12)
    expect_equal(antedep:::.count_params_inad(2, 12, 2, "nbinom", "nbinom"), 21 + 12 + 1 + 12)
    cinfo <- antedep:::.parse_constraint(c("theta", "kappa"), 1, 12, "nbinom")
    expect_equal(cinfo$df, 22)
    expect_equal(antedep:::.parse_constraint("all", 1, 12, "nbinom")$df, 32)
    expect_equal(antedep:::.parse_constraint("alpha1", 2, 12, "pois")$df, 10)
    expect_error(antedep:::.parse_constraint("kappa", 1, 12, "pois"), "nbinom")
})

test_that("Bell parameter inversion is accurate to machine precision", {
    mu <- c(0.05, 0.7, 3, 25)
    th <- antedep:::.bell_theta_from_mean(mu)
    expect_lt(max(abs(th * exp(th) - mu) / mu), 1e-12)
})

test_that("Louis information equals the numerical Hessian, with blocks and constraints", {
    skip_on_cran()
    skip_if_not_installed("numDeriv")
    set.seed(42)
    n_time <- 4
    cases <- list(
        list(order = 1, thinning = "nbinom", innovation = "nbinom", constrain = "alpha"),
        list(order = 1, thinning = "pois", innovation = "bell", constrain = NULL),
        list(order = 2, thinning = "binom", innovation = "nbinom", constrain = NULL)
    )
    for (cs in cases) {
        alpha <- if (cs$order == 1) c(0, 0.3, 0.4, 0.35) else
            rbind(c(0, 0.3, 0.3, 0.35), c(0, 0, 0.2, 0.15))
        theta <- if (cs$innovation == "bell") rep(0.8, n_time) else c(3, 2, 2.5, 2)
        blocks <- rep(1:2, each = 30)
        y <- simulate_inad(60, n_time, order = cs$order, thinning = cs$thinning,
                           innovation = cs$innovation, alpha = alpha, theta = theta,
                           nb_inno_size = 2, blocks = blocks, tau = 0.6)
        fit <- fit_inad(y, cs$order, cs$thinning, cs$innovation, blocks = blocks,
                        max_iter = 50, tol = 1e-8, constrain = cs$constrain)
        li <- antedep:::.louis_full_info(y, fit, blocks)
        ni <- antedep:::.inad_numeric_info(y, fit, blocks, step = 1e-3)
        # Numerical differentiation on the log scale is unreliable for estimates
        # very close to zero, so compare the parameters away from zero.
        ok <- li$fp$tab$est > 0.01 | li$fp$tab$group == "tau"
        expect_lt(max(abs(li$info[ok, ok] - ni$info[ok, ok])) / max(abs(ni$info[ok, ok])), 1e-3)
        ci <- ci_inad(y, fit, blocks = blocks)
        expect_s3_class(ci, "inad_ci")
        expect_true(all(ci$theta$lower < ci$theta$est & ci$theta$est < ci$theta$upper))
    }
})

test_that("constrained fits re-estimate NB innovation sizes and reach the MLE", {
    skip_on_cran()
    set.seed(7)
    blocks <- rep(1:2, each = 40)
    y <- simulate_inad(80, 5, order = 1, thinning = "nbinom", innovation = "nbinom",
                       alpha = c(0, 0.4, 0.4, 0.4, 0.4), theta = c(4, 2, 3, 2, 2.5),
                       nb_inno_size = c(1, 3, 0.8, 2, 1.5), blocks = blocks, tau = 0.8)
    fit_u <- fit_inad(y, 1, "nbinom", "nbinom", blocks = blocks, max_iter = 100, tol = 1e-8)
    fit_c <- fit_inad(y, 1, "nbinom", "nbinom", blocks = blocks, max_iter = 100, tol = 1e-8,
                      constrain = "alpha")
    expect_equal(fit_c$n_params, fit_u$n_params - 3)
    expect_false(isTRUE(all.equal(fit_c$nb_inno_size, fit_u$nb_inno_size)))
    expect_lte(fit_c$log_l, fit_u$log_l + 1e-6)
    li <- antedep:::.louis_full_info(y, fit_c, blocks)
    expect_lt(max(abs(li$score)), 1e-5)
    st <- test_stationarity_inad(y, 1, "nbinom", "nbinom", blocks = blocks,
                                 constrain = "innovation", fit_unconstrained = fit_u)
    expect_equal(st$df, 8)
    expect_equal(st$table$n_params[2], fit_u$n_params - 8)
})

test_that("EM is exact for negative binomial innovations and matches direct maximization", {
    skip_on_cran()
    set.seed(3)
    blocks <- rep(1:2, each = 40)
    y <- simulate_inad(80, 4, order = 1, thinning = "pois", innovation = "nbinom",
                       alpha = c(0, 0.5, 0.4, 0.6), theta = c(3, 2, 2, 2.5),
                       nb_inno_size = 1.5, blocks = blocks, tau = 0.5)
    em <- em_inad(y, 1, "pois", "nbinom", blocks = blocks, max_iter = 5000, tol = 1e-10)
    direct <- fit_inad(y, 1, "pois", "nbinom", blocks = blocks, max_iter = 200, tol = 1e-10)
    expect_true(all(diff(em$loglik_trace) > -1e-8))
    expect_equal(em$log_l, direct$log_l, tolerance = 1e-4)
    expect_equal(em$log_l,
                 logL_inad(y, 1, "pois", "nbinom", em$alpha, em$theta, em$nb_inno_size,
                           blocks, em$tau),
                 tolerance = 1e-10)
    expect_equal(length(em$nb_inno_size), 4)
})

test_that("direct fitting and EM use the same parameter space when alpha exceeds one", {
    skip_on_cran()
    # PT-PI conditional mean alpha * x + theta; alpha = 2, theta = 1 fits exactly.
    x <- rep(0:10, 3)
    y <- cbind(x, 2 * x + 1)
    truth <- logL_inad(y, 1, "pois", "pois", alpha = c(0, 2), theta = c(5, 1))
    direct <- fit_inad(y, order = 1, thinning = "pois", innovation = "pois")
    em <- em_inad(y, order = 1, thinning = "pois", innovation = "pois", max_iter = 10000,
                  tol = 1e-11, newton_polish = TRUE)
    expect_equal(direct$log_l, truth, tolerance = 1e-8)
    expect_equal(em$log_l, truth, tolerance = 1e-8)
    expect_equal(direct$alpha[2], 2, tolerance = 1e-5)

    set.seed(3)
    yy <- simulate_inad(80, 4, order = 1, thinning = "pois", innovation = "pois",
                        alpha = c(0, 1.6, 1.3, 0.5), theta = c(2, 1, 1, 2))
    bl <- rep(1:2, each = 40)
    d <- fit_inad(yy, 1, "pois", "pois", blocks = bl, max_iter = 200, tol = 1e-9)
    e <- em_inad(yy, 1, "pois", "pois", blocks = bl, max_iter = 20000, tol = 1e-11,
                 newton_polish = TRUE)
    expect_gt(max(d$alpha), 1)
    expect_equal(d$log_l, e$log_l, tolerance = 1e-6)
    expect_true(kkt_inad(d)$kkt_ok)
    fc <- fit_inad(yy, 1, "nbinom", "pois", blocks = bl, max_iter = 200, tol = 1e-9,
                   constrain = "alpha")
    expect_true(kkt_inad(fc)$kkt_ok)
})

test_that("binomial thinning keeps alpha in [0, 1] and zero predecessors give alpha = 0", {
    skip_on_cran()
    set.seed(11)
    y <- simulate_inad(50, 4, order = 1, thinning = "binom", innovation = "pois",
                       alpha = 0.5, theta = 2)
    y[, 2] <- 0
    d <- fit_inad(y, 1, "binom", "pois")
    e <- em_inad(y, 1, "binom", "pois", max_iter = 5000, tol = 1e-10)
    expect_true(all(d$alpha >= 0 & d$alpha <= 1))
    expect_equal(d$alpha[3], 0)
    expect_equal(e$alpha[3], 0)
    expect_equal(d$log_l, e$log_l, tolerance = 1e-5)
    k <- kkt_inad(e)
    expect_equal(k$n_unidentified, 1)
})

test_that("order tests report no uncalibrated p-value by default", {
    skip_on_cran()
    set.seed(5)
    y <- simulate_inad(40, 4, order = 1, thinning = "binom", innovation = "pois",
                       alpha = 0.4, theta = 2)
    bl <- rep(1:2, each = 20)
    out <- test_order_inad(y, 1, 2, "binom", "pois")
    expect_true(is.na(out$p.value))
    expect_true(is.na(out$p_value))
    expect_identical(out$reference, "none")
    ch <- test_order_inad(y, 1, 2, "binom", "pois", reference = "chibar",
                          fit_null = out$fit_null, fit_alt = out$fit_alt)
    expect_equal(ch$p_value, antedep:::.chi_bar_pvalue(out$lrt_stat, 2), tolerance = 1e-12)
    expect_warning(test_order_inad(y, 0, 1, "binom", "pois", blocks = bl, reference = "chibar",
                                   max_iter = 20),
                   "not been established")
    f1 <- fit_inad(y, 1, "binom", "pois")
    f2 <- fit_inad(y, 2, "binom", "pois")
    av <- anova(f1, f2)
    expect_true(is.na(av[["Pr(>Chisq)"]][2]))
})

test_that("direct order-2 fits reach the maximum when an innovation mean is on its boundary", {
    skip_on_cran()
    a1 <- c(0, 0.6, 0.4, 0.6, 0.4, 0.6)
    a2 <- c(0, 0, 0.2, 0.3, 0.2, 0.3)
    set.seed(1)
    y <- simulate_inad(30, 6, order = 2, thinning = "binom", innovation = "pois",
                       alpha = rbind(a1, a2), theta = 5 * (1 - a1 - a2))
    d <- fit_inad(y, 2, "binom", "pois")
    e <- em_inad(y, 2, "binom", "pois", max_iter = 20000, tol = 1e-10, newton_polish = TRUE)
    expect_equal(d$log_l, e$log_l, tolerance = 1e-7)
    k <- kkt_inad(d)
    expect_true(k$kkt_ok)
    expect_true(any(k$table$status == "lower bound" & grepl("^theta", k$table$param)))
    f1 <- fit_inad(y, 1, "binom", "pois")
    expect_gte(d$log_l, f1$log_l - 1e-8)
    fc <- fit_inad(y, 2, "binom", "pois", constrain = "alpha")
    expect_true(kkt_inad(fc)$kkt_ok)
    expect_lte(fc$log_l, d$log_l + 1e-8)
})

test_that("missing-data likelihood sums over thinned counts above the previous count", {
    skip_on_cran()
    # One missing value at occasion 2 of an INAD(1) series with Poisson and
    # negative binomial thinning: compare with direct summation over it.
    y <- matrix(c(2, NA, 5, 1, 3, 2), nrow = 2, byrow = TRUE)
    for (th in c("pois", "nbinom")) {
        al <- c(0, 1.5, 0.8); tt <- c(2, 1, 1.5)
        marg <- logL_inad(y, 1, th, "pois", alpha = al, theta = tt, na_action = "marginalize")
        brute <- log(sum(vapply(0:300, function(m) {
            exp(logL_inad(rbind(c(2, m, 5)), 1, th, "pois", alpha = al, theta = tt))
        }, numeric(1)))) + logL_inad(y[2, , drop = FALSE], 1, th, "pois", alpha = al, theta = tt)
        expect_equal(marg, brute, tolerance = 1e-8)
    }
})

test_that("Newton polish moves a coefficient just above zero onto the bound", {
    skip_on_cran()
    # First-order data: the second-order fit has its second-lag coefficients at 0.
    set.seed(5)
    y <- simulate_inad(40, 4, order = 1, thinning = "binom", innovation = "pois",
                       alpha = 0.4, theta = 2)
    d <- fit_inad(y, 2, "binom", "pois")
    expect_equal(d$alpha[3:4, 2], c(0, 0))
    expect_true(kkt_inad(d)$kkt_ok)
    # Perturb a boundary coefficient slightly into the interior, then polish.
    bumped <- d
    bumped$alpha[4, 2] <- 2e-8
    pol <- antedep:::.inad_newton_polish(y, bumped, nb_inno_size_ub = 50)
    expect_equal(unname(pol$fit$alpha[4, 2]), 0)
    expect_equal(pol$log_l, d$log_l, tolerance = 1e-9)
})

test_that("constrained fits reach the maximum by damped Newton steps from pooled values", {
    skip_on_cran()
    # A replicate of the time-varying second-order design (n = 12, N = 30) on which
    # undamped Newton steps stalled and the derivative-free search took minutes.
    a1 <- c(0, rep(c(0.6, 0.4), length.out = 11))
    a2 <- c(0, 0, rep(c(0.2, 0.3), length.out = 10))
    set.seed(20261024 + 100000 * 42 + 95)
    y <- simulate_inad(30, 12, order = 2, thinning = "binom", innovation = "pois",
                       alpha = rbind(a1, a2), theta = 5 * (1 - a1 - a2))
    fc <- fit_inad(y, 2, "binom", "pois", constrain = "alpha")
    expect_identical(fc$settings$path, "newton")
    expect_true(kkt_inad(fc)$kkt_ok)
    expect_equal(fc$log_l, -722.521188, tolerance = 1e-6)
})

test_that("the score at a degenerate bound is the one-sided derivative of the log-likelihood", {
    skip_on_cran()
    a1 <- c(0, rep(0.4, 4)); a2 <- c(0, 0, rep(0.2, 3))
    for (th in c("binom", "pois", "nbinom")) {
        set.seed(11)
        y <- simulate_inad(30, 5, order = 2, thinning = th, innovation = "pois",
                           alpha = rbind(a1, a2), theta = 2)
        fit <- fit_inad(y, 2, th, "pois")
        for (lag in 1:2) {
            f0 <- fit
            f0$alpha[4, lag] <- 0
            ll <- function(v) {
                A <- f0$alpha; A[4, lag] <- v
                logL_inad(y, 2, th, "pois", alpha = A, theta = f0$theta)
            }
            h <- 1e-6
            fd <- (-3 * ll(0) + 4 * ll(h) - ll(2 * h)) / (2 * h)
            s <- antedep:::.louis_full_info(y, f0)$score[[sprintf("alpha%d[4]", lag)]]
            expect_equal(s, fd, tolerance = 1e-6)
        }
    }
    # Binomial thinning at alpha = 1: the left derivative.
    set.seed(12)
    y <- simulate_inad(30, 4, order = 1, thinning = "binom", innovation = "pois",
                       alpha = 0.6, theta = 2)
    y[, 3] <- pmax(y[, 3], y[, 2])
    fit <- fit_inad(y, 1, "binom", "pois")
    fit$alpha[3] <- 1
    ll <- function(v) { a <- fit$alpha; a[3] <- v; logL_inad(y, 1, "binom", "pois", alpha = a, theta = fit$theta) }
    h <- 1e-6
    fd <- -(-3 * ll(1) + 4 * ll(1 - h) - ll(1 - 2 * h)) / (2 * h)
    expect_equal(antedep:::.louis_full_info(y, fit)$score[["alpha[3]"]], fd, tolerance = 1e-6)
})

test_that("the thinning pmf and the score stay accurate for a coefficient just above zero", {
    skip_on_cran()
    z <- 0:6
    for (a in c(1e-12, 1e-6, 0.3, 2)) {
        exact <- exp(lchoose(7 + z - 1, z) + z * log(a) - (7 + z) * log1p(a))
        # Each probability to relative accuracy, however small it is.
        expect_equal(antedep:::.thin_vec(z, 7, a, "nbinom") / exact, rep(1, length(z)),
                     tolerance = 1e-12)
    }
    # EM leaves a coefficient at its lower limit 1e-12. The score there, from
    # Fisher's identity, must agree with the one-sided derivative at zero.
    a1 <- c(0, rep(0.4, 4)); a2 <- c(0, 0, rep(0.2, 3))
    for (th in c("binom", "pois", "nbinom")) {
        set.seed(11)
        y <- simulate_inad(30, 5, order = 2, thinning = th, innovation = "pois",
                           alpha = rbind(a1, a2), theta = 2)
        fit <- fit_inad(y, 2, th, "pois")
        for (lag in 1:2) {
            f0 <- fit; f0$alpha[4, lag] <- 0
            f1 <- fit; f1$alpha[4, lag] <- 1e-12
            nm <- sprintf("alpha%d[4]", lag)
            expect_equal(antedep:::.louis_full_info(y, f1)$score[[nm]],
                         antedep:::.louis_full_info(y, f0)$score[[nm]], tolerance = 1e-8)
        }
    }
})

test_that("kkt_inad rejects a zero coefficient whose one-sided score is positive, and polish moves it", {
    x <- rep(0:10, 3)
    y <- cbind(x, 2 * x + 1)
    f <- fit_inad(y, order = 1, thinning = "pois", innovation = "pois")
    f$alpha <- c(0, 0)
    f$theta <- colMeans(y)
    k <- kkt_inad(f, y)
    # Right derivative sum_s x_s (y_s2 / 11 - 1) = 60.
    expect_equal(k$table$score[k$table$param == "alpha[2]"], 60, tolerance = 1e-8)
    expect_false(k$kkt_ok)
    pol <- antedep:::.inad_newton_polish(y, f, nb_inno_size_ub = 50)
    expect_gt(pol$fit$alpha[2], 0.01)
    expect_gt(pol$log_l, logL_inad(y, 1, "pois", "pois", alpha = c(0, 0), theta = f$theta) + 1)
})

test_that("kkt_inad does not accept a large score when the information is indefinite", {
    y <- matrix(c(4, 1, 5, 9, 4, 9, 4, 5, 1, 4, 5, 3, 5, 9, 7, 6, 5, 5, 8, 4, 5, 1, 2, 5, 5,
                  4, 4, 5, 7, 5, 4, 3, 6, 8, 6, 5, 6, 7, 7, 2, 8, 3, 3, 7, 4, 4, 3, 8, 6, 3,
                  5, 3, 2, 2, 7, 8, 6, 5, 5, 3, 5, 6, 4, 2, 8, 2, 5, 3, 4, 5, 5, 6, 4, 5, 7),
                ncol = 3)
    f <- fit_inad(y, 2, "binom", "nbinom")
    f$alpha <- matrix(c(0, 0.609216191619635, 0.237506282795221, 0, 0, 0.434348776517436), 3)
    f$theta <- c(1.36619792720303, 3.24480440441985, 0.788112287700642)
    f$nb_inno_size <- c(3.16484854647424, 3.41801985270577, 6.88928618320497)
    k <- kkt_inad(f, y)
    expect_false(k$info_pd)
    expect_false(k$kkt_ok)
    expect_gt(k$max_abs_interior_score, 40)
    expect_false(all(k$table$ok[k$table$status == "interior"]))
})
