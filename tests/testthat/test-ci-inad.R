test_that("ci_inad behaves correctly on bolus_inad", {
    skip_on_cran()
    skip_if_not_installed("nloptr")
    data(bolus_inad, package = "antedep")

    y <- bolus_inad$y
    blocks <- if (!is.null(bolus_inad$blocks)) bolus_inad$blocks else bolus_inad$bolus

    fit <- fit_inad(
        y = y,
        order = 1,
        thinning = "nbinom",
        innovation = "bell",
        blocks = blocks,
        max_iter = 25,
        tol = 1e-6,
        verbose = FALSE
    )

    ci <- ci_inad(
        y = y,
        fit = fit,
        blocks = blocks,
        level = 0.95,
        profile_maxeval = 1200,
        profile_xtol_rel = 1e-6
    )

    expect_s3_class(ci, "inad_ci")

    ## alpha and theta: either NULL or valid data.frame
    if (!is.null(ci$alpha)) {
        expect_true(is.data.frame(ci$alpha))
        expect_true(all(c("param", "est", "se", "lower", "upper", "width", "level") %in% names(ci$alpha)))
        expect_true(all(ci$alpha$lower >= 0))
        expect_true(all(ci$alpha$lower <= ci$alpha$upper))
    }

    if (!is.null(ci$theta)) {
        expect_true(is.data.frame(ci$theta))
        expect_true(all(c("param", "est", "se", "lower", "upper", "width", "level") %in% names(ci$theta)))
        expect_true(all(ci$theta$lower <= ci$theta$upper))
    }

    ## tau: profile CI, may be NULL if blocks degenerate
    if (!is.null(ci$tau)) {
        expect_true(is.data.frame(ci$tau))
        expect_true(all(c("param", "est", "se", "lower", "upper", "width", "level") %in% names(ci$tau)))
        finite_bounds <- is.finite(ci$tau$lower) & is.finite(ci$tau$upper)
        if (any(finite_bounds)) {
            expect_true(all(ci$tau$lower[finite_bounds] <= ci$tau$upper[finite_bounds]))
            expect_true(all(ci$tau$lower[finite_bounds] <= ci$tau$est[finite_bounds]))
            expect_true(all(ci$tau$est[finite_bounds] <= ci$tau$upper[finite_bounds]))
        }
    }
})

test_that("ci_inad never errors for fitted models", {
    skip_on_cran()
    data(bolus_inad, package = "antedep")

    fit <- fit_inad(
        y = bolus_inad$y,
        order = 1,
        thinning = "nbinom",
        innovation = "bell",
        blocks = bolus_inad$bolus
    )

    expect_silent(
        ci_inad(y = bolus_inad$y, fit = fit, blocks = bolus_inad$bolus)
    )
})
