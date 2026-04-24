test_that("bic_inad matches manual formula on bolus_inad", {
    skip_on_cran()
    skip_if_not_installed("nloptr")

    data(bolus_inad, package = "antedep")

    y <- bolus_inad$y
    blocks <- bolus_inad$bolus
    n_subjects <- nrow(y)

    fit <- fit_inad(
        y = y,
        order = 1,
        thinning = "nbinom",
        innovation = "bell",
        blocks = blocks,
        max_iter = 20,
        tol = 1e-6,
        verbose = FALSE
    )

    k_manual <- function(fit) {
        ord <- fit$settings$order
        innovation <- fit$settings$innovation
        N <- length(fit$theta)

        k <- N

        if (ord == 1) k <- k + (N - 1)
        if (ord == 2) k <- k + (2 * N - 3)

        if (!is.null(fit$tau)) {
            B <- length(fit$tau)
            k <- k + (B - 1)
        }

        if (!is.null(fit$nb_inno_size) && innovation == "nbinom") {
            if (length(fit$nb_inno_size) == 1L) k <- k + 1 else k <- k + N
        }

        k
    }

    bic_manual <- -2 * fit$log_l + k_manual(fit) * log(n_subjects)
    bic_pkg <- bic_inad(fit = fit, n_subjects = n_subjects)
    bic_pkg_infer <- bic_inad(fit = fit)

    expect_true(is.numeric(bic_pkg))
    expect_equal(length(bic_pkg), 1)
    expect_true(is.finite(bic_pkg))
    expect_equal(bic_pkg, bic_manual, tolerance = 0)
    expect_equal(bic_pkg_infer, bic_pkg, tolerance = 0)
})

test_that("bic_inad warns on legacy n_subjects fallback via settings$blocks", {
    skip_on_cran()
    fit_legacy <- list(
        log_l = -100,
        theta = c(1, 1, 1),
        tau = NULL,
        nb_inno_size = NULL,
        settings = list(
            order = 1,
            innovation = "pois",
            blocks = rep(1, 25)
        )
    )

    expect_warning(
        bic_legacy <- bic_inad(fit_legacy),
        "legacy fit\\$settings\\$blocks"
    )

    k_manual <- length(fit_legacy$theta) + (length(fit_legacy$theta) - 1)
    bic_manual <- -2 * fit_legacy$log_l + k_manual * log(length(fit_legacy$settings$blocks))
    expect_equal(as.numeric(bic_legacy), as.numeric(bic_manual), tolerance = 1e-12)
})
