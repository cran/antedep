test_that("standard methods work for Gaussian AD fits", {
  set.seed(11)
  y <- simulate_gau(n_subjects = 30, n_time = 4, order = 1, phi = 0.35)
  fit0 <- fit_gau(y, order = 0)
  fit1 <- fit_gau(y, order = 1)

  expect_s3_class(fit1, "ad_fit")
  expect_type(coef(fit1), "double")
  expect_true(length(coef(fit1)) > 0)
  expect_type(coef(fit1, type = "list"), "list")

  expect_s3_class(logLik(fit1), "logLik")
  expect_equal(nobs(fit1), nrow(y))
  expect_s3_class(BIC(fit0, fit1), "data.frame")

  ci <- confint(fit1)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("lower", "upper"))

  V <- vcov(fit1)
  expect_true(is.matrix(V))
  expect_equal(rownames(V), names(coef(fit1)))
  expect_equal(colnames(V), names(coef(fit1)))

  tab <- anova(fit0, fit1)
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 2)

  skip_if_not_installed("lmtest")
  expect_s3_class(lmtest::coeftest(fit1), "coeftest")
})

test_that("standard methods work for categorical AD fits", {
  set.seed(12)
  y <- simulate_cat(n_subjects = 40, n_time = 4, order = 1, n_categories = 2)
  fit <- fit_cat(y, order = 1)

  expect_s3_class(fit, "ad_fit")
  expect_type(coef(fit), "double")
  expect_type(coef(fit, type = "list"), "list")

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(colnames(ci), c("lower", "upper"))

  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_equal(rownames(V), names(coef(fit)))

  skip_if_not_installed("lmtest")
  expect_s3_class(lmtest::coeftest(fit), "coeftest")
})

test_that("categorical helper objects print without table access", {
  set.seed(121)
  y <- simulate_cat(n_subjects = 40, n_time = 4, order = 1, n_categories = 2)

  ord <- bic_order_cat(y, max_order = 1)
  expect_s3_class(ord, "bic_order_cat")
  ord_out <- capture.output(print(ord))
  expect_true(any(grepl("Best order", ord_out)))

  tests <- run_order_tests_cat(y, max_order = 1)
  expect_s3_class(tests, "cat_order_tests")
  tests_out <- capture.output(print(tests))
  expect_true(any(grepl("Selected order", tests_out)))

  stat <- run_stationarity_tests_cat(y, order = 1)
  expect_s3_class(stat, "stationarity_tests_cat")
  stat_out <- capture.output(print(stat))
  expect_true(any(grepl("Stationarity", stat_out)))
})

test_that("standard methods work for INAD fits where data-free methods apply", {
  set.seed(13)
  y <- simulate_inad(
    n_subjects = 12,
    n_time = 3,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    alpha = 0.3,
    theta = 2
  )
  fit0 <- fit_inad(y, order = 0, thinning = "binom", innovation = "pois", max_iter = 5)
  fit1 <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois", max_iter = 5)

  expect_s3_class(fit1, "ad_fit")
  expect_type(coef(fit1), "double")
  expect_type(coef(fit1, type = "list"), "list")
  expect_s3_class(logLik(fit1), "logLik")

  tab <- anova(fit0, fit1)
  expect_s3_class(tab, "anova")
  expect_equal(nrow(tab), 2)

  expect_error(vcov(fit1), "requires the original data matrix")
})

test_that("fit_ad dispatches to family-specific fitters", {
  set.seed(14)
  y <- simulate_gau(n_subjects = 20, n_time = 4, order = 1, phi = 0.25)
  fit <- fit_ad(y, family = "gaussian", order = 1)

  expect_s3_class(fit, "gau_fit")
  expect_s3_class(fit, "ad_fit")
})

test_that("alternative fit entry points retain ad_fit inheritance", {
  set.seed(15)
  y_g <- simulate_gau(n_subjects = 20, n_time = 4, order = 1, phi = 0.25)
  y_g[1, 2] <- NA
  expect_s3_class(em_gau(y_g, order = 1, max_iter = 2), "ad_fit")

  y_c <- simulate_cat(n_subjects = 30, n_time = 4, order = 1, n_categories = 2)
  y_c[1, 2] <- NA
  fit_c_em <- suppressWarnings(em_cat(y_c, order = 1, max_iter = 2))
  expect_s3_class(fit_c_em, "ad_fit")
})
