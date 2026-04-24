# test-em-cat.R - Tests for EM estimation in CAT models

test_that("em_cat returns closed-form fit when no missing data", {
  skip_on_cran()
  set.seed(1301)
  y <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 2)

  fit_em <- em_cat(y, order = 1)
  fit_cf <- fit_cat(y, order = 1)

  expect_s3_class(fit_em, "cat_fit")
  expect_equal(as.numeric(fit_em$log_l), as.numeric(fit_cf$log_l), tolerance = 1e-10)
  expect_equal(fit_em$settings$method, "closed_form")
})


test_that("em_cat fits missing-data order 0 model", {
  skip_on_cran()
  set.seed(1302)
  y <- simulate_cat(n_subjects = 40, n_time = 5, order = 0, n_categories = 3)
  y[sample(length(y), 25)] <- NA

  fit <- em_cat(y, order = 0, max_iter = 20, tol = 1e-7)

  expect_s3_class(fit, "cat_fit")
  expect_true(is.finite(fit$log_l))
  expect_equal(fit$settings$order, 0)
  expect_true(all(vapply(fit$marginal, function(x) abs(sum(x) - 1) < 1e-8, logical(1))))
  expect_false(is.null(fit$cell_counts))
  expect_equal(names(fit$cell_counts), paste0("t", 1:ncol(y)))
})


test_that("em_cat handles non-sequential block IDs and stores block levels", {
  skip_on_cran()
  set.seed(1303)
  y <- simulate_cat(n_subjects = 40, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 20)] <- NA
  blocks <- rep(c(2, 5, 2), length.out = nrow(y))

  fit <- em_cat(y, order = 1, blocks = blocks, homogeneous = FALSE, max_iter = 15)

  expect_equal(fit$settings$n_blocks, 2)
  expect_equal(sort(fit$settings$block_levels), c("2", "5"))
  expect_equal(length(fit$marginal), 2)
  expect_equal(length(fit$transition), 2)
  expect_equal(names(fit$cell_counts), c("block_1", "block_2"))
  expect_true("t1_to_t1" %in% names(fit$cell_counts$block_1))
})


test_that("em_cat returns fit_cat-style order 1 cell_counts usable by ci_cat", {
  skip_on_cran()
  set.seed(1307)
  y <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 25)] <- NA

  fit <- em_cat(y, order = 1, max_iter = 15, tol = 1e-6)

  expected_names <- c("t1_to_t1", paste0("t", 1:3, "_to_t", 2:4))
  expect_equal(names(fit$cell_counts), expected_names)
  expect_true(is.matrix(fit$cell_counts$t1_to_t2))

  ci <- ci_cat(fit)
  expect_s3_class(ci, "cat_ci")
})


test_that("em_cat order 1 does not report immediate convergence on first iteration", {
  skip_on_cran()
  set.seed(1304)
  y <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 15)] <- NA

  fit <- em_cat(y, order = 1, max_iter = 20, tol = 1e-6)

  expect_true(fit$settings$n_iter > 1)
})


test_that("em_cat transition probabilities are valid for order 1", {
  skip_on_cran()
  set.seed(1308)
  y <- simulate_cat(n_subjects = 40, n_time = 4, order = 1, n_categories = 3)
  y[sample(length(y), 20)] <- NA

  fit <- em_cat(y, order = 1, max_iter = 15)

  for (t in 2:ncol(y)) {
    trans_t <- fit$transition[[paste0("t", t)]]
    expect_true(all(trans_t >= 0))
    expect_equal(rowSums(trans_t), rep(1, nrow(trans_t)), tolerance = 1e-8)
  }
})


test_that("em_cat homogeneous and one-block heterogeneous fits agree", {
  skip_on_cran()
  set.seed(1309)
  y <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 25)] <- NA
  blocks_one <- rep(5, nrow(y))

  fit_h <- em_cat(y, order = 1, homogeneous = TRUE, max_iter = 15)
  fit_het_one <- em_cat(y, order = 1, blocks = blocks_one, homogeneous = FALSE, max_iter = 15)

  expect_equal(as.numeric(fit_h$log_l), as.numeric(fit_het_one$log_l), tolerance = 1e-8)
  expect_equal(as.numeric(fit_h$marginal$t1), as.numeric(fit_het_one$marginal$t1), tolerance = 1e-8)
})


test_that("em_cat accepts safeguard argument", {
  skip_on_cran()
  set.seed(1310)
  y <- simulate_cat(n_subjects = 40, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 18)] <- NA

  fit_safe <- em_cat(y, order = 1, safeguard = TRUE, max_iter = 15)
  fit_no_safe <- em_cat(y, order = 1, safeguard = FALSE, max_iter = 15)

  expect_s3_class(fit_safe, "cat_fit")
  expect_s3_class(fit_no_safe, "cat_fit")
})


test_that("em_cat order 2 reports not implemented", {
  skip_on_cran()
  set.seed(1305)
  y <- simulate_cat(n_subjects = 20, n_time = 4, order = 2, n_categories = 2)
  y[sample(length(y), 5)] <- NA

  expect_error(
    em_cat(y, order = 2),
    "not yet implemented"
  )
})


test_that("em_cat final log-likelihood is consistent when max_iter is reached", {
  skip_on_cran()
  set.seed(1306)
  y <- simulate_cat(n_subjects = 50, n_time = 4, order = 1, n_categories = 2)
  y[sample(length(y), 30)] <- NA

  expect_warning(
    fit <- em_cat(y, order = 1, max_iter = 1, tol = 1e-12),
    "did not converge"
  )

  ll_check <- logL_cat(
    y = y,
    order = 1,
    marginal = fit$marginal,
    transition = fit$transition,
    n_categories = 2,
    na_action = "marginalize"
  )

  expect_equal(as.numeric(fit$log_l), as.numeric(ll_check), tolerance = 1e-8)
  expect_equal(
    as.numeric(fit$aic),
    as.numeric(-2 * fit$log_l + 2 * fit$n_params),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(fit$bic),
    as.numeric(-2 * fit$log_l + log(nrow(y)) * fit$n_params),
    tolerance = 1e-8
  )
})
