test_that("predict.inad_fit returns mean and sample forecasts with stable shapes", {
  set.seed(20260520)
  y <- simulate_inad(
    n_subjects = 25,
    n_time = 5,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    alpha = 0.25,
    theta = 1.5
  )
  fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois")

  pred_mean <- predict(fit, h = 1, type = "mean", n_sims = 20, seed = 1)
  expect_equal(dim(pred_mean), c(nrow(y), 1L))
  expect_true(all(is.finite(pred_mean)))

  pred_sample <- predict(fit, newdata = y[, 1:3, drop = FALSE],
                         h = 2, type = "sample", n_sims = 10, seed = 2)
  expect_equal(dim(pred_sample), c(nrow(y), 2L, 10L))
  expect_true(all(pred_sample >= 0))
  expect_true(all(pred_sample == floor(pred_sample)))
})

test_that("predict.inad_fit is reproducible when a seed is supplied", {
  set.seed(20260520)
  y <- simulate_inad(
    n_subjects = 20,
    n_time = 4,
    order = 1,
    thinning = "pois",
    innovation = "pois",
    alpha = 0.2,
    theta = 1.2
  )
  fit <- fit_inad(y, order = 1, thinning = "pois", innovation = "pois")

  pred1 <- predict(fit, newdata = y[, 1:3, drop = FALSE],
                   h = 1, type = "sample", n_sims = 15, seed = 99)
  pred2 <- predict(fit, newdata = y[, 1:3, drop = FALSE],
                   h = 1, type = "sample", n_sims = 15, seed = 99)
  expect_equal(pred1, pred2)
})

test_that("predict.inad_fit validates horizons and histories", {
  set.seed(20260520)
  y <- simulate_inad(
    n_subjects = 12,
    n_time = 4,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    alpha = 0.3,
    theta = 1.4
  )
  fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "pois")

  expect_error(predict(fit, h = 4, n_sims = 5), "h must be smaller")
  bad_history <- y[, 1:3, drop = FALSE]
  bad_history[1, 1] <- NA
  expect_error(predict(fit, newdata = bad_history, h = 1, n_sims = 5),
               "cannot contain missing")
  expect_error(predict(fit, newdata = y[, 1:4, drop = FALSE], h = 1, n_sims = 5),
               "exceeds the fitted time grid")
})
