test_that("em_inad preserves original block labels in settings", {
  skip_on_cran()
  set.seed(77)
  y <- matrix(rpois(20 * 4, lambda = 2), nrow = 20, ncol = 4)
  blocks <- rep(c("grp_5", "grp_2"), each = 10)

  fit <- em_inad(
    y,
    order = 1,
    thinning = "binom",
    innovation = "pois",
    blocks = blocks,
    max_iter = 3,
    tol = 1e-5,
    verbose = FALSE
  )

  expect_s3_class(fit, "inad_fit")
  expect_equal(sort(fit$settings$block_levels), c("grp_2", "grp_5"))
  expect_equal(length(fit$tau), 2)
})
