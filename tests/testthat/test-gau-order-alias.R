test_that("test_one_sample_gau accepts order alias for p", {
  skip_on_cran()
  set.seed(1001)
  y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.3)
  mu0 <- rep(0, ncol(y))

  out_p <- test_one_sample_gau(y, mu0 = mu0, p = 1, use_modified = FALSE)
  out_order <- test_one_sample_gau(y, mu0 = mu0, order = 1, use_modified = FALSE)

  expect_equal(out_order$order, out_p$order)
  expect_equal(out_order$statistic, out_p$statistic)
  expect_equal(out_order$p_value, out_p$p_value)
})

test_that("test_two_sample_gau accepts order alias for p", {
  skip_on_cran()
  set.seed(1002)
  y <- simulate_gau(n_subjects = 50, n_time = 5, order = 1, phi = 0.25)
  blocks <- c(rep(1, 25), rep(2, 25))

  out_p <- test_two_sample_gau(y, blocks = blocks, p = 1, use_modified = FALSE)
  out_order <- test_two_sample_gau(y, blocks = blocks, order = 1, use_modified = FALSE)

  expect_equal(out_order$order, out_p$order)
  expect_equal(out_order$statistic, out_p$statistic)
  expect_equal(out_order$p_value, out_p$p_value)
})

test_that("test_order_gau accepts order_null and order_alt aliases", {
  skip_on_cran()
  set.seed(1004)
  y <- simulate_gau(n_subjects = 40, n_time = 5, order = 1, phi = 0.3)

  out_pq       <- test_order_gau(y, p = 0, q = 1, use_modified = FALSE)
  out_aliases  <- test_order_gau(y, order_null = 0, order_alt = 1,
                                 use_modified = FALSE)

  expect_equal(out_aliases$statistic, out_pq$statistic)
  expect_equal(out_aliases$p_value,   out_pq$p_value)
  expect_equal(out_aliases$order_null, 0L)
  expect_equal(out_aliases$order_alt,  1L)
})

test_that("test_*_gau p and order cannot both be supplied", {
  skip_on_cran()
  set.seed(1003)
  y <- simulate_gau(n_subjects = 30, n_time = 4, order = 1, phi = 0.2)
  blocks <- c(rep(1, 15), rep(2, 15))

  expect_error(
    test_one_sample_gau(y, mu0 = rep(0, ncol(y)), p = 1, order = 1),
    "Specify only one of 'p' or 'order'"
  )
  expect_error(
    test_two_sample_gau(y, blocks = blocks, p = 1, order = 1),
    "Specify only one of 'p' or 'order'"
  )
})
