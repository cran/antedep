## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(antedep)
set.seed(123)

## ----eval=FALSE---------------------------------------------------------------
# data("race_100km")
# dim(race_100km$y)
# colMeans(race_100km$y)

## ----warning=FALSE, message=FALSE, fig.width=10, fig.height=4.5, out.width='100%'----
data("bolus_inad")
plot_profile(
  bolus_inad$y,
  blocks = bolus_inad$blocks,
  block_labels = c("Group 1", "Group 2"),
  title = "Bolus Longitudinal Profiles"
)

## ----eval=FALSE---------------------------------------------------------------
# # More detailed dependence diagnostic
# plot_prism(bolus_inad$y)

## ----warning=FALSE, message=FALSE---------------------------------------------
y_gau <- simulate_gau(n_subjects = 20, n_time = 4, order = 1, phi = 0.5)
fit_gau1 <- fit_gau(y_gau, order = 1)
fit_gau2 <- fit_gau(y_gau, order = 2)

c(order1_logLik = fit_gau1$log_l, order2_logLik = fit_gau2$log_l)
bic_order_gau(y_gau, max_order = 2)$table

## -----------------------------------------------------------------------------
# Example AD order test
test_order_gau(y_gau, p = 0, q = 1)$p_value

## -----------------------------------------------------------------------------
y_gau_miss <- y_gau
y_gau_miss[sample(length(y_gau_miss), 8)] <- NA

fit_gau_em <- fit_gau(y_gau_miss, order = 1, na_action = "em", em_max_iter = 10)
fit_gau_cc <- fit_gau(y_gau_miss, order = 1, na_action = "complete")

c(em_logLik = fit_gau_em$log_l, complete_case_logLik = fit_gau_cc$log_l)
fit_gau_em$pct_missing

## -----------------------------------------------------------------------------
# One-sample and two-sample mean-profile tests
mu0 <- colMeans(y_gau)
blocks_gau <- rep(1:2, each = nrow(y_gau) / 2)

test_one_sample_gau(y_gau, mu0 = mu0, p = 1)$p_value
test_two_sample_gau(y_gau, blocks = blocks_gau, p = 1)$p_value

## -----------------------------------------------------------------------------
# Wald contrast test: H0 that adjacent means are equal
C <- diff(diag(ncol(y_gau)))
test_contrast_gau(y_gau, C = C, p = 1)$p_value

## -----------------------------------------------------------------------------
# Covariance homogeneity across groups
test_homogeneity_gau(y_gau, blocks = blocks_gau, p = 1)$p_value

## -----------------------------------------------------------------------------
y_inad <- simulate_inad(
  n_subjects = 20,
  n_time = 4,
  order = 1,
  thinning = "binom",
  innovation = "pois",
  alpha = 0.4,
  theta = 3
)

fit_inad1 <- fit_inad(y_inad, order = 1, thinning = "binom", innovation = "pois")
fit_inad1$log_l
bic_order_inad(y_inad, max_order = 2, thinning = "binom", innovation = "pois")$table

## -----------------------------------------------------------------------------
y_inad_miss <- y_inad
y_inad_miss[sample(length(y_inad_miss), 10)] <- NA

fit_inad_miss <- fit_inad(
  y_inad_miss,
  order = 1,
  thinning = "binom",
  innovation = "pois",
  na_action = "marginalize",
  max_iter = 10
)

fit_inad_miss$log_l
fit_inad_miss$pct_missing

## -----------------------------------------------------------------------------
# Missing-data LRTs (observed-data likelihood)
blocks_inad <- rep(1:2, each = nrow(y_inad_miss) / 2)
test_order_inad(y_inad_miss, order_null = 0, order_alt = 1, blocks = blocks_inad)$p_value
test_homogeneity_inad(y_inad_miss, blocks = blocks_inad, order = 1, test = "mean")$p_value

## -----------------------------------------------------------------------------
y_cat <- simulate_cat(n_subjects = 20, n_time = 4, order = 1, n_categories = 3)
fit_cat1 <- fit_cat(y_cat, order = 1)
fit_cat1$log_l

# Wald CIs (complete data)
ci_cat1 <- ci_cat(fit_cat1, parameters = "marginal")
head(ci_cat1$marginal, 3)

## -----------------------------------------------------------------------------
# Order tests for CAT
run_order_tests_cat(y_cat, max_order = 2)$table

## -----------------------------------------------------------------------------
y_cat_miss <- y_cat
y_cat_miss[sample(length(y_cat_miss), 12)] <- NA

fit_cat_miss <- fit_cat(y_cat_miss, order = 1, na_action = "marginalize")
fit_cat_miss$log_l
fit_cat_miss$settings$na_action_effective

## -----------------------------------------------------------------------------
# EM alternative for CAT missing data (orders 0/1)
fit_cat_miss_em <- fit_cat(y_cat_miss, order = 1, na_action = "em", em_max_iter = 10)
fit_cat_miss_em$log_l

## -----------------------------------------------------------------------------
# Missing-data order/homogeneity tests are supported
blocks_cat <- rep(1:2, each = nrow(y_cat_miss) / 2)
test_order_cat(y_cat_miss, order_null = 0, order_alt = 1)$p_value
test_homogeneity_cat(y_cat_miss, blocks = blocks_cat, order = 1)$p_value

## ----eval=FALSE---------------------------------------------------------------
# # INAD homogeneity/stationarity tools
# blocks <- rep(1:2, each = nrow(y_inad) / 2)
# run_homogeneity_tests_inad(y_inad, blocks = blocks, order = 1)
# run_stationarity_tests_inad(y_inad, order = 1, blocks = blocks)
# 
# # CAT stationarity/time-invariance tools
# test_timeinvariance_cat(y_cat, order = 1)
# test_stationarity_cat(y_cat, order = 1)

## ----eval=FALSE---------------------------------------------------------------
# Rscript scripts/check-package.R --as-cran --no-manual

