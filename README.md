# antedep

`antedep` fits antedependence models for longitudinal data:

- Gaussian AD models for continuous outcomes
- INAD models for count outcomes
- Categorical AD (CAT) models for discrete state outcomes

## Installation

```r
# install.packages("remotes")
remotes::install_github("TanchyKing/antedep")
```

## Production-Readiness Matrix

| Model | Data type | Complete-data fit/logLik | Missing-data fit/logLik | Notes |
|---|---|---|---|---|
| AD | Continuous | Ready | Ready (`fit_gau`, `logL_gau`) | Missing-data fit uses EM or observed-data likelihood modes |
| INAD | Counts | Ready | Ready (`fit_inad`, `logL_inad`) | Missing-data fit supports `na_action = "marginalize"` |
| CAT | Categorical states | Ready | Ready (`fit_cat`, `logL_cat`) | Missing-data fit supports orders 0, 1, 2 |

## Quick Start

Included datasets:

- `bolus_inad` (morphine bolus counts)
- `cattle_growth` (book companion cattle growth data; continuous)
- `cochlear_implant` (book companion speech recognition/cochlear data; continuous)
- `labor_force_cat` (labor-force categorical table expanded to subject sequences)
- `race_100km` (100km race split times; continuous)

### AD example (continuous)

```r
library(antedep)
set.seed(1)

y <- simulate_gau(n_subjects = 80, n_time = 5, order = 1)
fit <- fit_gau(y, order = 1)
fit$log_l
```

### Missing-data workflow

```r
library(antedep)
set.seed(1)

y <- simulate_inad(n_subjects = 60, n_time = 5, order = 1)
y[sample(length(y), 20)] <- NA

# Fit observed-data likelihood under MAR
fit_miss <- fit_inad(y, order = 1, na_action = "marginalize")
fit_miss$log_l
```

### CAT missing-data workflow

```r
library(antedep)
set.seed(1)

y_cat <- simulate_cat(n_subjects = 80, n_time = 5, order = 1, n_categories = 3)
y_cat[sample(length(y_cat), 30)] <- NA

# Observed-data likelihood (orders 0/1/2)
fit_cat_marg <- fit_cat(y_cat, order = 1, na_action = "marginalize")

# EM (orders 0/1) via explicit entry point or fit_cat dispatch
fit_cat_em1 <- em_cat(y_cat, order = 1, max_iter = 80)
fit_cat_em2 <- fit_cat(y_cat, order = 1, na_action = "em", em_max_iter = 80)
```

If EM becomes unstable or converges slowly, try increasing `epsilon`, increasing
`max_iter`, and using `safeguard = TRUE`.

## Known Limitations

- EM entry points: `em_gau` (Gaussian), `em_inad` (INAD), and `em_cat` (CAT, orders 0/1) are available; for CAT order 2 with missing data, use `fit_cat(na_action = "marginalize")`.
- For CAT models, `fit_cat()` supports both observed-data likelihood (`na_action = "marginalize"`) and EM (`na_action = "em"` for orders 0/1).
- Missing-data confidence intervals are not yet implemented (`ci_gau`, `ci_inad`, `ci_cat` require complete-data fits).
- AD missing-data EM for order 2 is currently available with simplified updates and should be treated as provisional.
- AD missing-data LRT/mean/covariance tests remain complete-data only.
- CAT missing-data stationarity/time-invariance tests (`lrt_stationarity_cat`, `lrt_timeinvariance_cat`, `run_stationarity_tests_cat`) remain complete-data only.
- INAD missing-data LRTs and CAT missing-data order/homogeneity LRTs are supported via observed-data likelihood (`na_action = "marginalize"`).
- Current model order support is up to order 2 for AD, INAD, and CAT.
- CAT missing-data marginalization currently supports orders 0, 1, and 2.

## Vignette

Source vignette: `vignettes/antedep-intro.Rmd`.
Rendered site article: `docs/articles/antedep-intro.html`.
Integrated function reference: `docs/reference/index.html`.

## Function Reference Site (pkgdown)

- Local build: `Rscript -e "pkgdown::build_site(preview = FALSE)"`
- CI deployment: `.github/workflows/pkgdown.yml` (GitHub Pages via Actions)

## Local Check Guide

- See `LOCAL_CHECKS.md` for tarball-mode vs directory-mode checks and expected NOTE messages.
