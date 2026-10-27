# antedep 0.3.1

## Documentation
- The examples of `ci_inad()`, `test_homogeneity_inad()`, and
  `run_homogeneity_tests_inad()` use the first four occasions of the bolus
  data and Poisson innovations, so that they run in seconds;
  `test_homogeneity_inad()` now shows how fits are reused across tests.
- The reference in `test_homogeneity_inad()` now cites the INAD manuscript as
  under review.
- No change in computations: every function returns the same results as in
  0.3.0.

# antedep 0.3.0

## New features
- `fit_inad()` gains a `constrain` argument for time-invariance constraints
  (`"alpha"`, `"alpha1"`, `"alpha2"`, `"theta"`, `"kappa"`, `"innovation"`,
  `"both"`, `"all"`). Constrained models are refitted by maximizing the
  observed-data log-likelihood jointly over every remaining free parameter,
  including negative binomial innovation sizes, and are finished with Newton
  steps using the exact score and observed information. `n_params`, `AIC`,
  and `BIC` count only the free parameters.
- `ci_inad()` now computes the observed information exactly by Louis' identity
  for all free parameters jointly: models of order 0, 1, and 2, with or without
  block effects, with or without time-invariance constraints, and including
  the negative binomial innovation size. `method = "hessian"` gives the
  numerical observed-data Hessian (fixed absolute steps) as a check. Intervals
  for positive parameters are formed on the log scale by default
  (`log_scale = TRUE`); the full covariance matrix is returned as `vcov`.
  Block effects receive Wald intervals in `tau`; for unconstrained fits,
  profile-likelihood intervals are also returned in `tau_profile`.
- `em_inad()` is now an exact EM algorithm for all nine thinning-innovation
  combinations. The E-step evaluates the full conditional distribution of the
  latent thinned counts, so the expected complete-data log-likelihood is exact
  even for negative binomial innovations, whose log-likelihood is not linear in
  the latent counts; negative binomial innovation sizes are now estimated.
  With block effects the innovation parameters are updated by conditional
  maximization (ECM). The returned fit includes `n_params`, `aic`, and `bic`.
- `test_stationarity_inad()` and `run_stationarity_tests_inad()` support
  constraints on negative binomial innovation sizes (`"kappa"`,
  `"innovation"`, `"all"`).
- Fit objects now inherit from a shared `ad_fit` parent class while retaining
  their family-specific classes (`gau_fit`, `cat_fit`, and `inad_fit`).
- Added `fit_ad()` as a unified fitting entry point for Gaussian, categorical,
  and INAD workflows.
- Added `anova()` methods for likelihood-ratio comparison of nested
  antedependence fits. For INAD fits of different orders no chi-square
  p-value is given, because the comparison is a boundary problem.
- New `kkt_inad()` evaluates the exact score of a fitted INAD model and checks
  the first-order (Karush-Kuhn-Tucker) conditions, reporting estimates on the
  boundary of the parameter space and on the negative binomial size cap.
  Interior stationarity requires a positive definite observed information,
  a Newton decrement and single-parameter quantities s_j^2 / I_jj no larger
  than the tolerance. A boundary estimate whose one-sided score points into
  the parameter space passes only if a line search on the observed-data
  log-likelihood along that parameter finds no increase above the tolerance.
- `fit_inad()` gains `newton_polish = TRUE`: complete-data fits are finished
  by Newton steps with the exact score and observed information, holding
  boundary estimates whose score points outward fixed. A coefficient exactly
  on a bound where its thinned count is degenerate, and whose one-sided score
  points inward, is first moved off the bound by a line search on the
  observed-data log-likelihood.

## Changes in behavior
- `test_order_inad()` no longer reports an ordinary chi-square p-value, which
  is not valid when the added thinning coefficients are on the boundary under
  the null hypothesis. By default (`reference = "none"`) it reports the
  likelihood ratio statistic and BIC values only; `reference = "chibar"` gives a
  working binomial-weight chi-bar-square p-value, with a warning for models
  with block effects, for which these weights have not been established.
  `use_chibar` is deprecated.

## Bug fixes
- Direct maximization in `fit_inad()` (models with and without block
  effects, constrained fits, and missing-data fits) restricted every thinning
  coefficient to be below 1. Only binomial thinning requires this; Poisson and
  negative binomial thinning coefficients may now exceed 1, as in `em_inad()`,
  so that both algorithms maximize over the same parameter space.
  `simulate_inad()` likewise accepts negative binomial thinning coefficients of
  1 or more.
- For models without block effects, `fit_inad()` eliminated the innovation
  mean through the sample moment equation and required the implied innovation
  mean to be positive, so maxima with the innovation mean on its boundary
  (zero) were missed. Each time point's parameters are now optimized jointly
  from several starting values (including the nested lower-order fit), so that
  an order-2 fit is never worse than the order-1 fit. For order-2 Bell models,
  the reported Bell parameter at the second occasion was the innovation mean.
- The missing-data likelihood (`na_action = "marginalize"`) restricted each
  thinned count to be no larger than the count being thinned, which is wrong
  for Poisson and negative binomial thinning; fits to incomplete data with
  these operators were therefore incorrect.
- A thinning coefficient whose predecessor counts are all zero is not
  identifiable; it is now set to zero by both `fit_inad()` and `em_inad()`.
- The final Newton steps could stall when a coefficient was just above zero
  with its score pointing out of the parameter space (as after EM, which
  approaches boundary values slowly), when a parameter on its bound had a
  Newton step pointing outward, or far from the maximum, where the observed
  information need not be positive definite. Such coefficients are now moved
  onto or held at the bound, and damped (Levenberg-Marquardt) steps are used
  when the Newton step fails. Constrained fits record in `settings$path`
  whether they were completed by these Newton steps or needed the
  derivative-free search.
- At a thinning coefficient of exactly 0 (or binomial 1) the thinned count is
  degenerate, and the conditional mean of the complete-data score is 0/0. The
  score routines replaced the coefficient by 1e-10, which gives minus the
  predecessor count instead of the one-sided derivative, and was also wrong
  for coefficients between 0 and 1e-10. As a result `kkt_inad()` could accept
  a fit with a zero coefficient that a small increase would improve, and the
  Newton steps held such coefficients at zero. The exact one-sided derivative
  is now used. `kkt_inad()` also accepted a negative Newton decrement when the
  observed information was indefinite; it now requires positive definiteness.
- The negative binomial thinning pmf was evaluated with `prob = 1/(1 + alpha)`,
  from which `dnbinom()` recovers `alpha/(1 + alpha)` by subtraction. For small
  alpha this loses about -log10(alpha) digits in the probabilities of positive
  thinned counts. At alpha = 1e-12, where EM leaves a coefficient at its lower
  limit, the score was off by about 0.4%; log-likelihoods were affected only by
  rounding. The pmf now uses the equivalent mean parameterization.
- Constrained fits in `test_stationarity_inad()` held negative binomial
  innovation sizes fixed at their unconstrained estimates instead of
  re-estimating them, giving incorrect likelihood ratio statistics and BIC
  values for negative binomial innovations.
- `.chi_bar_pvalue()` (used by `test_order_inad(use_chibar = TRUE)`) counted the
  point mass at zero as contributing 1 to the upper tail, inflating every
  chi-bar-square p-value by 2^-r.
- Parameter counts in the tables of `test_order_inad()`, `bic_order_inad()`,
  and `test_stationarity_inad()` now include one negative binomial innovation
  size per time point and the coefficient alpha_{2,1} of order-2 models, as in
  `fit_inad()`. For order 2, the `"alpha1"` constraint now includes
  alpha_{2,1} instead of fixing it at zero.
- `em_inad()` restricted the latent thinned counts to be no larger than the
  count being thinned, which is wrong for Poisson and negative binomial
  thinning, and capped their coefficients below 1.
- The Bell parameter is now obtained from the Bell mean to machine precision
  (previously `uniroot()` default tolerance), so Bell log-likelihoods with
  block effects are smooth in the parameters.
- Louis-based intervals in `ci_inad()` for order-2 models used only the lag-1
  coefficient, and Bell intervals with block effects omitted the chain rule
  from the innovation mean to the Bell parameter.

## Improvements
- Constrained fits start from pooled estimates and are completed by Newton
  steps; the derivative-free algorithm is used only when the first-order
  conditions are not met, which makes constrained fits much faster.
- Block-effect `tau` standard errors and confidence intervals are aligned in
  downstream tables (Wald in `tau`; profile likelihood in `tau_profile`).
- `coef()` now returns a named numeric vector by default for standard R
  compatibility; use `type = "list"` for the previous structured inspection
  format.
- `confint()` now returns a standard two-column matrix for all fit families;
  the structured `ci_gau()`, `ci_cat()`, and `ci_inad()` helpers remain
  available.
- Added broader `vcov()` support using confidence-interval standard errors
  where available, improving compatibility with tools such as
  `lmtest::coeftest()`.
- The family-specific `aic_*()` and `bic_*()` helpers are retained as thin
  compatibility wrappers for existing code, but examples and new documentation
  prefer the standard `AIC()` and `BIC()` generics.

# antedep 0.2.0

## New features
- New S3 methods `deviance()`, `confint()`, and `vcov()` for `gau_fit`,
  `cat_fit`, and `inad_fit` objects, improving compatibility with standard
  R model-fitting workflows.
- New `summary.partial_corr()` method prints a lag-by-lag table of
  intervenor-adjusted partial correlations (mean absolute value, range, and
  number of significant pairs at each lag).

## Improvements
- Matrix inversion in `ci_inad()`, `fit_gau_em()`, and `cat_test_stats.R`
  now uses `chol2inv(chol(.))` instead of `solve()` for symmetric
  positive-definite matrices, improving numerical stability and efficiency.

# antedep 0.1.0

## New features
- `fit_inad()` gains a `nb_inno_size_ub` argument (default 50) that caps the
  upper bound of the negative-binomial innovation size parameter during
  optimization, improving numerical stability for near-Poisson data.
- `test_order_gau()` accepts `order_null` and `order_alt` as convenience aliases
  for `p` and the absolute alternative order; both are also returned in the
  result object.

## Bug fixes
- `ci_inad()`: fixed a sign error in the observed Fisher information for the
  negative-binomial innovation size parameter; the Hessian term
  `(r + u) / (r + λ)²` was added instead of subtracted, producing confidence
  intervals that were too wide.
- `ci_inad()`: the numerical second derivative for `nb_inno_size` CIs now
  retries with progressively smaller step sizes (×0.1, ×0.01) before falling
  back to NA, avoiding spurious failures when the default step lands in a
  non-finite region.
- `test_homogeneity_inad()`: degrees of freedom for LRT tests involving
  `innovation = "nbinom"` are now computed from the actual number of NB size
  parameters in the fitted models rather than assuming a fixed count of 1.
  This corrects LRT statistics and p-values whenever `nb_inno_size` is fitted
  as a time-varying vector.
- `ci_inad()` tau profile CI: `nb_inno_size` (negative-binomial innovation
  dispersion) is now held fixed at its full-model MLE during profile refits,
  consistent with the constrained-fit paradigm used throughout the package.
  Previously it was re-optimised as a nuisance parameter, which could widen
  the interval to the point of crossing zero even when the LRT clearly rejects
  the null (Variant 1 vs Variant 2 fix).
- `ci_inad()` tau profile CI: the bracket search in `.ci_tau_profile_inad`
  no longer imposes an artificial upper cap (`max(|tau_mle| + 1, 1)`) on the
  search range. The maximum bracket iterations are increased from 20 to 50
  and the initial step size is set to `max(0.1, |tau_mle| * 0.2)`, preventing
  the search from stalling for large or near-zero MLEs.

## Initial release
- Initial CRAN submission candidate for Gaussian AD, categorical AD, and INAD workflows.
- Added/expanded examples for key user-facing modeling functions (`fit_*`, `em_*`, `simulate_*`, `logL_*`).
- Refreshed missing-data notes and harmonized documentation links and metadata.
- `logL_gau()` default missing-data behavior is now `na_action = "fail"` (previously
  marginalization-first in earlier drafts). For missing inputs, pass
  `na_action = "marginalize"` or `na_action = "complete"` explicitly.
- Added packaged datasets `labor_force_cat` (categorical labor-force sequences) and
  `race_100km` (continuous 100km race split times).
