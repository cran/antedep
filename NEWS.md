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
