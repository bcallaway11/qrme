# Changelog

## qrme 1.0.1

### New functionality

- **[`tsme_model_select()`](https://bcallaway11.github.io/qrme/reference/tsme_model_select.md)**
  — grid search over copula families (Gaussian, Clayton, Gumbel, Frank)
  and ME distributions (Gaussian, Laplace), returning AIC/BIC for each
  combination to guide model selection.
- **[`logLik.merr()`](https://bcallaway11.github.io/qrme/reference/logLik.merr.md)**,
  **`AIC.merr()`**, **`BIC.merr()`** — log-likelihood and information
  criteria for `merr` objects, enabling mixture-order selection via
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html).
- **[`print.tsme()`](https://bcallaway11.github.io/qrme/reference/print.tsme.md)**,
  **[`summary.tsme()`](https://bcallaway11.github.io/qrme/reference/summary.tsme.md)**
  — structured summaries of `tsme` results, including copula parameters,
  ME distributions, transition matrices, and upward mobility estimates.
- **[`autoplot.tsme()`](https://bcallaway11.github.io/qrme/reference/autoplot.tsme.md)**,
  **[`plot.tsme()`](https://bcallaway11.github.io/qrme/reference/plot.tsme.md)**
  — `ggplot2`-based visualisation of conditional quantile curves and
  poverty rates for ME-corrected, no-ME, and naive QR estimators.
- **Frank copula** added to
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) and
  [`qr2me()`](https://bcallaway11.github.io/qrme/reference/qr2me.md)
  (previously only Gaussian, Clayton, and Gumbel were supported).
- Bundled **`nlsy97`** (NLSY97 father-son income pairs) and
  **`nlsy97_tsme_fit`** (pre-computed
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md)
  result) as package datasets, used in the new vignettes and examples.

### Improvements in tuning parameters

- **Data-informed starting values** — `start_sigma` is now initialised
  so that the mixture standard deviation equals 0.5 × sd(Y) (instead of
  a fixed constant); `start_mu` is evenly spaced over ±0.25 × sd(Y) for
  `n_mix > 1` and 0 for `n_mix = 1`.
- **Adaptive `proposal_sd`** — the Metropolis–Hastings proposal standard
  deviation is initialised from the current ME parameters and updated
  each EM iteration to track the evolving ME scale, reducing the need to
  set it manually.
- **MH acceptance rate tracking** — the acceptance rate is reported at
  `verbose = 2` and a warning is issued if the final-iteration rate
  falls outside 10–90%, prompting the user to pass `proposal_sd`
  explicitly.
- **`conv_patience`** — new argument controlling how many consecutive
  iterations below `tol` are required before declaring convergence
  (default 1); raising to 2 is recommended when
  `conv_criterion = "loglik"`.
- **Bootstrap warm-start** — bootstrap replications now initialise ME
  parameters from the full-data fit rather than from scratch,
  substantially reducing bootstrap run time and improving stability.
- **Separate bootstrap tuning** —
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md)
  exposes `boot_mcmc_draws`, `boot_mcmc_burn_in`,
  `boot_n_copula_me_draws`, and `boot_tol` so that bootstrap
  replications can use lighter settings than the main fit without
  affecting point-estimate accuracy.

### Other changes

- Argument names in
  [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md)
  standardised to snake_case: `max_iter` → `max_em_iters`, `burn_in` →
  `mcmc_burn_in`, `n_draws` → `mcmc_draws`.
- Output fields in
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md)
  renamed to snake_case throughout (e.g., `meTmat` → `me_tmat`,
  `meCopParam` → `me_cop_param`).
- Conditional quantile curves in
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) now
  evaluated at the mean of the covariate matrix rather than a grid,
  matching the paper’s specification.
- Laplace M-step corrected to use the MLE scale estimator instead of the
  normal-mixture variance formula.
- Clayton copula conditional CDF clamped to avoid numerical instability
  at boundary quantiles.
- Added Emmanuel S. Tsyawo (<estsyawo@gmail.com>) as a listed author.
- Switched from `exportPattern` to explicit `@export` tags; C++
  internals are no longer part of the public API.
- Removed the `R/old/` legacy directory.
- Updated minimum R version to 4.0.0; added `URL` and `BugReports` to
  DESCRIPTION.
- Added two vignettes (Quarto): an introduction to
  [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md) and
  an application of
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) to
  intergenerational income mobility using the NLSY97 data.
- Added a pkgdown website (<https://bcallaway11.github.io/qrme/>).

## qrme 1.0.0

- Initial release.
