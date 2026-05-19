# tsme

function to compute nonlinear models with two sided measurement error

## Usage

``` r
tsme(
  data,
  y_formula,
  t_formula,
  tau,
  t_values,
  x_data = NULL,
  me_distribution = "gaussian",
  copula = "gaussian",
  mcmc_method = "MH",
  n_copula_me_draws = 100L,
  report_transition_matrix = TRUE,
  report_spearman = TRUE,
  report_upward_mobility = TRUE,
  report_poverty = TRUE,
  pov_line = log(20000),
  report_quantiles = c(0.1, 0.5, 0.9),
  y_n_mix = 1,
  t_n_mix = 1,
  tol = NULL,
  conv_criterion = "params",
  loglik_draws = 100L,
  conv_patience = 1L,
  max_em_iters = 100L,
  mcmc_draws = 200,
  mcmc_burn_in = 100,
  proposal_sd = NULL,
  start_sigma_y = NULL,
  start_mu_y = NULL,
  start_pi_y = NULL,
  start_sigma_t = NULL,
  start_mu_t = NULL,
  start_pi_t = NULL,
  n_copula_draws = 100L,
  ignore_me = FALSE,
  verbose = FALSE,
  se = FALSE,
  n_boot = 100,
  n_cores = 1,
  boot_n_copula_me_draws = NULL,
  boot_mcmc_draws = NULL,
  boot_mcmc_burn_in = NULL,
  boot_tol = NULL
)
```

## Arguments

- data:

  data.frame

- y_formula:

  formula for outcome model

- t_formula:

  formula for treatment model

- tau:

  values of tau to estimate first step quantile regressions for

- t_values:

  values of the treatment to compute conditional distribution-type
  parameters for

- x_data:

  matrix of values of covariates to average over for conditional
  distribution-type parameters. The default is NULL and in this case all
  covariates in the data will be averaged over. A main alternative would
  be to pass in a single row with particular values of covariates of
  interest.

- me_distribution:

  the distribution of the measurement error. "gaussian" is the default
  and supports a mixture of normals. "laplace" is also supported

- copula:

  what type of copula to use in second step. Options are "gaussian" (the
  default), "clayton", "gumbel", or "frank"

- mcmc_method:

  simulation step to use in the EM algorithm. Currently only `"MH"` is
  supported.

- n_copula_me_draws:

  Number of measurement-error draws used in the second-stage copula
  likelihood (default 100). Increase this for a less noisy copula
  likelihood at the cost of speed.

- report_transition_matrix:

  whether or not to report a transition matrix

- report_spearman:

  whether or not to report Spearman's rho (rank-rank correlation)

- report_upward_mobility:

  whether or not to report upward mobility parameters

- report_poverty:

  whether or not to report fraction of population below the poverty line
  as a function of parents' income

- pov_line:

  value of the poverty line (default log(20000))

- report_quantiles:

  quantiles of child's income as a function of parents' income to report
  (default is .1,.5,.9)

- y_n_mix:

  Number of Gaussian mixture components for the outcome (Y) measurement
  error distribution. Default is 1 (single Gaussian). Use
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html), and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) on a fitted
  [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md)
  object to select the appropriate value. `y_n_mix` and `t_n_mix` can
  differ – it is valid and common to use different complexity for each
  equation.

- t_n_mix:

  Number of Gaussian mixture components for the treatment (T)
  measurement error distribution. Default is 1. See `y_n_mix`.

- tol:

  Convergence tolerance. When `NULL` (default), a value is chosen
  automatically based on `conv_criterion`. See
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md)
  for details.

- conv_criterion:

  Convergence criterion passed to
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md):
  `"params"` (default) or `"loglik"`. See
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md).

- loglik_draws:

  Number of Monte Carlo draws for the log-likelihood convergence check
  when `conv_criterion = "loglik"` (default 100). Ignored when
  `conv_criterion = "params"`.

- conv_patience:

  Integer. Consecutive iterations that must satisfy the convergence
  criterion before stopping (default 1).

- max_em_iters:

  Maximum number of EM outer iterations per `qrme` call (default 100).

- mcmc_draws:

  total number of MCMC draws per EM step (default 200)

- mcmc_burn_in:

  number of MCMC draws to discard as burnin (default 100)

- proposal_sd:

  Standard deviation of the MH proposal used in the two first-stage
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) calls.
  When `NULL` (default), each first stage uses its own automatic scale,
  `sqrt(var(Y))` or `sqrt(var(T))`.

- start_sigma_y:

  Starting value(s) for the ME standard deviation(s) in the outcome (Y)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call. A
  vector of length `y_n_mix`. When `NULL` (default), `qrme` chooses its
  own starting value.

- start_mu_y:

  Starting value(s) for the ME mean(s) in the outcome (Y)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call. A
  vector of length `y_n_mix`. `NULL` uses the `qrme` default.

- start_pi_y:

  Starting mixture weight(s) for the outcome (Y)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call. A
  vector of length `y_n_mix` that sums to 1. `NULL` uses the `qrme`
  default.

- start_sigma_t:

  Starting value(s) for the ME standard deviation(s) in the
  treatment (T)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call. A
  vector of length `t_n_mix`. `NULL` uses the `qrme` default.

- start_mu_t:

  Starting value(s) for the ME mean(s) in the treatment (T)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call.
  `NULL` uses the `qrme` default.

- start_pi_t:

  Starting mixture weight(s) for the treatment (T)
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) call.
  `NULL` uses the `qrme` default.

- n_copula_draws:

  Number of copula draws per covariate row used to compute mobility
  summaries such as transition matrices and rank correlations (default
  100).

- ignore_me:

  whether or not to ignore measurement error (this is primarily a way to
  get speedy calculations using copula-based approach)

- verbose:

  Logical or nonnegative integer. If `FALSE` (default), suppresses
  progress output. If `TRUE` or `1`, reports major computational stages,
  EM iteration numbers, and convergence diagnostics. If `2`, also
  reports qrme-native details such as `pi`, `mu`, `sigma`, tolerance,
  finite-mixture summaries, and the no-measurement-error comparison
  copula step. If `3` or larger, also reports raw `normalmixEM()`
  output, which uses `lambda` for mixture probabilities.

- se:

  whether or not to estimate standard errors using the boostrap (default
  is FALSE)

- n_boot:

  if computing standard errors, the number of bootstrap iterations to
  use (default is 100)

- n_cores:

  allows for parallel processing in computing standard errors using the
  bootstrap (the default is 1)

- boot_n_copula_me_draws:

  Number of copula ME draws for each bootstrap draw. `NULL` (default)
  uses the same value as `n_copula_me_draws`. Reducing this (e.g.
  to 500) speeds up the bootstrap; per-draw MC noise averages out over
  bootstrap iterations.

- boot_mcmc_draws:

  Total MCMC draws per EM step in each bootstrap draw. `NULL` (default)
  uses the same value as `mcmc_draws`. With warm-starting from the
  full-data fit, fewer draws are needed per bootstrap draw; 200 is a
  reasonable value when the full fit uses 400.

- boot_mcmc_burn_in:

  MCMC burn-in draws per EM step in each bootstrap draw. `NULL`
  (default) uses the same value as `mcmc_burn_in`. Can be reduced
  alongside `boot_mcmc_draws` since warm-starting puts the chain near
  the posterior from the first iteration.

- boot_tol:

  Convergence tolerance for the EM algorithm in each bootstrap draw.
  `NULL` (default) uses the same value as `tol`. A slightly looser
  tolerance (e.g. 5 times the full-fit value) reduces iterations per
  draw; the resulting imprecision averages out over bootstrap
  iterations.

## Value

list of nonlinear measures of intergenerational income mobility adjusted
for measurement error

## Examples

``` r
if (FALSE) { # \dontrun{
  tau      <- seq(0.02, 0.98, length.out = 25)
  t_values <- quantile(nlsy97$lpi, probs = seq(0.1, 0.9, by = 0.1))
  fit <- tsme(
    data            = nlsy97,
    y_formula       = lci ~ ageC_97 + ageF,
    t_formula       = lpi ~ ageC_97 + ageF,
    tau             = tau,
    t_values        = t_values,
    copula          = "frank",
    me_distribution = "laplace",
    mcmc_draws      = 400L,
    mcmc_burn_in    = 20L,
    se              = TRUE,
    n_boot          = 200L,
    n_cores         = 10L
  )
  print(fit)
} # }
```
