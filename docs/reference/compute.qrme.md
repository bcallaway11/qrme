# compute.qrme

does the heavy lifting on computing quantile regression with left hand
side measurement error

## Usage

``` r
compute.qrme(
  formula,
  tau = 0.5,
  data,
  me_distribution = "gaussian",
  n_mix = 1,
  start_beta = NULL,
  start_mu = NULL,
  start_sigma = NULL,
  start_pi = NULL,
  mcmc_method = "MH",
  tol = NULL,
  conv_criterion = "params",
  loglik_draws = 100L,
  conv_patience = 1L,
  mcmc_draws = 200,
  mcmc_burn_in = 100,
  proposal_sd = NULL,
  n_cores = 1,
  max_em_iters = 100,
  verbose = FALSE
)
```

## Arguments

- formula:

  y ~ x

- tau:

  vector for which quantiles to compute quantile regression

- data:

  a data.frame that contains y and x

- me_distribution:

  Distribution of the measurement error: `"gaussian"` (default, supports
  a mixture of normals) or `"laplace"`.

- n_mix:

  The number of mixture components of the measurement error (default 1).
  Increase this for more flexible measurement-error distributions.

- start_beta:

  an LxK matrix of starting values for beta where L is the dimension of
  tau and K is the number of covariates (default is NULL and in this
  case, the starting values are set to be the QR coefficients coming
  from QR that ignores measurment error)

- start_mu:

  A vector of length n_mix of starting values for the mean of each
  mixture component. When `NULL` (default), set to 0 for `n_mix = 1` and
  to evenly-spaced values on `[-0.25 * sd(y), 0.25 * sd(y)]` for
  `n_mix > 1`, mean-centered to satisfy the mean-zero ME constraint.

- start_sigma:

  A vector of length n_mix of starting values for the standard deviation
  of each mixture component. When `NULL` (default), backed out from the
  mixture variance identity assuming equal weights and equal component
  SDs targeting a mixture SD of `0.5 * sd(y)`:
  `sigma = sqrt(max(target^2 - mean(start_mu^2), (target/sqrt(n_mix))^2))`.
  This gives `0.5 * sd(y)` for `n_mix = 1` and slightly larger values
  for `n_mix > 1` (since the component means absorb some variance).

- start_pi:

  A vector of length n_mix of starting values for the fraction of
  observations in each component of the mixture of normals distribution
  for the measurement error (default is NULL and in this case, the
  starting values are all set to be 1/n_mix)

- mcmc_method:

  The type of simulation step to use in the EM algorithm. Currently only
  `"MH"` for Metropolis-Hastings is supported.

- tol:

  Convergence tolerance. When `NULL` (default), a value is chosen
  automatically based on `conv_criterion`. See
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md)
  for details.

- conv_criterion:

  Convergence criterion: `"params"` (default) uses the Euclidean norm of
  parameter changes; `"loglik"` uses the relative change in the
  observed-data log-likelihood. See
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md).

- loglik_draws:

  Number of Monte Carlo draws for the log-likelihood convergence check
  when `conv_criterion = "loglik"` (default 100). Ignored when
  `conv_criterion = "params"`.

- conv_patience:

  Integer. Consecutive iterations that must satisfy the convergence
  criterion before stopping (default 1). Setting to 2 guards against
  false convergence from Monte Carlo noise. See
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md).

- mcmc_draws:

  Total number of MCMC draws per EM step (default 200)

- mcmc_burn_in:

  Number of MCMC draws to discard as burnin (default 100)

- proposal_sd:

  Standard deviation of the Metropolis-Hastings proposal (random-walk
  step size). When `NULL` (default), initialised to the SD of the ME
  mixture from the starting parameters and updated after each M-step to
  track the current ME scale. Pass a positive numeric to fix the
  proposal at that value for all iterations.

- n_cores:

  Number of cores for parallel bootstrap computation (default 1)

- max_em_iters:

  Maximum number of EM outer iterations. If convergence is not reached,
  the estimates from the final iteration are returned (default is 100)

- verbose:

  Logical or nonnegative integer. If `FALSE` (default), suppresses
  progress output. If `TRUE` or `1`, reports major computational stages,
  EM iteration numbers, and convergence diagnostics. If `2`, also
  reports qrme-native details such as `pi`, `mu`, `sigma`, tolerance,
  and finite-mixture summaries. If `3` or larger, also reports raw
  `normalmixEM()` output, which uses `lambda` for mixture probabilities.
