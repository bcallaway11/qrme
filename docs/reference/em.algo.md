# em.algo

A pseudo EM algorithm for quantile regression with measurement error.
The measurement error here follows a mixture of normals.

## Usage

``` r
em.algo(
  formula,
  data,
  betmatguess,
  tau,
  me_distribution = "gaussian",
  m = 1,
  piguess = 1,
  muguess = 0,
  sigguess = 1,
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

- data:

  a data.frame that contains y and x

- betmatguess:

  Initial values for the beta parameters. This should be an LxK matrix
  where L is the number of quantiles and K is the dimension of the
  covariates

- tau:

  An L-vector indicating which quantiles have been estimated

- me_distribution:

  The distribution of the measurement error. "gaussian" is the default
  and supports a mixture of normals. "laplace" is also supported.

- m:

  The number of components of the mixture distribution for the
  measurement error

- piguess:

  Starting value for the probabilities of each mixture distribution
  (should have length equal to k)

- muguess:

  Starting value for the mean of each mixture component (should have
  length equal to k)

- sigguess:

  Starting value for the standard deviation of each mixture component
  (should have length equal to k)

- mcmc_method:

  Simulation step to use in the EM algorithm. Currently only `"MH"` is
  supported.

- tol:

  Convergence tolerance. When `NULL` (default), `1e-2` is used for both
  criteria. For `"params"`, `tol` is compared against
  `max(|theta_new - theta_old| / (1 + |theta_old|))`, the maximum scaled
  parameter change across all parameters; this is dimensionless and
  application-invariant. For `"loglik"`, it is compared against the
  relative log-likelihood change `|ll_new - ll_old| / |ll_old|`. Tighten
  to `1e-3` or smaller for higher accuracy at the cost of more
  iterations.

- conv_criterion:

  Convergence criterion. `"params"` (default) stops when the maximum
  scaled parameter change falls below `tol` (see `tol` above).
  `"loglik"` stops when the relative change in the observed-data
  log-likelihood, `|ll_new - ll_old| / |ll_old|`, falls below `tol`; the
  log-likelihood is evaluated using Monte Carlo integration with
  `loglik_draws` draws per outer iteration.

- loglik_draws:

  Number of Monte Carlo draws used to evaluate the log-likelihood at
  each outer EM iteration when `conv_criterion = "loglik"`. More draws
  reduce noise in the convergence check at the cost of speed. Default is
  100.

- conv_patience:

  Integer. Number of consecutive iterations that must satisfy the
  convergence criterion before the algorithm stops. Default is 1
  (current behaviour). Setting this to 2 requires two back-to-back
  iterations below `tol`, guarding against false convergence caused by
  Monte Carlo noise in the log-likelihood or parameter updates.

- mcmc_draws:

  Total number of MCMC draws per EM step (default 200)

- mcmc_burn_in:

  Number of MCMC draws to discard as burnin (default 100)

- proposal_sd:

  Standard deviation of the Metropolis-Hastings proposal (random-walk
  step size). When `NULL` (default), initialised to the SD of the ME
  mixture from the starting parameters and updated after each M-step to
  track the current ME scale: `sqrt(sum(pi_k * (sig_k^2 + mu_k^2)))`.
  Pass a positive numeric to fix the proposal at that value for all
  iterations.

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

## Value

QRME object
