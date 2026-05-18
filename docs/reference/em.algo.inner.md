# Inner part of EM-algorithm for QR with measurement error

Does the heavy-lifting of the EM-algorithm for QR with measurment error

## Usage

``` r
em.algo.inner(
  formula,
  data,
  betmat,
  tau,
  me_distribution = "gaussian",
  m = 1,
  pi = 1,
  mu = 0,
  sig = 1,
  mcmc_method = "MH",
  mcmc_draws = 200,
  mcmc_burn_in = 100,
  proposal_sd = NULL,
  n_cores = 1,
  verbose = FALSE
)
```

## Arguments

- formula:

  y ~ x

- data:

  a data.frame that contains y and x

- betmat:

  LxK matrix of parameter values with L the number of quantiles and K
  the dimension of the covariates

- tau:

  An L-vector indicating which quantiles have been estimated

- me_distribution:

  The distribution of the measurement error. "gaussian" is the default
  and supports a mixture of normals. "laplace" is also supported.

- m:

  The number of components of the mixture distribution for the
  measurement error

- pi:

  The probability of each mixture component (should have length equal to
  m)

- mu:

  The mean of each mixture component (should have length equal to m)

- sig:

  The standard deviation of each mixture component (should have length
  equal to m)

- mcmc_method:

  Simulation step to use in the EM algorithm. Currently only `"MH"` is
  supported.

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

- verbose:

  Logical or nonnegative integer. If `FALSE` (default), suppresses
  progress output. If `TRUE` or `1`, reports major computational stages,
  EM iteration numbers, and convergence diagnostics. If `2`, also
  reports qrme-native details such as `pi`, `mu`, `sigma`, tolerance,
  and finite-mixture summaries. If `3` or larger, also reports raw
  `normalmixEM()` output, which uses `lambda` for mixture probabilities.

## Value

A list of QR parameters and parameters for mixture of normals for the
measurement error term
