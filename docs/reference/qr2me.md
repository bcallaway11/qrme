# QR with 2-sided measurement error

Estimates QR parameters in the case with measurement error in the
outcome and measurement error in a particular continuous "treatment"
variable.

## Usage

``` r
qr2me(
  y_name,
  t_name,
  x_formula,
  tau,
  data,
  x_data = NULL,
  t_values = NULL,
  me_distribution = "gaussian",
  copula = "gaussian",
  Qyx,
  Qtx,
  return_fytx_list = FALSE,
  n_copula_me_draws = 100L,
  n_copula_draws = 100L,
  verbose = FALSE
)
```

## Arguments

- y_name:

  name of the outcome in the passed in data

- t_name:

  name of the treatment in the passed in data

- x_formula:

  a one-sided formula for additional covariates (assumed not to be
  measured with error)

- tau:

  a vector containing particular quantiles that have been estimated ??

- data:

  a data.frame containing the data used for estimation

- x_data:

  If you want conditional distributions to be returned, pass in the
  value of the distribution here; otherwise the default behavior is to
  return a single distribution that averages over all values of X in the
  dataset

- t_values:

  a vector of values of the treatment variable at which to compute the
  conditional distribution of Y given X and T

- me_distribution:

  which type of measurement error distribution to use (default is
  "gaussian"), "laplace" is also supported

- copula:

  which type of copula to use. Options are "gaussian" (the default),
  "clayton", "gumbel", or "frank"

- Qyx:

  quantile regression estimates (can be adjusted for measurement error)
  of Y on X

- Qtx:

  quantile regression estimates (can be adjusted for measurement error)
  of T on X

- return_fytx_list:

  whether or not to return the conditional distribution for every value
  of x in x_data (default is FALSE because this can take up a lot of
  room in memory)

- n_copula_me_draws:

  Number of measurement-error draws used in the second-stage copula
  likelihood (default 100). Increase this for a less noisy copula
  likelihood at the cost of speed.

- n_copula_draws:

  Number of copula draws per covariate row used to compute mobility
  summaries such as transition matrices and rank correlations (default
  100).

- verbose:

  Logical or nonnegative integer. If `FALSE` (default), suppresses
  progress output. If `TRUE` or `1`, reports major computational stages,
  EM iteration numbers, and convergence diagnostics. If `2`, also
  reports qrme-native details such as `pi`, `mu`, `sigma`, tolerance,
  finite-mixture summaries, and the no-measurement-error comparison
  copula step. If `3` or larger, also reports raw `normalmixEM()`
  output, which uses `lambda` for mixture probabilities.
