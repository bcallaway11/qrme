# tsme_model_select

Fit [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md) over
a grid of copula families, measurement error distributions, and mixture
component counts, then rank the specifications by AIC and BIC.

The log-likelihood is computed as the sum of three components:

- `ll_y`: log-likelihood of the outcome ME model

- `ll_t`: log-likelihood of the treatment ME model

- `ll_cop`: copula log-likelihood

\$\$AIC = -2 \ell + 2k, \quad BIC = -2 \ell + k \log n\$\$ where \\k =
k_Y + k_T + 1\\ is computed per grid row, with \\k_m = 2\\ for \\m = 1\\
and \\k_m = 3m - 1\\ for \\m \ge 2\\, plus one copula parameter.

Laplace ME does not support mixture models (`n_mix > 1`); any grid row
combining Laplace with `y_n_mix > 1` or `t_n_mix > 1` is dropped
automatically with a message.

## Usage

``` r
tsme_model_select(
  data,
  y_formula,
  t_formula,
  tau,
  t_values,
  copulas = c("gaussian", "clayton", "gumbel", "frank"),
  me_distributions = c("gaussian", "laplace"),
  y_n_mix = 1L,
  t_n_mix = 1L,
  n_cores = 1L,
  return_fits = FALSE,
  ...
)
```

## Arguments

- data:

  data.frame passed to
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md)

- y_formula:

  formula for the outcome model

- t_formula:

  formula for the treatment model

- tau:

  quantile levels for the first-stage QR (passed to
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md))

- t_values:

  treatment values for conditional distribution summaries

- copulas:

  character vector of copula families to try. Any subset of
  `c("gaussian","clayton","gumbel","frank")` (default: all four).

- me_distributions:

  character vector of ME distributions to try. Any subset of
  `c("gaussian","laplace")` (default: both).

- y_n_mix:

  integer or integer vector of ME mixture component counts to try for
  the outcome equation (default `1L`). When a vector is supplied, all
  values are crossed with the other grid dimensions.

- t_n_mix:

  integer or integer vector of ME mixture component counts to try for
  the treatment equation (default `1L`).

- n_cores:

  number of parallel workers (default 1)

- return_fits:

  logical; if `TRUE` the fitted tsme objects are returned alongside the
  table (default `FALSE`)

- ...:

  additional arguments passed to
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md) for
  every cell in the grid (e.g. `n_copula_me_draws`, `mcmc_draws`)

## Value

a list with element `table` (a data.frame sorted by BIC with columns
`copula`, `me_distribution`, `y_n_mix`, `t_n_mix`, `k_params`, `ll_y`,
`ll_t`, `ll_cop`, `ll`, `aic`, `bic`) and, if `return_fits = TRUE`,
element `fits` (a list of tsme objects in the same order as `table`)

## Examples

``` r
if (FALSE) { # \dontrun{
  tau      <- seq(0.05, 0.95, length.out = 15)
  t_values <- quantile(nlsy97$lpi, probs = seq(0.1, 0.9, by = 0.1))
  sel <- tsme_model_select(
    data             = nlsy97,
    y_formula        = lci ~ ageC_97 + ageF,
    t_formula        = lpi ~ ageC_97 + ageF,
    tau              = tau,
    t_values         = t_values,
    y_n_mix          = 1:2,
    t_n_mix          = 1:2,
    n_cores          = 4L,
    mcmc_draws       = 200L,
    mcmc_burn_in     = 20L
  )
  print(sel)
} # }
```
