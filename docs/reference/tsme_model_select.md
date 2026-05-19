# tsme_model_select

Fit [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md) over
a grid of copula families and measurement error distributions and rank
the specifications by AIC and BIC.

The log-likelihood is computed as the sum of three components:

- `ll_y`: log-likelihood of the outcome ME model (`logLik(fit$me_qyx)`)

- `ll_t`: log-likelihood of the treatment ME model
  (`logLik(fit$me_qtx)`)

- `ll_cop`: copula log-likelihood (`fit$me_cop_loglik`)

\$\$AIC = -2 \ell + 2k, \quad BIC = -2 \ell + k \log n\$\$ where \\k\\
is the number of free parameters specified by `k_params`. For
`y_n_mix = t_n_mix = 1` the free parameters are \\\sigma_Y\\,
\\\sigma_T\\, and \\\theta\_{\text{cop}}\\, giving `k_params = 3`.

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
  k_params = 3L,
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

- k_params:

  number of free parameters used in AIC/BIC (default 3, appropriate for
  `y_n_mix = t_n_mix = 1`)

- return_fits:

  logical; if `TRUE` the fitted tsme objects are returned alongside the
  table (default `FALSE`)

- ...:

  additional arguments passed to
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md) for
  every cell in the grid (e.g. `n_copula_me_draws`, `mcmc_draws`,
  `y_n_mix`, `t_n_mix`)

## Value

a list with element `table` (a data.frame sorted by BIC) and, if
`return_fits = TRUE`, element `fits` (a list of tsme objects in the same
order as `table`)

## Examples

``` r
if (FALSE) { # \dontrun{
  tau      <- seq(0.05, 0.95, length.out = 15)
  t_values <- quantile(nlsy97$lpi, probs = seq(0.1, 0.9, by = 0.1))
  sel <- tsme_model_select(
    data            = nlsy97,
    y_formula       = lci ~ ageC_97 + ageF,
    t_formula       = lpi ~ ageC_97 + ageF,
    tau             = tau,
    t_values        = t_values,
    me_distribution = "laplace",
    mcmc_draws      = 200L,
    mcmc_burn_in    = 20L,
    n_cores         = 4L
  )
  print(sel)
} # }
```
