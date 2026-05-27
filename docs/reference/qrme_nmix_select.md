# qrme_nmix_select

Fit [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) with
a range of mixture component counts and rank the specifications by AIC
and BIC. Use this to choose the number of ME mixture components for a
single equation.

The number of free ME parameters is:

- `n_mix = 0`: \\k = 0\\ (no ME; LL from standard QR)

- `n_mix = 1`: \\k = 2\\ (\\\mu\\, \\\sigma\\)

- `n_mix >= 2`: \\k = 3m - 1\\ (\\\mu_1,\ldots,\mu_m\\,
  \\\sigma_1,\ldots,\sigma_m\\, \\\pi_1,\ldots,\pi\_{m-1}\\)

## Usage

``` r
qrme_nmix_select(
  formula,
  data,
  tau,
  n_mix_vals = 0:3,
  return_fits = FALSE,
  ...
)
```

## Arguments

- formula:

  formula for the outcome model (passed to
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md))

- data:

  data.frame

- tau:

  quantile levels for QR

- n_mix_vals:

  integer vector of component counts to evaluate (default `0:3`; `0` is
  the no-ME baseline)

- return_fits:

  logical; if `TRUE` the fitted qrme objects are included in the output
  (default `FALSE`)

- ...:

  additional arguments passed to
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) for
  every fit (e.g. `me_distribution`, `mcmc_draws`, `mcmc_burn_in`).
  `n_mix` and `se` are set internally and must not be passed.

## Value

a list with:

- `table`:

  data.frame sorted by BIC with columns `n_mix`, `k_me`, `ll`, `AIC`,
  `BIC`

- `fits`:

  (only when `return_fits = TRUE`) named list of qrme objects; the
  `n_mix = 0` entry is `NULL`

## Examples

``` r
if (FALSE) { # \dontrun{
  set.seed(1)
  n <- 300; X <- runif(n)
  Y <- 1 + 2 * X + rnorm(n) + rnorm(n, sd = 0.5)
  sel <- qrme_nmix_select(
    Y ~ X, data.frame(Y, X),
    tau          = c(0.25, 0.5, 0.75),
    n_mix_vals   = 0:3,
    mcmc_draws   = 100L,
    mcmc_burn_in = 50L
  )
  sel$table
} # }
```
