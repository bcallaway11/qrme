# qrme_start_search

Run [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) from
multiple random starting values and return the fit with the highest
log-likelihood. Use this to assess whether the EM algorithm converges to
a unique solution or is sensitive to starting values.

For each restart, `start_sigma` is drawn uniformly from `sigma_range`,
`start_mu` is drawn from \\N(0, 0.01)\\, and mixture weights `start_pi`
are drawn from a symmetric Dirichlet (i.e., uniform on the simplex).

## Usage

``` r
qrme_start_search(
  formula,
  data,
  tau,
  n_starts = 10L,
  seed = NULL,
  sigma_range = c(0.1, 2),
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

- n_starts:

  number of random restarts (default 10)

- seed:

  optional integer seed for reproducibility (default `NULL`)

- sigma_range:

  length-2 numeric vector giving the \\\[\min, \max\]\\ interval from
  which `start_sigma` is drawn for each component (default
  `c(0.1, 2.0)`)

- return_fits:

  logical; if `TRUE` the fitted qrme objects for all successful restarts
  are returned (default `FALSE`)

- ...:

  additional arguments passed to
  [`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md) for
  every restart (e.g. `n_mix`, `me_distribution`, `mcmc_draws`).
  Arguments `start_sigma`, `start_mu`, `start_pi`, and `se` are set
  internally and should not be passed via `...`.

## Value

a list with:

- `table`:

  data.frame of start index and log-likelihood for each successful
  restart, sorted descending by log-likelihood

- `best_fit`:

  the qrme fit with the highest log-likelihood

- `fits`:

  (only when `return_fits = TRUE`) list of all successful qrme fits
