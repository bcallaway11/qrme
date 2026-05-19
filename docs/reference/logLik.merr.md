# logLik.merr

Observed-data log-likelihood for a `merr` object via Monte Carlo
integration over the measurement error distribution. Enables `AIC(fit)`
and `BIC(fit)` automatically.

## Usage

``` r
# S3 method for class 'merr'
logLik(object, ndraws = 500, ...)
```

## Arguments

- object:

  A `merr` object returned by
  [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md).

- ndraws:

  Number of Monte Carlo draws for the integral (default 500).

- ...:

  Currently unused.

## Value

A scalar of class `"logLik"` with attributes `df` (number of free
parameters) and `nobs` (sample size).

## Examples

``` r
# \donttest{
  set.seed(42)
  n <- 200; X <- runif(n)
  Y <- 1 + 2 * X + rnorm(n) + rnorm(n, sd = 0.5)
  fit1 <- qrme(Y ~ X, data.frame(Y, X), tau = 0.5,
               n_mix = 1L, mcmc_draws = 80L, mcmc_burn_in = 40L,
               max_em_iters = 10L, verbose = FALSE)
  fit2 <- qrme(Y ~ X, data.frame(Y, X), tau = 0.5,
               n_mix = 2L, mcmc_draws = 80L, mcmc_burn_in = 40L,
               max_em_iters = 10L, verbose = FALSE)
#> Warning: Too many fixups:  doubling m
#> Warning: EM algorithm failed to converge after 10 iterations
  data.frame(n_mix = 1:2,
             AIC   = c(AIC(fit1), AIC(fit2)),
             BIC   = c(BIC(fit1), BIC(fit2)))
#>   n_mix      AIC      BIC
#> 1     1 745.0586 756.4428
#> 2     2 746.3905 766.1827
# }
```
