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
