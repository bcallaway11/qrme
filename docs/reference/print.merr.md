# print.merr

Print method for merr objects returned by
[`qrme`](https://bcallaway11.github.io/qrme/reference/qrme.md). Shows
the measurement error distribution type and its fitted parameters
(mixture weights, means, and standard deviations) with bootstrap
standard errors when available.

## Usage

``` r
# S3 method for class 'merr'
print(x, ...)
```

## Arguments

- x:

  a merr object

- ...:

  unused

## Examples

``` r
# \donttest{
  set.seed(1)
  n <- 200; X <- runif(n)
  Y <- 1 + X + rnorm(n) + rnorm(n, sd = 0.5)
  fit <- qrme(Y ~ X, data.frame(Y, X), tau = 0.5,
              mcmc_draws = 80L, mcmc_burn_in = 40L,
              max_em_iters = 10L, verbose = FALSE)
#> Warning: Too many fixups:  doubling m
#> Warning: Too many fixups:  doubling m
#> Warning: Too many fixups:  doubling m
#> Warning: EM algorithm failed to converge after 10 iterations
  print(fit)
#> Measurement error: gaussian
#>  Pi Mu Sigma
#>   1  0 0.443
# }
```
