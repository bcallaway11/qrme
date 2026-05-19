# fvyxC

Computes density of measurement error conditional on y and x given QR
estimates and distribution of measurement error

## Usage

``` r
fvyxC(v, betmat, me_distribution, m, pi, mu, sig, y, x, tau)
```

## Arguments

- v:

  value to estimate the density at

- betmat:

  matrix of quantile regression coefficients

- me_distribution:

  the distribution of the measurement error. "gaussian" is the default
  and supports a mixture of normals. "laplace" is also supported.

- m:

  number of mixture components

- pi:

  vector of mixture probabilities

- mu:

  vector of mixture means

- sig:

  vector of mixture standard deviations

- y:

  a particular value of y

- x:

  a vector of x's

- tau:

  vector of values where QR were estimated

## Value

estimate of density of measurement error conditional on y and x
