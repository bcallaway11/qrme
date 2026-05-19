# fvC

Computes the density of the measurement error term given it follows a
mixture of normals distribution

## Usage

``` r
fvC(v, me_distribution, m, pi, mu, sig)
```

## Arguments

- v:

  value to estimate the density at

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

## Value

estimated density of measurement error at v
