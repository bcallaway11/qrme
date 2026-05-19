# imp_sampC

return a vector of weights to be used in importance sampling note that,
unlike mh_mcmcC, here the measurement error vector has already been
drawn and all we need to do is compute weights

## Usage

``` r
imp_sampC(
  Y,
  X,
  V,
  mcmc_draws,
  proposal_sd,
  betmat,
  me_distribution,
  m,
  pi,
  mu,
  sig,
  tau
)
```

## Arguments

- Y:

  vector of outcomes

- X:

  matrix of covariates

- V:

  vector of measurement errors

- mcmc_draws:

  total number of Monte Carlo draws

- proposal_sd:

  standard deviation of the random-walk MH proposal MH algorithm

- betmat:

  matrix of QR parameters

- me_distribution:

  the distribution of the measurement error. "gaussian" is the default
  and supports a mixture of normals. "laplace" is also supported.

- m:

  number of mixture components for measurement error

- pi:

  mixture probabilities

- mu:

  means of mixture components

- sig:

  standard deviations of mixture components

- tau:

  which values QR's have been estimated for

## Value

vector of weights to be used in importance sampling
