# mh_mcmcC

mh_mcmcC

## Usage

``` r
mh_mcmcC(
  Y,
  X,
  startval,
  mcmc_draws,
  mcmc_burn_in,
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

- startval:

  starting value for the markov chain

- mcmc_draws:

  total number of Monte Carlo draws

- mcmc_burn_in:

  number of initial draws to discard as burnin

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
