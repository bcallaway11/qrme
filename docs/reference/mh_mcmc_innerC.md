# mh_mcmc_innerC

Inner part of MCMC algorithm

## Usage

``` r
mh_mcmc_innerC(
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
  y,
  x,
  tau
)
```

## Arguments

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

- y:

  outcome value (for particular observation)

- x:

  vector of covariates (for particular observation)

- tau:

  which values QR's have been estimated for

## Value

vector of MCMC draws of measurement error
