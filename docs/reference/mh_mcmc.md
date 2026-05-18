# mh_mcmc

A Metropolis-Hastings algorithm for drawing measurment errors.

## Usage

``` r
mh_mcmc(
  startval = 0,
  mcmc_draws = 200,
  mcmc_burn_in = 100,
  proposal_sd = NULL,
  betmat,
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

  The first value in the markov chain

- mcmc_draws:

  The total number of measurement error draws to make

- mcmc_burn_in:

  The number of draws to drop

- proposal_sd:

  Standard deviation of the random-walk MH proposal. Passed from
  `em.algo` after automatic scaling; see
  [`em.algo`](https://bcallaway11.github.io/qrme/reference/em.algo.md).

- betmat:

  LxK matrix of parameter values with L the number of quantiles and K
  the dimension of the covariates

- m:

  The dimension of the measurement error

- pi:

  The probability of each mixture component (should have length equal to
  m)

- mu:

  The mean of each mixture component (should have length equal to m)

- sig:

  The standard deviation of each mixture component (should have length
  equal to m)

- y:

  particular value of y

- x:

  particular value of x

- tau:

  an L-vector of all the quantiles where betas were estimated

## Value

vector of draws of measurement error
