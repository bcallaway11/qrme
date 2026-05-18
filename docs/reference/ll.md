# ll

this is the likelihood function for estimating the copula parameter

## Usage

``` r
ll(
  params,
  y,
  t,
  x,
  copula = "gumbel",
  Fyx,
  Ftx,
  fyx,
  ftx,
  Us,
  Vs,
  ndraws = 100,
  eps = .Machine$double.eps
)
```

## Arguments

- params:

  a value of the parameters to calculate the log likelihood for

- y:

  should be a vector of outcomes

- t:

  should be a vector of treatments

- x:

  a matrix of values of x

- copula:

  the type of copula, should be supported by cop.pdf

- Fyx:

  function that contains the cdf of y conditional on x

- Ftx:

  function that contains the cdf of t conditional on x

- fyx:

  function that contains the pdf of y conditional on x

- ftx:

  function that contains the pdf of t conditional on x

- Us:

  draws for the measurement error in the outcome equation (these are
  made in the function that calls the likelihood function)

- Vs:

  draws for the measurement error in the treatment equation (these are
  made in the function that calls the likelihood function)

## Value

scalar negative value of log likelihood function
