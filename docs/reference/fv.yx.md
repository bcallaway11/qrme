# fv.yx

Return the density of the measurement error conditional on y and x; this
takes as given some QR parameters from Y\* (the true outcome)
conditional on X. Here, we also presume that the distribution of the
measurement error is a mixture of normal distributions

## Usage

``` r
fv.yx(v, betmat, m, pi, mu, sig, Y, X, tau)
```

## Arguments

- v:

  A particular value of the measurement error to estimate f(v\|y,x)

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

- Y:

  An nx1 vector of outcomes

- X:

  An nxK matrix of covariates

- tau:

  an L-vector of all the quantiles where betas were estimated

## Value

n x 1 vector of f(v\|Y,X)
