# fy.x

return the density of Y (for particular value y) conditional on X (which
can include n observations) when Q(Y\|X) has been estimated using QR.
This is used in our simulation approach.

## Usage

``` r
fy.x(y, betmat, X, tau)
```

## Arguments

- y:

  particular value of y to estimate f(y\|x)

- betmat:

  LxK matrix of parameter values with L the number of quantiles and K
  the dimension of the covariates

- X:

  An nxK matrix with n the number of observations of X

- tau:

  an L-vector containing the quantile at which Q(Y\|X) was estimated

## Value

An nx1 vector that contains f(y\|X)
