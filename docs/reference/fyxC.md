# fyxC

Converts quantile regression estimates into density estimates

## Usage

``` r
fyxC(y, betmat, X, tau)
```

## Arguments

- y:

  vector of outcomes

- betmat:

  matrix of quantile regression coefficients

- X:

  vector of covariates (should be same dimension as betamat)

- tau:

  vector of values where QR were estimated

## Value

value of density at y and X given QR estimates
