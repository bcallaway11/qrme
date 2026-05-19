# fYXmatC

Takes n observations of X and a vector of Y's and return an n x Y.size()
matrix of fyx evaluated at those points

## Usage

``` r
fYXmatC(Y, betmat, X, tau)
```

## Arguments

- Y:

  vector of outcomes

- betmat:

  matrix of QR parameters

- X:

  matrix of covariates

- tau:

  which values QR's have been estimated for

## Value

matrix of conditional density estimates
