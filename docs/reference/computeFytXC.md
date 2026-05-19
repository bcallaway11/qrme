# computeFytXC

convert first step quantile regression estimates and second step copula
estimates into the distribution of Y conditional on T and X this is not
currently used (we replaced this numerical approach with a direct
calculation inside the R code). This also only works currently for a
Gumbel copula.

## Usage

``` r
computeFytXC(yvals, t_values, Qyxpreds, Ftxpreds, tau, copula, copParam)
```

## Arguments

- yvals:

  values to compute conditional distribution for

- t_values:

  values of treatment to compute conditional distribution for

- Qyxpreds:

  matrix of quantile regression predictions

- Ftxpreds:

  matrix of conditional distribution (T conditional on X) predictions

- tau:

  values at which QR's were estimated

- copula:

  which type of copula was estimated (currently ignored and imposes that
  the copula is Gumbel

- copParam:

  previously estimated copula parameter

## Value

cube corresponding to conditional distribution as a function of y, t,
and each row of X in the data
