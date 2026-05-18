# parms2coppar

convert parameters from optimization routine into copula parameters,
parameters are between 0 and 1

## Usage

``` r
parms2coppar(params, copula = "gumbel", x)
```

## Arguments

- params:

  a vector of parameters

- copula:

  the type of copula, should be supported by cop.pdf

- x:

  a matrix of values of x

## Value

a vector of copula parameters
