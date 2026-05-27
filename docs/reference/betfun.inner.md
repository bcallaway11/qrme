# betfun.inner

Does the heavy lifting for betfun. Basically, betfun is just a wrapper
for this that can handle a matrix of values of parameters. This function
does the work, but only for a single vector of betas.

## Usage

``` r
betfun.inner(betvec, tau, isconst = FALSE)
```

## Arguments

- betvec:

  vector of parameter values

- tau:

  corresponding vector of quantiles where beta was estimated

- isconst:

  logical; if `TRUE` the function is treated as a constant (intercept)
  term and tail extrapolation uses a log-linear correction rather than a
  flat extension

## Value

function that takes argument from (0,1)
