# fv

The distribution of the measurement error using a mixture of normal
distributions

## Usage

``` r
fv(v, m = 1, pi = 1, mu = 0, sig = 1)
```

## Arguments

- v:

  A particular value of the measurement error to estimate f(v\|y,x)

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

## Value

scalar f(v)
