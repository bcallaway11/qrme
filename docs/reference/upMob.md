# upMob

Upward Mobility Paramters

## Usage

``` r
upMob(Y, T, amount = 0, qcutoffs = c(0.25, 0.5, 0.75))
```

## Arguments

- Y:

  vector of outcomes

- T:

  vector of treatments

- amount:

  How much child's rank has to increase in order to count for upward
  mobility

- qcutoffs:

  the quantiles where to draw the cutoff points in the transition
  matrix; default is at the quartiles

## Value

vector of amount of upward mobility by parents' income quartile (or
other specified cutoffs)
