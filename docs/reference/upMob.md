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

## Examples

``` r
upMob(nlsy97$lci, nlsy97$lpi)
#> [1] 0.8302583 0.5725191 0.3867596 0.1666667
```
