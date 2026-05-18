# tmat

A function for computing transition matrices

## Usage

``` r
tmat(Y, T, qcutoffs = c(0.25, 0.5, 0.75))
```

## Arguments

- Y:

  vector of outcomes

- T:

  vector of treatments

- qcutoffs:

  the quantiles where to draw the cutoff points in the transition
  matrix; default is at the quartiles

## Value

A matrix containing values of the transition matrix. The columns contain
different groups by value of the treatment, the rows contain different
values for the outcome (it is ordered from bottom to top though as is
standard for a transition matrix)
