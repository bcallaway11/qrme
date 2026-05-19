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

## Examples

``` r
tmat(nlsy97$lci, nlsy97$lpi)
#>           [,1]      [,2]      [,3]      [,4]
#> [1,] 0.1605839 0.2335766 0.2554745 0.3759124
#> [2,] 0.1970803 0.2153285 0.3029197 0.2664234
#> [3,] 0.3284672 0.2299270 0.3138686 0.1751825
#> [4,] 0.3029197 0.2773723 0.1751825 0.1897810
```
