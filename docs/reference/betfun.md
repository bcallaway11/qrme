# betfun

Turns matrices of QR parameter values into functions. The returned
function will be a map from (0,1) to a vector of dimension K where K is
the number of elements in X (also it is the number of columns in the
passed in matrix of QR parameter values). It works by linear
interpolation of parameter values (for interior values of tau) and based
on some assumptions to handle values in the tails.

## Usage

``` r
betfun(betmat, tau)
```

## Arguments

- betmat:

  An L x K matrix of parameter values where L is the number of tau and K
  is the dimension of X

- tau:

  The particular quantiles that the parameters were estimated at. It
  should be an L-dimensional vector.

## Value

A function that can be called on any element from 0 to 1 and returns a
vector of parameter values
