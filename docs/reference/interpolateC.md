# interpolateC

Returns interpolated value at x from parallel arrays ( xData, yData )
Assumes that xData has at least two elements, is sorted and is strictly
monotonic increasing

## Usage

``` r
interpolateC(x, y, xval, extrapolate)
```

## Arguments

- x:

  vector of x's

- y:

  vector of y's

- xval:

  a particular value to interpolate for

- extrapolate:

  whether or not to linearly extrapolate beyond endpoints of x

## Value

interpolated value
