# rlaplace

Draw samples from the Laplace distribution. Note that the function takes
the standard deviation as an argument, but it converts it to the scale
of the Laplace distribution, which is the typical parameterization.

## Usage

``` r
rlaplace(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  number of samples to draw

- mu:

  location parameter (default is 0)

- sigma:

  standard deviation, function will convert it to the scale

## Value

vector of n samples from the Laplace distribution

## Examples

``` r
set.seed(1)
x <- rlaplace(500, mu = 0, sigma = 1)
c(mean = mean(x), sd = sd(x))
#>       mean         sd 
#> 0.02228648 0.93280121 
```
