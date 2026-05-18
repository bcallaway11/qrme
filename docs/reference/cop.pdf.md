# cop.pdf

returns the pdf of several copula functions evaluated at u and v; the
available copula functions are: Gumbel (the default), Frank, (should
implement Joe and Gaussian)

## Usage

``` r
cop.pdf(u, v, type = "gumbel", delt, eps = .Machine$double.eps)
```

## Arguments

- u:

  first copula argument (u and v can be vectors)

- v:

  second copula argument

- type:

  one of c("gumbel","frank")

- delt:

  the copula paramter

- eps:

  A small positive number to trim out zeros

## Value

vector of copula pdf values
