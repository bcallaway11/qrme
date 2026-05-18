# makeRQS

Take the results of the optimization and convert them into a quantile
regression object, so we can use all the tools from the quantreg package
(e.g. inverting the quantiles). The key step here is rearrangement,
because the objective function doesn't impose any ordering – see the
discussion in HLLP. We follow HLLP's recommendation and order the QR
parameters by what makes the quantiles to be increasing for the mean
values of the x's. This means that for any particular value of the x's,
the quantile function may not necessarily be increasing in tau. However,
we can separately rearrange those as needed. But this gives a way to
report the QR parameters.

## Usage

``` r
makeRQS(params, formla, data, tau)
```

## Arguments

- params:

  an LxK matrix of QR parameters where L is the number of quantiles that
  parameters have been estimated at and K is the dimension of the
  covariates.

- formla:

  y ~ x, a formula for the outcome on the regressors

- data:

  a data.frame containing y and x

- tau:

  the vector of quantiles where QR was estimated

## Value

rqs object
