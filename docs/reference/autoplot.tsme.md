# autoplot.tsme

ggplot2-based plot method for tsme objects. Returns a `ggplot` object
that can be further customized.

## Usage

``` r
# S3 method for class 'tsme'
autoplot(
  object,
  type = c("cond_quant", "pov_rate"),
  which = c("me", "nome", "qr", "all"),
  ci = TRUE,
  level = 0.9,
  ...
)
```

## Arguments

- object:

  a tsme object returned by
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md)

- type:

  plot type: `"cond_quant"` (conditional quantile curves, default) or
  `"pov_rate"` (poverty rates)

- which:

  which estimator to display: `"me"` (ME-corrected, default), `"nome"`
  (no-ME), `"qr"` (standard QR), or `"all"` (all three on one plot; only
  valid when `type = "pov_rate"`)

- ci:

  logical; add confidence bands when bootstrap SEs are available
  (default `TRUE`)

- level:

  confidence level for the bands (default 0.90)

- ...:

  unused

## Value

a ggplot object

## Examples

``` r
# \donttest{
  ggplot2::autoplot(nlsy97_tsme_fit, type = "cond_quant", which = "me")

  ggplot2::autoplot(nlsy97_tsme_fit, type = "pov_rate",   which = "all", ci = FALSE)

# }
```
