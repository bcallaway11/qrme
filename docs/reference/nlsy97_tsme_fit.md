# Pre-computed tsme() fit for the NLSY97 intergenerational mobility application

A trimmed \`tsme\` object from fitting the two-sided measurement error
model to the
[`nlsy97`](https://bcallaway11.github.io/qrme/reference/nlsy97.md) data,
as reported in Callaway, Li, Murtazashvili, and Tsyawo (2026). The model
uses a Frank copula, Laplace measurement error distribution, and one
mixture component for both the outcome (son's log income) and treatment
(father's log income) equations. Bootstrap standard errors are based on
200 replications.

## Usage

``` r
nlsy97_tsme_fit
```

## Format

A list of class `tsme` with fields documented in
[`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md). Key
fields:

- me_cop_param:

  Estimated Frank copula parameter (ME-corrected)

- me_tmat:

  ME-corrected 4x4 intergenerational transition matrix

- me_up_mob:

  ME-corrected upward mobility by father's income quartile

- me_spearman:

  ME-corrected Spearman rank correlation

- me_cond_quant:

  ME-corrected conditional quantile curves at 10th, 50th, and 90th
  percentiles across nine values of father's log income

- me_qyx:

  Fitted `merr` object for the child income equation

- me_qtx:

  Fitted `merr` object for the father income equation

## Source

Derived from the
[`nlsy97`](https://bcallaway11.github.io/qrme/reference/nlsy97.md)
dataset; see that entry for provenance. Model estimated by
[`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md).

## Details

The large intermediate `rqs` objects (`nome_qyx`, `nome_qtx`, `qrytx`)
have been dropped to keep the dataset compact; all fields needed by
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and `autoplot()` are
retained.

## References

Callaway, B., Li, T., Murtazashvili, I., and Tsyawo, E. S. (2026).
Distributional Effects with Two-Sided Measurement Error: An Application
to Intergenerational Income Mobility.
[doi:10.48550/arXiv.2107.09235](https://doi.org/10.48550/arXiv.2107.09235)
