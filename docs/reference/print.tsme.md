# print.tsme

Print method for tsme objects. Shows formulas, copula parameters,
measurement error standard deviations, and Spearman's rho for all three
estimators. Standard errors are shown when bootstrap results are
available.

## Usage

``` r
# S3 method for class 'tsme'
print(x, ...)
```

## Arguments

- x:

  a tsme object returned by
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md)

- ...:

  unused

## Value

`x`, invisibly.

## Examples

``` r
print(nlsy97_tsme_fit)
#> Two-Sided Measurement Error Model (tsme)
#> ------------------------------------------
#> Outcome:    lci ~ ageC_97 + ageF 
#> Treatment:  lpi ~ ageC_97 + ageF 
#> 
#> Copula: frank      theta(ME) = 2.1152 (SE 0.6443) | theta(No-ME) = 1.3026 (SE 0.2153)
#> 
#> ME parameters:
#>   Y (laplace): sigma = 0.4980 (SE 0.0057)
#>   T (laplace): sigma = 0.4837 (SE 0.0075)
#> 
#> Spearman's rho:  ME = 0.3664 (SE 0.0798) | No-ME = 0.2295 (SE 0.0330) | Obs = 0.2085 (SE 0.0289)
```
