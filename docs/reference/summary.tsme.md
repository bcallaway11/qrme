# summary.tsme

Summary method for tsme objects. Calls
[`print.tsme`](https://bcallaway11.github.io/qrme/reference/print.tsme.md)
then appends the upward mobility table and ME-corrected transition
matrix.

## Usage

``` r
# S3 method for class 'tsme'
summary(object, ...)
```

## Arguments

- object:

  a tsme object returned by
  [`tsme`](https://bcallaway11.github.io/qrme/reference/tsme.md)

- ...:

  unused

## Value

`object`, invisibly.

## Examples

``` r
summary(nlsy97_tsme_fit)
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
#> 
#> Upward Mobility:
#>                 ME  ME_se   NoME NoME_se Observed Obs_se
#> Q1 (bottom) 0.7836 0.0211 0.8264  0.0084   0.8303 0.0181
#> Q2          0.5756 0.0112 0.5897  0.0058   0.5725 0.0253
#> Q3          0.4218 0.0119 0.4054  0.0060   0.3868 0.0288
#> Q4 (top)    0.2071 0.0207 0.1771  0.0081   0.1667 0.0248
#> 
#> Transition Matrix (ME-corrected, rows=father, cols=son):
#>                 Q1     Q2     Q3     Q4
#> Q1 (bottom) 0.4231 0.2809 0.1848 0.1113
#> Q2          0.2802 0.2889 0.2479 0.1830
#> Q3          0.1836 0.2476 0.2865 0.2823
#> Q4 (top)    0.1131 0.1826 0.2808 0.4235
```
