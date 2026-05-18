# NLSY97 intergenerational income mobility sample

A father-son analysis sample from the National Longitudinal Survey of
Youth 1997 (NLSY97), constructed for use in Callaway, Li, Murtazashvili,
and Tsyawo (2026). The sample consists of father-son pairs with
non-missing father's education, observed in 2014-2015. Log income
variables are pre-computed for convenience.

## Usage

``` r
nlsy97
```

## Format

A data frame with 1,066 rows and 16 variables:

- oldest_son:

  1 if the son is the oldest son in the household, 0 otherwise

- ageC_97:

  son's age in 1997

- ageC_14:

  son's age in 2014

- whiteC:

  1 if son is non-Hispanic white, 0 otherwise

- blackC:

  1 if son is non-Hispanic Black, 0 otherwise

- hispanicC:

  1 if son is Hispanic, 0 otherwise

- incF:

  father's annual income (dollars)

- incS:

  son's annual income (dollars)

- ageF:

  father's age

- educF:

  father's years of schooling

- ageM:

  mother's age

- educM:

  mother's years of schooling

- msa:

  1 if the family lives in a metropolitan statistical area, 0 otherwise

- asvab:

  son's ASVAB percentile score

- lci:

  log son's income (`log(incS)`)

- lpi:

  log father's income (`log(incF)`)

## Source

Bureau of Labor Statistics, National Longitudinal Survey of Youth 1997
(<https://www.bls.gov/nls/nlsy97.htm>). Father-son analysis sample
constructed by Irina Murtazashvili.

## References

Callaway, B., Li, T., Murtazashvili, I., and Tsyawo, E. S. (2026).
Distributional Effects with Two-Sided Measurement Error: An Application
to Intergenerational Income Mobility.
[doi:10.48550/arXiv.2107.09235](https://doi.org/10.48550/arXiv.2107.09235)
