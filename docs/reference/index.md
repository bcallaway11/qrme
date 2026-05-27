# Package index

## Package overview

- [`qrme-package`](https://bcallaway11.github.io/qrme/reference/qrme-package.md)
  : qrme: Quantile regression with measurement error

## Main estimators

The two primary user-facing estimators.
[`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md) handles
measurement error in the outcome;
[`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) handles
measurement error in both the outcome and a continuous treatment.

- [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md) :
  qrme
- [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) :
  tsme

## Model selection and diagnostics

Tools for choosing among copula families, ME distributions, and mixture
orders, and for checking convergence sensitivity.

- [`tsme_model_select()`](https://bcallaway11.github.io/qrme/reference/tsme_model_select.md)
  : tsme_model_select
- [`qrme_nmix_select()`](https://bcallaway11.github.io/qrme/reference/qrme_nmix_select.md)
  : qrme_nmix_select
- [`qrme_start_search()`](https://bcallaway11.github.io/qrme/reference/qrme_start_search.md)
  : qrme_start_search
- [`logLik(`*`<merr>`*`)`](https://bcallaway11.github.io/qrme/reference/logLik.merr.md)
  : logLik.merr

## Mobility measures

Functions for computing intergenerational (or rank-rank) mobility
parameters from a joint distribution.

- [`tmat()`](https://bcallaway11.github.io/qrme/reference/tmat.md) :
  tmat
- [`upMob()`](https://bcallaway11.github.io/qrme/reference/upMob.md) :
  upMob

## S3 methods

Print, summarise, and plot fitted objects.

- [`print(`*`<merr>`*`)`](https://bcallaway11.github.io/qrme/reference/print.merr.md)
  : print.merr
- [`print(`*`<tsme>`*`)`](https://bcallaway11.github.io/qrme/reference/print.tsme.md)
  : print.tsme
- [`summary(`*`<tsme>`*`)`](https://bcallaway11.github.io/qrme/reference/summary.tsme.md)
  : summary.tsme
- [`autoplot(`*`<tsme>`*`)`](https://bcallaway11.github.io/qrme/reference/autoplot.tsme.md)
  : autoplot.tsme
- [`plot(`*`<tsme>`*`)`](https://bcallaway11.github.io/qrme/reference/plot.tsme.md)
  : plot.tsme

## Data

NLSY97 father-son analysis sample and a pre-computed
[`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) result
for use in the vignettes and examples.

- [`nlsy97`](https://bcallaway11.github.io/qrme/reference/nlsy97.md) :
  NLSY97 intergenerational income mobility sample
- [`nlsy97_tsme_fit`](https://bcallaway11.github.io/qrme/reference/nlsy97_tsme_fit.md)
  : Pre-computed tsme() fit for the NLSY97 intergenerational mobility
  application

## Lower-level functions

Exported for advanced users and reproducibility. Most users will not
call these directly.

- [`compute.qrme()`](https://bcallaway11.github.io/qrme/reference/compute.qrme.md)
  : compute.qrme
- [`compute.tsme()`](https://bcallaway11.github.io/qrme/reference/compute.tsme.md)
  : compute.tsme
- [`em.algo()`](https://bcallaway11.github.io/qrme/reference/em.algo.md)
  : em.algo
- [`em.algo.inner()`](https://bcallaway11.github.io/qrme/reference/em.algo.inner.md)
  : Inner part of EM-algorithm for QR with measurement error
- [`betfun()`](https://bcallaway11.github.io/qrme/reference/betfun.md) :
  betfun
- [`betfun.inner()`](https://bcallaway11.github.io/qrme/reference/betfun.inner.md)
  : betfun.inner
- [`makeRQS()`](https://bcallaway11.github.io/qrme/reference/makeRQS.md)
  : makeRQS
- [`rlaplace()`](https://bcallaway11.github.io/qrme/reference/rlaplace.md)
  : rlaplace
- [`mh_mcmc()`](https://bcallaway11.github.io/qrme/reference/mh_mcmc.md)
  : mh_mcmc
- [`fv()`](https://bcallaway11.github.io/qrme/reference/fv.md) : fv
- [`fv.yx()`](https://bcallaway11.github.io/qrme/reference/fv.yx.md) :
  fv.yx
- [`fy.x()`](https://bcallaway11.github.io/qrme/reference/fy.x.md) :
  fy.x
- [`cop.pdf()`](https://bcallaway11.github.io/qrme/reference/cop.pdf.md)
  : cop.pdf
- [`ll()`](https://bcallaway11.github.io/qrme/reference/ll.md) : ll
- [`parms2coppar()`](https://bcallaway11.github.io/qrme/reference/parms2coppar.md)
  : parms2coppar
- [`qr2meobj()`](https://bcallaway11.github.io/qrme/reference/qr2meobj.md)
  : qr2meobj
- [`qr2me()`](https://bcallaway11.github.io/qrme/reference/qr2me.md) :
  QR with 2-sided measurement error
- [`addplot()`](https://bcallaway11.github.io/qrme/reference/addplot.md)
  : addplot
