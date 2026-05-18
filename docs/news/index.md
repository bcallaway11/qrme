# Changelog

## qrme 1.0.1

- Added Emmanuel S. Tsyawo (<estsyawo@gmail.com>) as a listed author in
  DESCRIPTION.
- Bundled the NLSY97 father-son analysis sample (`nlsy97`) as package
  data.
- Switched from `exportPattern` to explicit `@export` tags; C++
  internals and helper functions are no longer part of the public API.
- Added `@export` to
  [`tmat()`](https://bcallaway11.github.io/qrme/reference/tmat.md),
  which was inadvertently unexported.
- Removed the `R/old/` legacy directory.
- Updated minimum R version to 4.0.0 and added `URL` / `BugReports`
  fields to DESCRIPTION.
- Added two vignettes: an introduction to
  [`qrme()`](https://bcallaway11.github.io/qrme/reference/qrme.md) and
  an application of
  [`tsme()`](https://bcallaway11.github.io/qrme/reference/tsme.md) to
  intergenerational income mobility using NLSY97 data.

## qrme 1.0.0

- Initial release.
