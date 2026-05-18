# =============================================================================
# Title: Package dataset documentation
# Description: Roxygen2 documentation for datasets bundled with qrme.
# Author: Brant Callaway
# Last update: 2026-05-18
# Date created: 2026-05-18
# =============================================================================

#' NLSY97 intergenerational income mobility sample
#'
#' A father-son analysis sample from the National Longitudinal Survey of Youth
#' 1997 (NLSY97), constructed for use in Callaway, Li, Murtazashvili, and
#' Tsyawo (2026).  The sample consists of father-son pairs with non-missing
#' father's education, observed in 2014-2015.  Log income variables are
#' pre-computed for convenience.
#'
#' @format A data frame with 1,066 rows and 16 variables:
#' \describe{
#'   \item{oldest_son}{1 if the son is the oldest son in the household, 0 otherwise}
#'   \item{ageC_97}{son's age in 1997}
#'   \item{ageC_14}{son's age in 2014}
#'   \item{whiteC}{1 if son is non-Hispanic white, 0 otherwise}
#'   \item{blackC}{1 if son is non-Hispanic Black, 0 otherwise}
#'   \item{hispanicC}{1 if son is Hispanic, 0 otherwise}
#'   \item{incF}{father's annual income (dollars)}
#'   \item{incS}{son's annual income (dollars)}
#'   \item{ageF}{father's age}
#'   \item{educF}{father's years of schooling}
#'   \item{ageM}{mother's age}
#'   \item{educM}{mother's years of schooling}
#'   \item{msa}{1 if the family lives in a metropolitan statistical area, 0 otherwise}
#'   \item{asvab}{son's ASVAB percentile score}
#'   \item{lci}{log son's income (\code{log(incS)})}
#'   \item{lpi}{log father's income (\code{log(incF)})}
#' }
#'
#' @source Bureau of Labor Statistics, National Longitudinal Survey of Youth
#'   1997 (\url{https://www.bls.gov/nls/nlsy97.htm}).  Father-son analysis
#'   sample constructed by Irina Murtazashvili.
#'
#' @references Callaway, B., Li, T., Murtazashvili, I., and Tsyawo, E. S.
#'   (2026).  Distributional Effects with Two-Sided Measurement Error: An
#'   Application to Intergenerational Income Mobility.
#'   \doi{10.48550/arXiv.2107.09235}
"nlsy97"

#' Pre-computed tsme() fit for the NLSY97 intergenerational mobility application
#'
#' A trimmed `tsme` object from fitting the two-sided measurement error model
#' to the \code{\link{nlsy97}} data, as reported in Callaway, Li, Murtazashvili,
#' and Tsyawo (2026).  The model uses a Frank copula, Laplace measurement error
#' distribution, and one mixture component for both the outcome (son's log
#' income) and treatment (father's log income) equations.  Bootstrap standard
#' errors are based on 200 replications.
#'
#' The large intermediate \code{rqs} objects (\code{nome_qyx}, \code{nome_qtx},
#' \code{qrytx}) have been dropped to keep the dataset compact; all fields
#' needed by \code{print()}, \code{summary()}, and \code{autoplot()} are
#' retained.
#'
#' @format A list of class \code{tsme} with fields documented in
#'   \code{\link{tsme}}.  Key fields:
#' \describe{
#'   \item{me_cop_param}{Estimated Frank copula parameter (ME-corrected)}
#'   \item{me_tmat}{ME-corrected 4x4 intergenerational transition matrix}
#'   \item{me_up_mob}{ME-corrected upward mobility by father's income quartile}
#'   \item{me_spearman}{ME-corrected Spearman rank correlation}
#'   \item{me_cond_quant}{ME-corrected conditional quantile curves at 10th, 50th,
#'     and 90th percentiles across nine values of father's log income}
#'   \item{me_qyx}{Fitted \code{merr} object for the child income equation}
#'   \item{me_qtx}{Fitted \code{merr} object for the father income equation}
#' }
#'
#' @source Derived from the \code{\link{nlsy97}} dataset; see that entry for
#'   provenance.  Model estimated by \code{\link{tsme}}.
#'
#' @references Callaway, B., Li, T., Murtazashvili, I., and Tsyawo, E. S.
#'   (2026).  Distributional Effects with Two-Sided Measurement Error: An
#'   Application to Intergenerational Income Mobility.
#'   \doi{10.48550/arXiv.2107.09235}
"nlsy97_tsme_fit"
