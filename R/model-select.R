# =============================================================================
# Title: Model Selection for tsme
# Description: tsme_model_select() fits tsme() over a grid of copula families
#   and ME distributions and returns an AIC/BIC comparison table.
# Author: Brant Callaway
# Last update: 2026-05-17
# Date created: 2026-05-17
# =============================================================================

#' @title tsme_model_select
#' @description Fit \code{\link{tsme}} over a grid of copula families and
#'   measurement error distributions and rank the specifications by AIC and BIC.
#'
#'   The log-likelihood is computed as the sum of three components:
#'   \itemize{
#'     \item \code{ll_y}: log-likelihood of the outcome ME model (\code{logLik(fit$me_qyx)})
#'     \item \code{ll_t}: log-likelihood of the treatment ME model (\code{logLik(fit$me_qtx)})
#'     \item \code{ll_cop}: copula log-likelihood (\code{fit$me_cop_loglik})
#'   }
#'   \deqn{AIC = -2 \ell + 2k, \quad BIC = -2 \ell + k \log n}
#'   where \eqn{k} is the number of free parameters specified by \code{k_params}.
#'   For \code{y_n_mix = t_n_mix = 1} the free parameters are \eqn{\sigma_Y},
#'   \eqn{\sigma_T}, and \eqn{\theta_{\text{cop}}}, giving \code{k_params = 3}.
#'
#' @param data data.frame passed to \code{\link{tsme}}
#' @param y_formula formula for the outcome model
#' @param t_formula formula for the treatment model
#' @param tau quantile levels for the first-stage QR (passed to \code{\link{tsme}})
#' @param t_values treatment values for conditional distribution summaries
#' @param copulas character vector of copula families to try.  Any subset of
#'   \code{c("gaussian","clayton","gumbel","frank")} (default: all four).
#' @param me_distributions character vector of ME distributions to try.  Any
#'   subset of \code{c("gaussian","laplace")} (default: both).
#' @param k_params number of free parameters used in AIC/BIC (default 3,
#'   appropriate for \code{y_n_mix = t_n_mix = 1})
#' @param return_fits logical; if \code{TRUE} the fitted tsme objects are
#'   returned alongside the table (default \code{FALSE})
#' @param ... additional arguments passed to \code{\link{tsme}} for every cell
#'   in the grid (e.g. \code{n_copula_me_draws}, \code{mcmc_draws},
#'   \code{y_n_mix}, \code{t_n_mix})
#'
#' @return a list with element \code{table} (a data.frame sorted by BIC) and,
#'   if \code{return_fits = TRUE}, element \code{fits} (a list of tsme objects
#'   in the same order as \code{table})
#'
#' @export
tsme_model_select <- function(data, y_formula, t_formula, tau, t_values,
                              copulas          = c("gaussian", "clayton",
                                                   "gumbel",   "frank"),
                              me_distributions = c("gaussian", "laplace"),
                              k_params         = 3L,
                              return_fits      = FALSE,
                              ...) {
  grid <- expand.grid(
    copula          = copulas,
    me_distribution = me_distributions,
    stringsAsFactors = FALSE
  )
  n <- nrow(data)

  results <- lapply(seq_len(nrow(grid)), function(i) {
    cop  <- grid$copula[i]
    dist <- grid$me_distribution[i]
    message(sprintf("  Fitting: copula=%-10s me_distribution=%s", cop, dist))
    fit <- tryCatch(
      tsme(data = data, y_formula = y_formula, t_formula = t_formula,
           tau = tau, t_values = t_values,
           copula = cop, me_distribution = dist,
           se = FALSE, ...),
      error = function(e) {
        message(sprintf("    FAILED: %s", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    ll_y   <- as.numeric(logLik(fit$me_qyx))
    ll_t   <- as.numeric(logLik(fit$me_qtx))
    ll_cop <- fit$me_cop_loglik
    ll     <- ll_y + ll_t + ll_cop
    list(
      copula          = cop,
      me_distribution = dist,
      ll_y            = ll_y,
      ll_t            = ll_t,
      ll_cop          = ll_cop,
      ll              = ll,
      aic             = -2 * ll + 2 * k_params,
      bic             = -2 * ll + k_params * log(n),
      fit             = fit
    )
  })

  ok  <- !sapply(results, is.null)
  res <- results[ok]

  tbl <- data.frame(
    copula          = sapply(res, `[[`, "copula"),
    me_distribution = sapply(res, `[[`, "me_distribution"),
    ll_y            = sapply(res, `[[`, "ll_y"),
    ll_t            = sapply(res, `[[`, "ll_t"),
    ll_cop          = sapply(res, `[[`, "ll_cop"),
    ll              = sapply(res, `[[`, "ll"),
    aic             = sapply(res, `[[`, "aic"),
    bic             = sapply(res, `[[`, "bic"),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  tbl <- tbl[order(tbl$bic), ]
  rownames(tbl) <- NULL

  out <- list(table = tbl)
  if (return_fits) out$fits <- lapply(res, `[[`, "fit")
  out
}
