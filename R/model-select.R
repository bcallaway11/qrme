# =============================================================================
# Title: Model Selection for qrme and tsme
# Description: qrme_nmix_select() selects the number of ME mixture components
#   via AIC/BIC. tsme_model_select() fits tsme() over a grid of copula families
#   and ME distributions and returns an AIC/BIC comparison table.
# Author: Brant Callaway
# Last update: 2026-05-26
# Date created: 2026-05-17
# =============================================================================

#' @title qrme_nmix_select
#' @description Fit \code{\link{qrme}} with a range of mixture component counts
#'   and rank the specifications by AIC and BIC.  Use this to choose the number
#'   of ME mixture components for a single equation.
#'
#'   The number of free ME parameters is:
#'   \itemize{
#'     \item \code{n_mix = 0}: \eqn{k = 0} (no ME; LL from standard QR)
#'     \item \code{n_mix = 1}: \eqn{k = 2} (\eqn{\mu}, \eqn{\sigma})
#'     \item \code{n_mix >= 2}: \eqn{k = 3m - 1}
#'       (\eqn{\mu_1,\ldots,\mu_m}, \eqn{\sigma_1,\ldots,\sigma_m},
#'       \eqn{\pi_1,\ldots,\pi_{m-1}})
#'   }
#'
#' @param formula formula for the outcome model (passed to \code{\link{qrme}})
#' @param data data.frame
#' @param tau quantile levels for QR
#' @param n_mix_vals integer vector of component counts to evaluate
#'   (default \code{0:3}; \code{0} is the no-ME baseline)
#' @param return_fits logical; if \code{TRUE} the fitted qrme objects are
#'   included in the output (default \code{FALSE})
#' @param ... additional arguments passed to \code{\link{qrme}} for every fit
#'   (e.g. \code{me_distribution}, \code{mcmc_draws}, \code{mcmc_burn_in}).
#'   \code{n_mix} and \code{se} are set internally and must not be passed.
#'
#' @return a list with:
#'   \describe{
#'     \item{\code{table}}{data.frame sorted by BIC with columns \code{n_mix},
#'       \code{k_me}, \code{ll}, \code{AIC}, \code{BIC}}
#'     \item{\code{fits}}{(only when \code{return_fits = TRUE}) named list of
#'       qrme objects; the \code{n_mix = 0} entry is \code{NULL}}
#'   }
#'
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   n <- 300; X <- runif(n)
#'   Y <- 1 + 2 * X + rnorm(n) + rnorm(n, sd = 0.5)
#'   sel <- qrme_nmix_select(
#'     Y ~ X, data.frame(Y, X),
#'     tau          = c(0.25, 0.5, 0.75),
#'     n_mix_vals   = 0:3,
#'     mcmc_draws   = 100L,
#'     mcmc_burn_in = 50L
#'   )
#'   sel$table
#' }
#'
#' @export
qrme_nmix_select <- function(formula, data, tau,
                             n_mix_vals  = 0:3,
                             return_fits = FALSE,
                             ...) {
  n          <- nrow(data)
  dots_clean <- list(...)[!names(list(...)) %in% c("n_mix", "se")]

  # free ME parameter count per component count
  k_me <- function(m) {
    if (m == 0L) return(0L)
    if (m == 1L) return(2L)
    3L * m - 1L
  }

  # n_mix = 0: no-ME log-likelihood via QR asymmetric Laplace LL (sig -> 0)
  ll_nome <- function() {
    mf     <- model.frame(formula, data)
    y      <- model.response(mf)
    x      <- model.matrix(formula, data)
    betmat <- t(coef(quantreg::rq(formula, tau = tau, data = data)))
    loglik_raw(y = y, x = x, betmat = betmat, tau = tau,
               pi = 1, mu = 0, sig = 1e-10,
               me_distribution = "normal", ndraws = 1L)
  }

  results <- lapply(n_mix_vals, function(m) {
    if (m == 0L) {
      ll <- tryCatch(ll_nome(), error = function(e) NA_real_)
      return(list(n_mix = 0L, k = 0L, ll = ll, fit = NULL))
    }
    fit <- tryCatch(
      do.call(qrme, c(
        list(formula = formula, data = data, tau = tau, n_mix = m, se = FALSE),
        dots_clean
      )),
      error = function(e) {
        message(sprintf("  n_mix = %d failed: %s", m, conditionMessage(e)))
        NULL
      }
    )
    ll <- if (is.null(fit)) {
      NA_real_
    } else {
      tryCatch(as.numeric(logLik(fit)), error = function(e) NA_real_)
    }
    list(n_mix = m, k = k_me(m), ll = ll, fit = fit)
  })

  tbl <- data.frame(
    n_mix = sapply(results, `[[`, "n_mix"),
    k_me  = sapply(results, `[[`, "k"),
    ll    = sapply(results, `[[`, "ll"),
    stringsAsFactors = FALSE
  )
  tbl$AIC <- -2 * tbl$ll + 2 * tbl$k_me
  tbl$BIC <- -2 * tbl$ll + tbl$k_me * log(n)
  tbl <- tbl[order(tbl$BIC), ]
  rownames(tbl) <- NULL

  out <- list(table = tbl)
  if (return_fits) {
    fits        <- lapply(results, `[[`, "fit")
    names(fits) <- paste0("n_mix", sapply(results, `[[`, "n_mix"))
    out$fits    <- fits
  }
  out
}

#' @title tsme_model_select
#' @description Fit \code{\link{tsme}} over a grid of copula families and
#'   measurement error distributions and rank the specifications by AIC and BIC.
#'
#'   The log-likelihood is computed as the sum of three components:
#'   \itemize{
#'     \item \code{ll_y}: log-likelihood of the outcome ME model
#'     \item \code{ll_t}: log-likelihood of the treatment ME model
#'     \item \code{ll_cop}: copula log-likelihood
#'   }
#'   \deqn{AIC = -2 \ell + 2k, \quad BIC = -2 \ell + k \log n}
#'   where \eqn{k} is computed automatically from \code{y_n_mix} and
#'   \code{t_n_mix}: \eqn{k = k_Y + k_T + 1}, with \eqn{k_m = 2} for
#'   \eqn{m = 1} and \eqn{k_m = 3m - 1} for \eqn{m \ge 2}, plus one
#'   copula parameter.
#'
#' @param data data.frame passed to \code{\link{tsme}}
#' @param y_formula formula for the outcome model
#' @param t_formula formula for the treatment model
#' @param tau quantile levels for the first-stage QR (passed to
#'   \code{\link{tsme}})
#' @param t_values treatment values for conditional distribution summaries
#' @param copulas character vector of copula families to try.  Any subset of
#'   \code{c("gaussian","clayton","gumbel","frank")} (default: all four).
#' @param me_distributions character vector of ME distributions to try.  Any
#'   subset of \code{c("gaussian","laplace")} (default: both).
#' @param y_n_mix number of ME mixture components for the outcome equation
#'   (default 1)
#' @param t_n_mix number of ME mixture components for the treatment equation
#'   (default 1)
#' @param n_cores number of parallel workers (default 1); set to
#'   \code{length(copulas) * length(me_distributions)} to run one worker per
#'   grid cell
#' @param return_fits logical; if \code{TRUE} the fitted tsme objects are
#'   returned alongside the table (default \code{FALSE})
#' @param ... additional arguments passed to \code{\link{tsme}} for every cell
#'   in the grid (e.g. \code{n_copula_me_draws}, \code{mcmc_draws})
#'
#' @return a list with element \code{table} (a data.frame sorted by BIC) and,
#'   if \code{return_fits = TRUE}, element \code{fits} (a list of tsme objects
#'   in the same order as \code{table})
#'
#' @examples
#' \dontrun{
#'   tau      <- seq(0.05, 0.95, length.out = 15)
#'   t_values <- quantile(nlsy97$lpi, probs = seq(0.1, 0.9, by = 0.1))
#'   sel <- tsme_model_select(
#'     data             = nlsy97,
#'     y_formula        = lci ~ ageC_97 + ageF,
#'     t_formula        = lpi ~ ageC_97 + ageF,
#'     tau              = tau,
#'     t_values         = t_values,
#'     y_n_mix          = 1L,
#'     t_n_mix          = 1L,
#'     n_cores          = 4L,
#'     mcmc_draws       = 200L,
#'     mcmc_burn_in     = 20L
#'   )
#'   print(sel)
#' }
#'
#' @export
tsme_model_select <- function(data, y_formula, t_formula, tau, t_values,
                              copulas          = c("gaussian", "clayton",
                                                   "gumbel",   "frank"),
                              me_distributions = c("gaussian", "laplace"),
                              y_n_mix          = 1L,
                              t_n_mix          = 1L,
                              n_cores          = 1L,
                              return_fits      = FALSE,
                              ...) {
  grid <- expand.grid(
    copula          = copulas,
    me_distribution = me_distributions,
    stringsAsFactors = FALSE
  )
  n <- nrow(data)

  # free ME params per equation: k=2 for n_mix=1, k=3m-1 for n_mix>=2
  k_me <- function(m) if (m == 1L) 2L else 3L * m - 1L
  k_params <- k_me(y_n_mix) + k_me(t_n_mix) + 1L  # +1 for copula parameter

  fit_one <- function(i) {
    cop  <- grid$copula[i]
    dist <- grid$me_distribution[i]
    message(sprintf("  Fitting: copula=%-10s me_distribution=%s", cop, dist))
    fit <- tryCatch(
      tsme(data = data, y_formula = y_formula, t_formula = t_formula,
           tau = tau, t_values = t_values,
           copula = cop, me_distribution = dist,
           y_n_mix = y_n_mix, t_n_mix = t_n_mix,
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
  }

  results <- if (n_cores > 1L) {
    parallel::mclapply(seq_len(nrow(grid)), fit_one, mc.cores = n_cores)
  } else {
    lapply(seq_len(nrow(grid)), fit_one)
  }

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

  out <- list(table = tbl, k_params = k_params)
  if (return_fits) out$fits <- lapply(res, `[[`, "fit")
  out
}
