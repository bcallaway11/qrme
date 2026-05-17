# =============================================================================
# Title: Start Value Search for qrme
# Description: qrme_start_search() runs qrme() from multiple random starting
#   values and returns the best fit by log-likelihood, helping diagnose
#   convergence sensitivity.
# Author: Brant Callaway
# Last update: 2026-05-17
# Date created: 2026-05-17
# =============================================================================

#' @title qrme_start_search
#' @description Run \code{\link{qrme}} from multiple random starting values
#'   and return the fit with the highest log-likelihood.  Use this to assess
#'   whether the EM algorithm converges to a unique solution or is sensitive
#'   to starting values.
#'
#'   For each restart, \code{start_sigma} is drawn uniformly from
#'   \code{sigma_range}, \code{start_mu} is drawn from \eqn{N(0, 0.01)},
#'   and mixture weights \code{start_pi} are drawn from a symmetric Dirichlet
#'   (i.e., uniform on the simplex).
#'
#' @param formula formula for the outcome model (passed to \code{\link{qrme}})
#' @param data data.frame
#' @param tau quantile levels for QR
#' @param n_starts number of random restarts (default 10)
#' @param seed optional integer seed for reproducibility (default \code{NULL})
#' @param sigma_range length-2 numeric vector giving the \eqn{[\min, \max]}
#'   interval from which \code{start_sigma} is drawn for each component
#'   (default \code{c(0.1, 2.0)})
#' @param return_fits logical; if \code{TRUE} the fitted qrme objects for all
#'   successful restarts are returned (default \code{FALSE})
#' @param ... additional arguments passed to \code{\link{qrme}} for every
#'   restart (e.g. \code{n_mix}, \code{me_distribution}, \code{mcmc_draws}).
#'   Arguments \code{start_sigma}, \code{start_mu}, \code{start_pi}, and
#'   \code{se} are set internally and should not be passed via \code{...}.
#'
#' @return a list with:
#'   \describe{
#'     \item{\code{table}}{data.frame of start index and log-likelihood for
#'       each successful restart, sorted descending by log-likelihood}
#'     \item{\code{best_fit}}{the qrme fit with the highest log-likelihood}
#'     \item{\code{fits}}{(only when \code{return_fits = TRUE}) list of all
#'       successful qrme fits}
#'   }
#'
#' @export
qrme_start_search <- function(formula, data, tau,
                              n_starts    = 10L,
                              seed        = NULL,
                              sigma_range = c(0.1, 2.0),
                              return_fits = FALSE,
                              ...) {
  if (!is.null(seed)) set.seed(seed)
  dots  <- list(...)
  n_mix <- if (!is.null(dots$n_mix)) dots$n_mix else 1L

  # strip start values and se from ... — set internally per restart
  reserved <- c("start_sigma", "start_mu", "start_pi", "se")
  dots_clean <- dots[!names(dots) %in% reserved]

  run_one <- function(i) {
    start_sigma <- sort(stats::runif(n_mix, sigma_range[1L], sigma_range[2L]))
    start_mu    <- stats::rnorm(n_mix, mean = 0, sd = 0.1)
    start_pi    <- if (n_mix > 1L) {
      p <- stats::runif(n_mix)
      p / sum(p)
    } else {
      1
    }
    tryCatch(
      do.call(qrme, c(
        list(formula     = formula,
             data        = data,
             tau         = tau,
             start_sigma = start_sigma,
             start_mu    = start_mu,
             start_pi    = start_pi,
             se          = FALSE),
        dots_clean
      )),
      error = function(e) {
        message(sprintf("  Start %d failed: %s", i, conditionMessage(e)))
        NULL
      }
    )
  }

  fits <- lapply(seq_len(n_starts), run_one)
  ok   <- !sapply(fits, is.null)
  if (!any(ok)) stop("All starts failed; check the model specification.")

  ok_idx <- seq_len(n_starts)[ok]
  lls    <- sapply(ok_idx, function(i) as.numeric(logLik(fits[[i]])))

  tbl <- data.frame(
    start = ok_idx,
    ll    = lls,
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(-tbl$ll), ]
  rownames(tbl) <- NULL

  best_fit <- fits[ok][[which.max(lls)]]
  out      <- list(table = tbl, best_fit = best_fit)
  if (return_fits) out$fits <- fits[ok]
  out
}
