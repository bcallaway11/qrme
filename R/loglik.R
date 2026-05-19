# Internal helper: compute observed-data log-likelihood from raw parameters via
# Monte Carlo integration. Called by logLik.merr (ndraws=500) and by em.algo
# at each outer iteration for convergence checking (ndraws=100).
loglik_raw <- function(y, x, betmat, tau, pi, mu, sig, me_distribution, ndraws = 100L) {
  m <- length(sig)
  n <- length(y)

  if (me_distribution == "laplace") {
    vdraws <- rlaplace(ndraws, mu = 0, sigma = sig[1])
  } else {
    comp   <- sample(seq_len(m), ndraws, replace = TRUE, prob = pi)
    vdraws <- rnorm(ndraws, mean = mu[comp], sd = sig[comp])
  }

  log_fy <- vapply(seq_len(n), function(i) {
    xi       <- as.matrix(x[i, ])
    fyx_vals <- vapply(y[i] - vdraws, function(ys) fyxC(ys, betmat, xi, tau), numeric(1))
    log(max(mean(fyx_vals), .Machine$double.eps))
  }, numeric(1))

  sum(log_fy)
}


#' @title logLik.merr
#'
#' @description Observed-data log-likelihood for a \code{merr} object via
#'   Monte Carlo integration over the measurement error distribution. Enables
#'   \code{AIC(fit)} and \code{BIC(fit)} automatically.
#'
#' @param object A \code{merr} object returned by \code{qrme()}.
#' @param ndraws Number of Monte Carlo draws for the integral (default 500).
#' @param ... Currently unused.
#'
#' @return A scalar of class \code{"logLik"} with attributes \code{df}
#'   (number of free parameters) and \code{nobs} (sample size).
#'
#' @examples
#' \donttest{
#'   set.seed(42)
#'   n <- 200; X <- runif(n)
#'   Y <- 1 + 2 * X + rnorm(n) + rnorm(n, sd = 0.5)
#'   fit1 <- qrme(Y ~ X, data.frame(Y, X), tau = 0.5,
#'                n_mix = 1L, mcmc_draws = 80L, mcmc_burn_in = 40L,
#'                max_em_iters = 10L, verbose = FALSE)
#'   fit2 <- qrme(Y ~ X, data.frame(Y, X), tau = 0.5,
#'                n_mix = 2L, mcmc_draws = 80L, mcmc_burn_in = 40L,
#'                max_em_iters = 10L, verbose = FALSE)
#'   data.frame(n_mix = 1:2,
#'              AIC   = c(AIC(fit1), AIC(fit2)),
#'              BIC   = c(BIC(fit1), BIC(fit2)))
#' }
#'
#' @export
logLik.merr <- function(object, ndraws = 500, ...) {
  y       <- object$y
  x       <- object$x
  betmat  <- object$bet
  tau     <- object$tau
  pi      <- object$pi
  mu      <- object$mu
  sig     <- object$sig
  me_distribution <- object$me_distribution
  m       <- length(sig)
  n       <- length(y)

  ll_val <- loglik_raw(y, x, betmat, tau, pi, mu, sig, me_distribution, ndraws = ndraws)

  ## df: L*K (beta) + (m-1) (pi, sums-to-1) + (m-1) (mu, mean-zero) + m (sigma)
  L  <- length(tau)
  K  <- ncol(x)
  df <- L * K + 3 * m - 2

  structure(ll_val, df = df, nobs = n, class = "logLik")
}
