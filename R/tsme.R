# =============================================================================
# Title: Two-Sided Measurement Error Models
# Description: Functions for nonlinear models with two-sided measurement error
#   via copula (tsme) and its internal workhorse (compute.tsme).
# Author: Brant Callaway
# Last update: 2026-05-15
# Date created: 2026-05-14
# =============================================================================

#' @title compute.tsme
#' @description does the heavy lifting for computing nonlinear models with measurement error
#'
#' @inheritParams tsme
#'
#' @keywords internal
#' @export
compute.tsme <- function(data, y_formula, t_formula, tau, t_values, x_data=NULL,
                         me_distribution="gaussian", copula="gaussian",
                         mcmc_method="MH", n_copula_me_draws=100L,
                         report_transition_matrix=TRUE, report_spearman=TRUE, report_upward_mobility=TRUE,
                         report_poverty=TRUE,
                         pov_line=log(20000), report_quantiles=c(.1,.5,.9),
                         y_n_mix=1, t_n_mix=1, tol=NULL, conv_criterion="params",
                         loglik_draws=100L, conv_patience=1L, max_em_iters=100L,
                         mcmc_draws=200, mcmc_burn_in=100, proposal_sd=NULL,
                         start_sigma_y=NULL, start_mu_y=NULL, start_pi_y=NULL,
                         start_sigma_t=NULL, start_mu_t=NULL, start_pi_t=NULL,
                         n_copula_draws=100L, ignore_me=FALSE,
                         verbose=FALSE) {
  
  y_name <- lhs.vars(y_formula)
  t_name <- lhs.vars(t_formula)

  Qyx <- NULL
  Qtx <- NULL
  meres <- NULL
  
  if (!ignore_me) {
    qrme_progress(verbose, "tsme first stage: outcome measurement-error model")
    Qyx <- qrme(y_formula, data=data, tau=tau, me_distribution=me_distribution, n_mix=y_n_mix, mcmc_method=mcmc_method,
                tol=tol, conv_criterion=conv_criterion, loglik_draws=loglik_draws, conv_patience=conv_patience,
                mcmc_draws=mcmc_draws, mcmc_burn_in=mcmc_burn_in, proposal_sd=proposal_sd,
                start_sigma=start_sigma_y, start_mu=start_mu_y, start_pi=start_pi_y,
                max_em_iters=max_em_iters, verbose=verbose, se=FALSE) # don't bootstrap these
    qrme_progress(verbose, "tsme first stage: treatment measurement-error model")
    Qtx  <- qrme(t_formula, data=data, tau=tau, me_distribution=me_distribution, n_mix=t_n_mix, mcmc_method=mcmc_method,
                 tol=tol, conv_criterion=conv_criterion, loglik_draws=loglik_draws, conv_patience=conv_patience,
                 mcmc_draws=mcmc_draws, mcmc_burn_in=mcmc_burn_in, proposal_sd=proposal_sd,
                 start_sigma=start_sigma_t, start_mu=start_mu_t, start_pi=start_pi_t,
                 max_em_iters=max_em_iters, verbose=verbose, se=FALSE)

    # now get joint distribution
    qrme_progress(verbose, "tsme second stage: measurement-error copula model")
    meres <- qr2me(y_name=y_name,
                   t_name=t_name,
                   x_formula=y_formula,
                   tau=tau,
                   data=data,
                   t_values=t_values,
                   x_data=x_data,
                   me_distribution=me_distribution,
                   copula=copula,
                   Qyx=Qyx,
                   Qtx=Qtx,
                   n_copula_me_draws=n_copula_me_draws,
                   n_copula_draws=n_copula_draws,
                   return_fytx_list=FALSE,
                   verbose=verbose)
  }
  
 
 
  # get results without measurement error
  rqyx <- rq(y_formula, tau=tau, data=data)
  rqtx <- rq(t_formula, tau=tau, data=data)
  class(rqyx) <- c("merr", "rqs")
  class(rqtx) <- c("merr", "rqs")
  rqyx$pi <- rqtx$pi <- 1
  rqyx$mu <- rqtx$mu <- 0
  rqyx$sig <- rqtx$sig <- 0

  nomeres <- qr2me(y_name, t_name, y_formula, tau=tau, data=data,
                   t_values=t_values, x_data=x_data, n_copula_me_draws=n_copula_me_draws,
                   n_copula_draws=n_copula_draws,
                   copula=copula, 
                   Qyx=rqyx, Qtx=rqtx, return_fytx_list=FALSE,
                   verbose=qrme_verbose_level(verbose) >= 2L)



  if (report_transition_matrix) {
    meTmat <- meres$t_mat#tmat(Qyx$Ystar, Qtx$Ystar)
    obsTmat <- tmat(data[,y_name], data[,t_name])
    nomeTmat <- nomeres$t_mat
  }

  if (report_spearman) {
    mePs <- meres$Ps
    nomePs <- nomeres$Ps
    obsPs  <- cor(data[,y_name], data[,t_name], method="spearman")
  }

  if (report_upward_mobility) {
    meUm <- meres$up_mob
    nomeUm <- nomeres$up_mob
    obsUm <- upMob(data[,y_name], data[,t_name])
  }


  # results just using quantile regression
  tau <- seq(0,1,length.out=100)
  qr_formula <- toformula(y_name, c(t_name, rhs.vars(y_formula)))
  qrytx <- rq(qr_formula, tau=tau, data=data)
  qrytx$Fyt <- lapply(1:length(t_values), function(i) {
    if (is.null(x_data)) newdta <- data else newdta <- x_data
    newdta[,t_name] <- t_values[i]
    if (nrow(newdta) == 1) {
      Qytx <- predict(qrytx, newdata=newdta)
      Fytx <- BMisc::makeDist(Qytx[1,], tau, rearrange=TRUE, method="linear")
      return(Fytx)
    } else {
      Fytx <- predict(qrytx, newdata=newdta, type="Fhat", stepfun=TRUE)
      Fytx <- rearrange(Fytx)
      combineDfs(seq(min(data[,y_name]), max(data[,y_name]), length.out=500), Fytx)
    }
  })
  
  if (report_poverty) {
    if (ignore_me) mePovrate <- NULL else mePovrate <- sapply(1:(length(meres$t_values)), function(i) meres$Fyt[[i]](pov_line))
    nomePovrate <- sapply(1:(length(nomeres$t_values)), function(i) nomeres$Fyt[[i]](pov_line))
    qrPovrate <- sapply(1:(length(nomeres$t_values)), function(i) qrytx$Fyt[[i]](pov_line))
  }

  meresQ <- if(ignore_me) meresQ <- NULL else meresQ <-  sapply(meres$Fyt, function(Fy) quantile(Fy, type=1, probs=report_quantiles))
  nomeresQ <- sapply(nomeres$Fyt, function(Fy) quantile(Fy, type=1, probs=report_quantiles))
  qrytxQ <- sapply(qrytx$Fyt, function(Fy) quantile(Fy, type=1, probs=report_quantiles))

  meCopParam   <- meres$cop.param
  nomeCopParam <- nomeres$cop.param
  meCopLoglik   <- meres$cop_loglik
  nomeCopLoglik <- nomeres$cop_loglik

  out <- list(y_formula=y_formula, t_formula=t_formula, tau=tau, t_values=t_values, copula=copula,
              meCopParam=meCopParam, nomeCopParam=nomeCopParam,
              meCopLoglik=meCopLoglik, nomeCopLoglik=nomeCopLoglik,
              mix_loglik=c(Y=Qyx$mix_loglik, T=Qtx$mix_loglik),
              y_n_mix=y_n_mix, t_n_mix=t_n_mix, report_quantiles=report_quantiles,
              tol=c(Y=Qyx$tol, T=Qtx$tol),
              conv_criterion=c(Y=Qyx$conv_criterion, T=Qtx$conv_criterion),
              conv_criteria=c(Y=Qyx$conv_criteria, T=Qtx$conv_criteria),
              conv_converged=c(Y=Qyx$conv_converged, T=Qtx$conv_converged),
              mix_n_iter=c(Y=Qyx$mix_n_iter, T=Qtx$mix_n_iter),
              mix_converged=c(Y=Qyx$mix_converged, T=Qtx$mix_converged),
              meQyx=Qyx, meQtx=Qtx, meresQ=meresQ,
              nomeQyx=rqyx, nomeQtx=rqtx, nomeresQ=nomeresQ,
              qrytxQ=qrytxQ, qrytx=qrytx,
              meTmat=meTmat, nomeTmat=nomeTmat, obsTmat=obsTmat,
              mePs=mePs, nomePs=nomePs, obsPs=obsPs,
              meUm=meUm, nomeUm=nomeUm, obsUm=obsUm,
              mePovrate=mePovrate, nomePovrate=nomePovrate, qrPovrate=qrPovrate)

  out

}



#' @title tsme
#' @description function to compute nonlinear models with two sided measurement error
#' 
#' @param data data.frame
#' @param y_formula formula for outcome model
#' @param t_formula formula for treatment model
#' @param tau values of tau to estimate first step quantile regressions for
#' @param t_values values of the treatment to compute conditional distribution-type
#'  parameters for
#' @param x_data matrix of values of covariates to average over for conditional
#'  distribution-type parameters.   The default is NULL and in this case
#'  all covariates in the data will be averaged over.  A main alternative
#'  would be to pass in a single row with particular values of covariates
#'  of interest.
#' #' @param me_distribution the distribution of the measurement error.  "gaussian" is the
#'  default and supports a mixture of normals.  "laplace" is also supported
#' @param copula what type of copula to use in second step.  Options are
#'  "gaussian" (the default), "clayton", "gumbel", or "frank"
#' @param mcmc_method simulation step to use in the EM algorithm. Currently
#'  only \code{"MH"} is supported.
#' @param n_copula_me_draws Number of measurement-error draws used in the
#'  second-stage copula likelihood (default 100). Increase this for a less
#'  noisy copula likelihood at the cost of speed.
#' @param report_transition_matrix whether or not to report a transition matrix
#' @param report_spearman whether or not to report Spearman's rho (rank-rank correlation)
#' @param report_upward_mobility whether or not to report upward mobility parameters
#' @param report_poverty whether or not to report fraction of population below
#'  the poverty line as a function of parents' income
#' @param pov_line value of the poverty line (default log(20000))
#' @param report_quantiles quantiles of child's income as a function of parents' income
#'  to report (default is .1,.5,.9)
#' @param y_n_mix Number of Gaussian mixture components for the outcome (Y)
#'   measurement error distribution. Default is 1 (single Gaussian). Use
#'   \code{logLik()}, \code{AIC()}, and \code{BIC()} on a fitted
#'   \code{qrme()} object to select the appropriate value. \code{y_n_mix} and
#'   \code{t_n_mix} can differ — it is valid and common to use different
#'   complexity for each equation.
#' @param t_n_mix Number of Gaussian mixture components for the treatment (T)
#'   measurement error distribution. Default is 1. See \code{y_n_mix}.
#' @param tol Convergence tolerance.  When \code{NULL} (default), a value is
#'  chosen automatically based on \code{conv_criterion}.  See \code{\link{em.algo}}
#'  for details.
#' @param conv_criterion Convergence criterion passed to \code{\link{qrme}}:
#'  \code{"params"} (default) or \code{"loglik"}.  See \code{\link{em.algo}}.
#' @param loglik_draws Number of Monte Carlo draws for the log-likelihood
#'  convergence check when \code{conv_criterion = "loglik"} (default 100).
#'  Ignored when \code{conv_criterion = "params"}.
#' @param conv_patience Integer. Consecutive iterations that must satisfy the
#'  convergence criterion before stopping (default 1).
#' @param max_em_iters Maximum number of EM outer iterations per \code{qrme} call
#'  (default 100).
#' @param mcmc_draws total number of MCMC draws per EM step (default 200)
#' @param mcmc_burn_in number of MCMC draws to discard as burnin (default 100)
#' @param proposal_sd Standard deviation of the MH proposal used in the two
#'  first-stage \code{\link{qrme}} calls. When \code{NULL} (default), each
#'  first stage uses its own automatic scale, \code{sqrt(var(Y))} or
#'  \code{sqrt(var(T))}.
#' @param start_sigma_y Starting value(s) for the ME standard deviation(s) in
#'  the outcome (Y) \code{\link{qrme}} call. A vector of length \code{y_n_mix}.
#'  When \code{NULL} (default), \code{qrme} chooses its own starting value.
#' @param start_mu_y Starting value(s) for the ME mean(s) in the outcome (Y)
#'  \code{\link{qrme}} call. A vector of length \code{y_n_mix}. \code{NULL}
#'  uses the \code{qrme} default.
#' @param start_pi_y Starting mixture weight(s) for the outcome (Y)
#'  \code{\link{qrme}} call. A vector of length \code{y_n_mix} that sums to 1.
#'  \code{NULL} uses the \code{qrme} default.
#' @param start_sigma_t Starting value(s) for the ME standard deviation(s) in
#'  the treatment (T) \code{\link{qrme}} call. A vector of length \code{t_n_mix}.
#'  \code{NULL} uses the \code{qrme} default.
#' @param start_mu_t Starting value(s) for the ME mean(s) in the treatment (T)
#'  \code{\link{qrme}} call. \code{NULL} uses the \code{qrme} default.
#' @param start_pi_t Starting mixture weight(s) for the treatment (T)
#'  \code{\link{qrme}} call. \code{NULL} uses the \code{qrme} default.
#' @param n_copula_draws Number of copula draws per covariate row used
#'  to compute mobility summaries such as transition matrices and rank
#'  correlations (default 100).
#' @param ignore_me whether or not to ignore measurement error (this is primarily
#'  a way to get speedy calculations using copula-based approach)
#' @param verbose Logical or nonnegative integer. If \code{FALSE} (default),
#'  suppresses progress output. If \code{TRUE} or \code{1}, reports major
#'  computational stages, EM iteration numbers, and convergence diagnostics.
#'  If \code{2}, also reports qrme-native details such as \code{pi},
#'  \code{mu}, \code{sigma}, tolerance, finite-mixture summaries, and the
#'  no-measurement-error comparison copula step. If \code{3} or larger, also
#'  reports raw \code{normalmixEM()} output, which uses \code{lambda} for
#'  mixture probabilities.
#' @param se whether or not to estimate standard errors using the boostrap
#'  (default is FALSE)
#' @param n_boot if computing standard errors, the number of bootstrap iterations
#'  to use (default is 100)
#' @param n_cores allows for parallel processing in computing standard errors using
#'  the bootstrap (the default is 1)
#'
#' @return list of nonlinear measures of intergenerational income mobility
#'  adjusted for measurement error
#' @export
tsme <- function(data, y_formula, t_formula, tau, t_values, x_data=NULL,
                 me_distribution="gaussian", copula="gaussian",
                 mcmc_method="MH", n_copula_me_draws=100L,
                 report_transition_matrix=TRUE, report_spearman=TRUE, report_upward_mobility=TRUE,
                 report_poverty=TRUE, pov_line=log(20000), report_quantiles=c(.1,.5,.9),
                 y_n_mix=1, t_n_mix=1, tol=NULL, conv_criterion="params",
                 loglik_draws=100L, conv_patience=1L, max_em_iters=100L,
                 mcmc_draws=200, mcmc_burn_in=100, proposal_sd=NULL,
                 start_sigma_y=NULL, start_mu_y=NULL, start_pi_y=NULL,
                 start_sigma_t=NULL, start_mu_t=NULL, start_pi_t=NULL,
                 n_copula_draws=100L, ignore_me=FALSE, verbose=FALSE,
                 se=FALSE, n_boot=100, n_cores=1) {

  if (qrme_verbose_level(verbose) >= 1L) {
    message("tsme method")
    message("----------------------")
    message(
      "Citation: Callaway, Brantly, Tong Li, Irina Murtazashvili, ",
      "and Emmanuel Tsyawo, Distributional Effects with Two-Sided ",
      "Measurement Error: An Application to Intergenerational Income ",
      "Mobility, Working Paper, 2025."
    )
    message("----------------------")
  }

  res <- compute.tsme(data=data,
                      y_formula=y_formula,
                      t_formula=t_formula,
                      tau=tau,
                      t_values=t_values,
                      x_data=x_data,
                      me_distribution=me_distribution,
                      copula=copula,
                      mcmc_method=mcmc_method,
                      n_copula_me_draws=n_copula_me_draws,
                      report_transition_matrix=report_transition_matrix,
                      report_spearman=report_spearman,
                      report_upward_mobility=report_upward_mobility,
                      report_poverty=report_poverty,
                      pov_line=pov_line,
                      report_quantiles=report_quantiles,
                      y_n_mix=y_n_mix,
                      t_n_mix=t_n_mix,
                      tol=tol,
                      conv_criterion=conv_criterion,
                      loglik_draws=loglik_draws,
                      conv_patience=conv_patience,
                      max_em_iters=max_em_iters,
                      mcmc_draws=mcmc_draws,
                      mcmc_burn_in=mcmc_burn_in,
                      proposal_sd=proposal_sd,
                      start_sigma_y=start_sigma_y,
                      start_mu_y=start_mu_y,
                      start_pi_y=start_pi_y,
                      start_sigma_t=start_sigma_t,
                      start_mu_t=start_mu_t,
                      start_pi_t=start_pi_t,
                      n_copula_draws=n_copula_draws,
                      ignore_me=ignore_me,
                      verbose=verbose)
  

  if (se) {

    eachIter <- pbapply::pblapply(1:n_boot, function(b) {
      n <- nrow(data)
      brows <- sample(1:n, size=n, replace=TRUE)
      bdata <- data[brows,]
      tryCatch({
        out <- compute.tsme(data=bdata,
                            y_formula=y_formula,
                            t_formula=t_formula,
                            tau=tau,
                            t_values=t_values,
                            x_data=x_data,
                            me_distribution=me_distribution,
                            copula=copula,
                            mcmc_method=mcmc_method,
                            n_copula_me_draws=n_copula_me_draws,
                            report_transition_matrix=report_transition_matrix,
                            report_spearman=report_spearman,
                            report_upward_mobility=report_upward_mobility,
                            report_poverty=report_poverty,
                            pov_line=pov_line,
                            y_n_mix=y_n_mix,
                            t_n_mix=t_n_mix,
                            tol=tol,
                            conv_criterion=conv_criterion,
                            loglik_draws=loglik_draws,
                            conv_patience=conv_patience,
                            max_em_iters=max_em_iters,
                            mcmc_draws=mcmc_draws,
                            mcmc_burn_in=mcmc_burn_in,
                            proposal_sd=proposal_sd,
                            start_sigma_y=start_sigma_y,
                            start_mu_y=start_mu_y,
                            start_pi_y=start_pi_y,
                            start_sigma_t=start_sigma_t,
                            start_mu_t=start_mu_t,
                            start_pi_t=start_pi_t,
                            n_copula_draws=n_copula_draws,
                            verbose=FALSE)
        out
      }, error=function(cond) {
        return(NULL) # use this as code for error on that bootstrap iteration
        #return(cond)
      })
    }, cl=n_cores)

    
    # drop list elements where bootstrap failed
    eachIter <- eachIter[!sapply(eachIter, is.null)]

    # first step estimators
    # outcome
    res$meQyx$sig.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$sig)), 2, sd)
    res$meQyx$mu.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$mu)), 2, sd)
    res$meQyx$pi.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQyx$pi)), 2, sd)
    res$meQyx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) e$meQyx$bet)), 1:2, sd)
    # treatment
    res$meQtx$sig.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$sig)), 2, sd)
    res$meQtx$mu.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$mu)), 2, sd)
    res$meQtx$pi.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meQtx$pi)), 2, sd)
    res$meQtx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) e$meQtx$bet)), 1:2, sd)
    # qr
    res$qrytx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) coef(e$qrytx))), 1:2, sd)
    # nome
    res$nomeQyt$bet.se <- apply(simplify2array(lapply(eachIter, function(e) t(coef(e$nomeQyx)))), 1:2, sd)
    res$nomeQtx$bet.se <- apply(simplify2array(lapply(eachIter, function(e) t(coef(e$nomeQtx)))), 1:2, sd)
    
    # transition matrix
    res$meTmat.se <- apply(simplify2array(lapply(eachIter, function(e) e$meTmat)), 1:2, sd)
    res$nomeTmat.se <- apply(simplify2array(lapply(eachIter, function(e) e$nomeTmat)), 1:2, sd)
    res$obsTmat.se <- apply(simplify2array(lapply(eachIter, function(e) e$obsTmat)), 1:2, sd)

    # spearmans rho
    res$mePs.se <- sd(sapply(eachIter, function(e) e$mePs))
    res$nomePs.se <- sd(sapply(eachIter, function(e) e$nomePs))
    res$obsPs.se <- sd(sapply(eachIter, function(e) e$obsPs))

    # upward mobility
    res$meUm.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$meUm)), 2, sd)
    res$nomeUm.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$nomeUm)), 2, sd)
    res$obsUm.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$obsUm)), 2, sd)

    # poverty rate
    res$mePovrate.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$mePovrate)), 2, sd)
    res$nomePovrate.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$nomePovrate)), 2, sd)
    res$qrPovrate.se  <- apply(do.call(rbind, lapply(eachIter, function(e) e$qrPovrate)), 2, sd)

    # quantiles
    res$meresQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$meresQ)), 1:2, sd)
    res$nomeresQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$nomeresQ)), 1:2, sd)
    res$qrytxQ.se <- apply(simplify2array(lapply(eachIter, function(e) e$qrytxQ)), 1:2, sd)

    res$n_boot <- length(eachIter)
    
  }

  res
}
