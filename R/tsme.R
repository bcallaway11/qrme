#-----------------------------------------------------------------------------
# functions for nonlinear models with measurement error
#-----------------------------------------------------------------------------

#' @title compute.tsme
#' @description does the heavy lifting for computing nonlinear models with measurement error
#'
#' @inheritParams tsme
#'
#' @keywords internal
#' @export
compute.tsme <- function(data, Yformla, Tformla, tau, tvals, xdf=NULL,
                         me_dist=me_dist, copType="gaussian",
                         mcmc_method="MH", copula_me_draws=100L,
                         reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE,
                         reportPov=TRUE,
                         povline=log(20000), reportQ=c(.1,.5,.9),
                         Ynmix=1, Tnmix=1, tol=NULL, conv_crit="params",
                         ndraws_ll=100L, conv_patience=1L, maxit=100L,
                         mcmc_draws=200, mcmc_burnin=100, proposal_sd=NULL,
                         mobility_copula_draws=100L, ignore_me=FALSE,
                         verbose=FALSE) {
  
  yname <- lhs.vars(Yformla)
  tname <- lhs.vars(Tformla)

  Qyx <- NULL
  Qtx <- NULL
  meres <- NULL
  
  if (!ignore_me) {
    qrme_progress(verbose, "tsme first stage: outcome measurement-error model")
    Qyx <- qrme(Yformla, data=data, tau=tau, me_dist=me_dist, nmix=Ynmix, mcmc_method=mcmc_method,
                tol=tol, conv_crit=conv_crit, ndraws_ll=ndraws_ll, conv_patience=conv_patience,
                mcmc_draws=mcmc_draws, mcmc_burnin=mcmc_burnin, proposal_sd=proposal_sd,
                maxit=maxit, verbose=verbose, se=FALSE) # don't bootstrap these
    qrme_progress(verbose, "tsme first stage: treatment measurement-error model")
    Qtx  <- qrme(Tformla, data=data, tau=tau, me_dist=me_dist, nmix=Tnmix, mcmc_method=mcmc_method,
                 tol=tol, conv_crit=conv_crit, ndraws_ll=ndraws_ll, conv_patience=conv_patience,
                 mcmc_draws=mcmc_draws, mcmc_burnin=mcmc_burnin, proposal_sd=proposal_sd,
                 maxit=maxit, verbose=verbose, se=FALSE)

    # now get joint distribution
    qrme_progress(verbose, "tsme second stage: measurement-error copula model")
    meres <- qr2me(yname=yname,
                   tname=tname,
                   xformla=Yformla,
                   tau=tau,
                   data=data,
                   tvals=tvals,
                   xdf=xdf,
                   me_dist=me_dist,
                   copula=copType,
                   Qyx=Qyx,
                   Qtx=Qtx,
                   copula_me_draws=copula_me_draws,
                   mobility_copula_draws=mobility_copula_draws,
                   retFytxlist=FALSE,
                   verbose=verbose)
  }
  
 
 
  # get results without measurement error
  rqyx <- rq(Yformla, tau=tau, data=data)
  rqtx <- rq(Tformla, tau=tau, data=data)
  class(rqyx) <- c("merr", "rqs")
  class(rqtx) <- c("merr", "rqs")
  rqyx$pi <- rqtx$pi <- 1
  rqyx$mu <- rqtx$mu <- 0
  rqyx$sig <- rqtx$sig <- 0

  nomeres <- qr2me(yname, tname, Yformla, tau=tau, data=data,
                   tvals=tvals, xdf=xdf, copula_me_draws=copula_me_draws,
                   mobility_copula_draws=mobility_copula_draws,
                   copula=copType, 
                   Qyx=rqyx, Qtx=rqtx, retFytxlist=FALSE,
                   verbose=qrme_verbose_level(verbose) >= 2L)



  if (reportTmat) {
    meTmat <- meres$t_mat#tmat(Qyx$Ystar, Qtx$Ystar)
    obsTmat <- tmat(data[,yname], data[,tname])
    nomeTmat <- nomeres$t_mat
  }

  if (reportSP) {
    mePs <- meres$Ps
    nomePs <- nomeres$Ps
    obsPs  <- cor(data[,yname], data[,tname], method="spearman")
  }

  if (reportUM) {
    meUm <- meres$up_mob
    nomeUm <- nomeres$up_mob
    obsUm <- upMob(data[,yname], data[,tname])
  }


  # results just using quantile regression
  tau <- seq(0,1,length.out=100)
  qrformla <- toformula(yname, c(tname, rhs.vars(Yformla)))
  qrytx <- rq(qrformla, tau=tau, data=data)
  qrytx$Fyt <- lapply(1:length(tvals), function(i) {
    if (is.null(xdf)) newdta <- data else newdta <- xdf
    newdta[,tname] <- tvals[i]
    if (nrow(newdta) == 1) {
      Qytx <- predict(qrytx, newdata=newdta)
      Fytx <- BMisc::makeDist(Qytx[1,], tau, rearrange=TRUE, method="linear")
      return(Fytx)
    } else {
      Fytx <- predict(qrytx, newdata=newdta, type="Fhat", stepfun=TRUE)
      Fytx <- rearrange(Fytx)
      combineDfs(seq(min(data[,yname]), max(data[,yname]), length.out=500), Fytx)
    }
  })
  
  if (reportPov) {
    if (ignore_me) mePovrate <- NULL else mePovrate <- sapply(1:(length(meres$tvals)), function(i) meres$Fyt[[i]](povline))
    nomePovrate <- sapply(1:(length(nomeres$tvals)), function(i) nomeres$Fyt[[i]](povline))
    qrPovrate <- sapply(1:(length(nomeres$tvals)), function(i) qrytx$Fyt[[i]](povline))
  }

  meresQ <- if(ignore_me) meresQ <- NULL else meresQ <-  sapply(meres$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))
  nomeresQ <- sapply(nomeres$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))
  qrytxQ <- sapply(qrytx$Fyt, function(Fy) quantile(Fy, type=1, probs=reportQ))

  meCopParam <- meres$cop.param
  nomeCopParam <- nomeres$cop.param
  
  out <- list(Yformla=Yformla, Tformla=Tformla, tau=tau, tvals=tvals, copType=copType,
              meCopParam=meCopParam, nomeCopParam=nomeCopParam, Ynmix=Ynmix, Tnmix=Tnmix, reportQ=reportQ,
              tol=c(Y=Qyx$tol, T=Qtx$tol),
              conv_crit=c(Y=Qyx$conv_crit, T=Qtx$conv_crit),
              conv_criteria=c(Y=Qyx$conv_criteria, T=Qtx$conv_criteria),
              conv_converged=c(Y=Qyx$conv_converged, T=Qtx$conv_converged),
              mix_n_iter=c(Y=Qyx$mix_n_iter, T=Qtx$mix_n_iter),
              mix_loglik=c(Y=Qyx$mix_loglik, T=Qtx$mix_loglik),
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
#' @param Yformla formula for outcome model
#' @param Tformla formula for treatment model
#' @param tau values of tau to estimate first step quantile regressions for
#' @param tvals values of the treatment to compute conditional distribution-type
#'  parameters for
#' @param xdf matrix of values of covariates to average over for conditional
#'  distribution-type parameters.   The default is NULL and in this case
#'  all covariates in the data will be averaged over.  A main alternative
#'  would be to pass in a single row with particular values of covariates
#'  of interest.
#' #' @param me_dist the distribution of the measurement error.  "gaussian" is the
#'  default and supports a mixture of normals.  "laplace" is also supported
#' @param copType what type of copula to use in second step.  Options are
#'  "gaussian" (the default), "clayton", or "gumbel"
#' @param mcmc_method simulation step to use in the EM algorithm. Currently
#'  only \code{"MH"} is supported.
#' @param copula_me_draws Number of measurement-error draws used in the
#'  second-stage copula likelihood (default 100). Increase this for a less
#'  noisy copula likelihood at the cost of speed.
#' @param reportTmat whether or not to report a transition matrix
#' @param reportSP whether or not to report Spearman's rho (rank-rank correlation)
#' @param reportUM whether or not to report upward mobility parameters
#' @param reportPov whether or not to report fraction of population below
#'  the poverty line as a function of parents' income
#' @param povline value of the poverty line (default log(20000))
#' @param reportQ quantiles of child's income as a function of parents' income
#'  to report (default is .1,.5,.9)
#' @param Ynmix number of mixture components for outcome measurement error
#'  model
#' @param Tnmix number of mixture components for treatment measurement error
#'  model
#' @param tol Convergence tolerance.  When \code{NULL} (default), a value is
#'  chosen automatically based on \code{conv_crit}.  See \code{\link{em.algo}}
#'  for details.
#' @param conv_crit Convergence criterion passed to \code{\link{qrme}}:
#'  \code{"params"} (default) or \code{"loglik"}.  See \code{\link{em.algo}}.
#' @param ndraws_ll Number of Monte Carlo draws for the log-likelihood
#'  convergence check when \code{conv_crit = "loglik"} (default 100).
#'  Ignored when \code{conv_crit = "params"}.
#' @param conv_patience Integer. Consecutive iterations that must satisfy the
#'  convergence criterion before stopping (default 1).
#' @param maxit Maximum number of EM outer iterations per \code{qrme} call
#'  (default 100).
#' @param mcmc_draws total number of MCMC draws per EM step (default 200)
#' @param mcmc_burnin number of MCMC draws to discard as burnin (default 100)
#' @param proposal_sd Standard deviation of the MH proposal used in the two
#'  first-stage \code{\link{qrme}} calls. When \code{NULL} (default), each
#'  first stage uses its own automatic scale, \code{sqrt(var(Y))} or
#'  \code{sqrt(var(T))}.
#' @param mobility_copula_draws Number of copula draws per covariate row used
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
#' @param ncores allows for parallel processing in computing standard errors using
#'  the bootstrap (the default is 1)
#'
#' @return list of nonlinear measures of intergenerational income mobility
#'  adjusted for measurement error
#' @export
tsme <- function(data, Yformla, Tformla, tau, tvals, xdf=NULL,
                 me_dist="gaussian", copType="gaussian",
                 mcmc_method="MH", copula_me_draws=100L,
                 reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE,
                 reportPov=TRUE, povline=log(20000), reportQ=c(.1,.5,.9),
                 Ynmix=1, Tnmix=1, tol=NULL, conv_crit="params",
                 ndraws_ll=100L, conv_patience=1L, maxit=100L,
                 mcmc_draws=200, mcmc_burnin=100, proposal_sd=NULL,
                 mobility_copula_draws=100L, ignore_me=FALSE, verbose=FALSE,
                 se=FALSE, n_boot=100, ncores=1) {

  res <- compute.tsme(data=data,
                      Yformla=Yformla,
                      Tformla=Tformla,
                      tau=tau,
                      tvals=tvals,
                      xdf=xdf,
                      me_dist=me_dist,
                      copType=copType,
                      mcmc_method=mcmc_method,
                      copula_me_draws=copula_me_draws,
                      reportTmat=reportTmat,
                      reportSP=reportSP,
                      reportUM=reportUM,
                      reportPov=reportPov,
                      povline=povline,
                      reportQ=reportQ,
                      Ynmix=Ynmix,
                      Tnmix=Tnmix,
                      tol=tol,
                      conv_crit=conv_crit,
                      ndraws_ll=ndraws_ll,
                      conv_patience=conv_patience,
                      maxit=maxit,
                      mcmc_draws=mcmc_draws,
                      mcmc_burnin=mcmc_burnin,
                      proposal_sd=proposal_sd,
                      mobility_copula_draws=mobility_copula_draws,
                      ignore_me=ignore_me,
                      verbose=verbose)
  

  if (se) {

    eachIter <- pbapply::pblapply(1:n_boot, function(b) {
      n <- nrow(data)
      brows <- sample(1:n, size=n, replace=TRUE)
      bdata <- data[brows,]
      tryCatch({
        out <- compute.tsme(data=bdata,
                            Yformla=Yformla,
                            Tformla=Tformla,
                            tau=tau,
                            tvals=tvals,
                            xdf=xdf,
                            me_dist=me_dist,
                            copType=copType,
                            mcmc_method=mcmc_method,
                            copula_me_draws=copula_me_draws,
                            reportTmat=reportTmat,
                            reportSP=reportSP,
                            reportUM=reportUM,
                            reportPov=reportPov,
                            povline=povline,
                            Ynmix=Ynmix,
                            Tnmix=Tnmix,
                            tol=tol,
                            conv_crit=conv_crit,
                            ndraws_ll=ndraws_ll,
                            conv_patience=conv_patience,
                              maxit=maxit,
                              mcmc_draws=mcmc_draws,
                              mcmc_burnin=mcmc_burnin,
                              proposal_sd=proposal_sd,
                              mobility_copula_draws=mobility_copula_draws,
                              verbose=FALSE)
        out
      }, error=function(cond) {
        return(NULL) # use this as code for error on that bootstrap iteration
        #return(cond)
      })
    }, cl=ncores)

    
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
