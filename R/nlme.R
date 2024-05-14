#-----------------------------------------------------------------------------
# functions for nonlinear models with measurement error
#-----------------------------------------------------------------------------

#' @title compute.nlme
#' @description does the heavy lifting for computing nonlinear models with measurement error
#'
#' @inheritParams nlme
#'
#' @keywords internal
#' @export
compute.nlme <- function(data, Yformla, Tformla, tau, tvals, xdf=NULL,
                         copType="gaussian", simstep="MH", ndraws=250,
                         reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE,
                         reportPov=TRUE,
                         povline=log(20000), reportQ=c(.1,.5,.9),
                         Ynmix=1, Tnmix=1, tol=1, iters=400,
                         burnin=200, drawsd=4, ignore_me=FALSE,
                         messages=FALSE) {
  
  yname <- lhs.vars(Yformla)
  tname <- lhs.vars(Tformla)

  Qyx <- NULL
  Qtx <- NULL
  meres <- NULL
  
  if (!ignore_me) {
    Qyx <- qrme(Yformla, data=data, tau=tau, nmix=Ynmix, simstep=simstep, 
                tol=tol, iters=iters, burnin=burnin, drawsd=drawsd,
                messages=messages, se=FALSE) # don't bootstrap these
    Qtx  <- qrme(Tformla, data=data, tau=tau, nmix=Tnmix, simstep=simstep,
                 tol=tol, drawsd=drawsd, messages=messages, se=FALSE)

    # now get joint distribution
    meres <- qr2me(yname=yname,
                   tname=tname,
                   xformla=Yformla,
                   tau=tau,
                   data=data,
                   tvals=tvals,
                   xdf=xdf,
                   copula=copType,
                   Qyx=Qyx,
                   Qtx=Qtx,
                   ndraws=ndraws,
                   retFytxlist=FALSE,
                   messages=FALSE)
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
                   tvals=tvals, xdf=xdf, ndraws=ndraws,
                   copula=copType, 
                   Qyx=rqyx, Qtx=rqtx, retFytxlist=FALSE,
                   messages=FALSE)



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
              meQyx=Qyx, meQtx=Qtx, meresQ=meresQ, 
              nomeQyx=rqyx, nomeQtx=rqtx, nomeresQ=nomeresQ, 
              qrytxQ=qrytxQ, qrytx=qrytx,
              meTmat=meTmat, nomeTmat=nomeTmat, obsTmat=obsTmat, 
              mePs=mePs, nomePs=nomePs, obsPs=obsPs, 
              meUm=meUm, nomeUm=nomeUm, obsUm=obsUm,
              mePovrate=mePovrate, nomePovrate=nomePovrate, qrPovrate=qrPovrate)

  out

}



#' @title nlme
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
#' @param copType what type of copula to use in second step.  Options are
#'  "gaussian" (the default), "clayton", or "gumbel"
#' @param simstep whether to use an MH algorithm ("MH") or an importance
#'  sampling algorithm ("ImpSamp")
#' @param ndraws number of draws to use in MH algorithm to estimate first
#'  step quantile regressions (default 250)
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
#' @param tol tolerance for first step quantile regression model to converge
#'  (default is 1).  Note that convergence  will be sensitive to \code{length(tau)}
#'  and the number of mixture components included in the model
#' @param iters the number of MCMC iterations (default is 400)
#' @param burnin the number of MCMC iterations to drop (default is 200)
#' @param drawsd starting value of standard deviations of mixture components
#' @param ignore_me whether or not to ignore measurement error (this is primarily
#'  a way to get speedy calculations using copula-based approach)
#' @param messages whether or not to report details of computation as they
#'  occur (default is \code{FALSE})
#' @param se whether or not to estimate standard errors using the boostrap
#'  (default is FALSE)
#' @param biters if computing standard errors, the number of bootstrap iterations
#'  to use (default is 100)
#' @param cl allows for parallel processing in computing standard errors using
#'  the bootstrap (the default is 1)
#'
#' @return list of nonlinear measures of intergenerational income mobility
#'  adjusted for measurement error
#' @export
nlme <- function(data, Yformla, Tformla, tau, tvals, xdf=NULL, copType="gaussian",
                 simstep="MH", ndraws=250,
                 reportTmat=TRUE, reportSP=TRUE, reportUM=TRUE,
                 reportPov=TRUE, povline=log(20000), reportQ=c(.1,.5,.9),
                 Ynmix=1, Tnmix=1, tol=1, iters=400, burnin=200,
                 drawsd=4, ignore_me=FALSE, messages=FALSE,
                 se=FALSE, biters=100, cl=1) {

  res <- compute.nlme(data=data,
                      Yformla=Yformla,
                      Tformla=Tformla,
                      tau=tau,
                      tvals=tvals,
                      xdf=xdf,
                      copType=copType,
                      simstep=simstep,
                      ndraws=ndraws,
                      reportTmat=reportTmat,
                      reportSP=reportSP,
                      reportUM=reportUM,
                      reportPov=reportPov,
                      povline=povline,
                      reportQ=reportQ,
                      Ynmix=Ynmix,
                      Tnmix=Tnmix,
                      tol=tol,
                      iters=iters,
                      burnin=burnin,
                      drawsd=drawsd,
                      ignore_me=ignore_me,
                      messages=messages)
  

  if (se) {

    eachIter <- pbapply::pblapply(1:biters, function(b) {
      n <- nrow(data)
      brows <- sample(1:n, size=n, replace=TRUE)
      bdata <- data[brows,]
      tryCatch({
        out <- compute.nlme(data=bdata,
                            Yformla=Yformla,
                            Tformla=Tformla,
                            tau=tau,
                            tvals=tvals,
                            xdf=xdf,
                            copType=copType,
                            ndraws=ndraws,
                            reportTmat=reportTmat,
                            reportSP=reportSP,
                            reportUM=reportUM,
                            reportPov=reportPov,
                            povline=povline,
                            Ynmix=Ynmix,
                            Tnmix=Tnmix,
                            tol=tol,
                            drawsd=drawsd,
                            messages=messages)
        out
      }, error=function(cond) {
        return(NULL) # use this as code for error on that bootstrap iteration
        #return(cond)
      })
    }, cl=cl)

    
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

    res$biters <- length(eachIter)
    
  }

  res
}
