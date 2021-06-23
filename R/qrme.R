#' @title compute.qrme
#' @description does the heavy lifting on computing quantile regression with
#'  left hand side measurement error
#'
#' @param compute.qrme
#'
#' @keywords internal
#' @export
compute.qrme <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL,
                         startsig=NULL, startpi=NULL, simstep="MH", tol=1, iters=400,
                         burnin=200, drawsd=4, cl=1, messages=FALSE) {
  xformla <- formla
  xformla[[2]] <- NULL ## drop y variable
  x <- model.matrix(xformla, data)
  yname <- as.character(formla[[2]])
  y <- data[,yname]
  k <- ncol(x) ## number of x variables
  m <- nmix
  if (is.null(startbet)) {
    ## this defaults the start values of the beta_0 to be
    ## the observed quantiles of the outcome and the other
    ## betas to be equal to 0
    betvals <- t(coef(quantreg::rq(formla, tau=tau, data=data)))
    
  } else {
    betvals <- startbet
  }
  if (is.null(startmu)) {
    muvals1 <- seq(-m, m, length.out=m)
    muvals <- muvals1 - mean(muvals1)
  } else {
    muvals <- startmu
  }
  if (is.null(startsig)) {
    sigvals <- rep(1,m)
  } else {
    sigvals <- startsig
  }
  if (is.null(startpi)) {
    pivals <- rep(1/m, m)
  } else {
    pivals <- startpi
  }
  qrparams <- list(bet=betvals, m=m, pi=pivals, mu=muvals, sig=sigvals)


  ## Estimate QR model with measurement error using EM algorithm
  res <- em.algo(formla, data,
                 betmatguess=betvals, tau=tau,
                 m=m, piguess=pivals, muguess=muvals,
                 sigguess=sigvals, simstep=simstep, tol=tol,
                 iters=iters, burnin=burnin, drawsd=drawsd, cl=cl, messages=messages)

  
  out <- makeRQS(res, formla, data, tau=tau)

  class(out) <- c("merr", class(out))

  ## set final parameters for measurement error class
  out$bet <- res$bet
  out$pi <- res$pi
  out$mu <- res$mu
  out$sig <- res$sig
  #out$Ystar <- res$Ystar

  out
}


#' @title qrme
#'
#' @description Quantile Regression with measurement error in the dependent
#'  variable using an EM-algorithm.  In practice, this function assumes
#'  that the measurement error distribution is a mixture of normal
#'  distributions.
#'
#' @param formla y ~ x
#' @param tau vector for which quantiles to compute quantile regression
#' @param data a data.frame that contains y and x
#' @param nmix The number of mixture components of the measurement error
#' @param startbet an LxK matrix of starting values for beta where
#'  L is the dimension of tau and K is the number of covariates (default is
#'  NULL and in this case, the starting values are set to be the QR
#'  coefficients coming from QR that ignores measurment error)
#' @param startmu A vector of length nmix of starting values for the mean
#'  of the mixture of normals distribution for the measurment error (default
#'  is NULL and in this case, the starting values are basically set to be
#'  equally spaced from -nmix to nmix but forced to have mean 0)
#' @param startsig A vector of length nmix of starting values for the
#'  standard deviation of the mixture of normals distribution for the
#'  measurement error (default is NULL and in this case, the starting values
#'  are all set to be 1)
#' @param startpi A vector of length nmix of starting values for the fraction
#'  of observations in each component of the mimxture of normals distribution
#'  for the measurement error (default is NULL and in this case, the starting
#'  values are all set to be 1/nmix)
#' @param simstep The type of simulation step to use in the EM algorithm.
#'  The default is "MH" for Metropolis-Hasting.  The alternative is
#'  "ImpSamp" for importance sampling. 
#' @param tol This is the convergence criteria.  When the change in the
#'  Euclidean distance between the new parameters (at each iteration) and
#'  the old parameters (from the previous iteration) is smaller than tol,
#'  the algorithm concludes.  In general, larger values for tol will result
#'  in a fewer number of iterations and smaller values will result in more
#'  accurate estimates.
#' @param iters How many iterations to use in the simulation step (default is
#'  400)
#' @param burnin How many iterations to drop in the simulation step (default
#'  is 200)
#' @param drawsd The starting standard deviation for the measurement error
#'  term.
#' @param cl The numbe of clusters to use for parallel computation (default
#'  is 1 so that computation is not done in parallel)
#' @param se Whether or not to compute standard errors using the bootstrap
#'  (default is FALSE)
#' @param biters Number of bootstrap iterations to use.  Only is considered
#'  in the case where computing standard errors (default is 100)
#' @param messages Whether or not to report details of estimation procedure
#'  (default is FALSE)
#' 
#' @return an object of class "merr"
#'
#' @export
qrme <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL,
                 startsig=NULL, startpi=NULL, simstep="MH", tol=1, iters=400,
                 burnin=200, drawsd=4, cl=1, se=FALSE, biters=100, messages=FALSE) {



  res <- compute.qrme(formla=formla,
                      tau=tau,
                      data=data,
                      nmix=nmix,
                      startbet=startbet,
                      startmu=startmu,
                      startsig=startsig,
                      startpi=startpi,
                      simstep=simstep,
                      tol=tol,
                      iters=iters,
                      burnin=burnin,
                      drawsd=drawsd,
                      cl=cl, messages=messages)

  if (se) {
    eachIter <- pbapply::pblapply(1:biters, function(b) {
      n <- nrow(data)
      brows <- sample(1:n, size=n, replace=TRUE)
      bdata <- data[brows,]
      tryCatch({
        out <- compute.qrme(formla=formla,
                            tau=tau,
                            data=bdata,
                            nmix=nmix,
                            startbet=res$bet,
                            startmu=res$mu,
                            startsig=res$sig,
                            startpi=res$pi,
                            simstep=simstep,
                            tol=tol,
                            iters=iters,
                            burnin=burnin,
                            drawsd=drawsd,
                            cl=1)
        out$Ystar=NULL ## just drop this because it takes up a lot of memory
        out
      }, error=function(cond) {
        return(NULL)
      })
    }, cl=cl)

    ## drop list elements where bootstrap failed
    eachIter <- eachIter[!sapply(eachIter, is.null)]

    ## only works if these are scalar...now only matrix
    res$sig.se <- apply(do.call(rbind, lapply(eachIter, function(e) e$sig)), 2, sd)
    res$mu.se <- apply(t(sapply(eachIter, function(e) e$mu)), 2, sd)##sd(sapply(eachIter, function(e) e$mu))
    res$pi.se <- apply(t(sapply(eachIter, function(e) e$pi)), 2, sd)
    res$bet.se <- apply(simplify2array(lapply(eachIter, function(e) e$bet)), 1:2, sd) ## element-wise standard deviation of list of matrices
    
  }

  res
}




#' @title QR with 2-sided measurement error
#'
#' @description Estimates QR parameters in the case with measurement error
#'  in the outcome and measurement error in a particular continuous "treatment"
#'  variable.
#'
#' @param yname name of the outcome in the passed in data
#' @param tname name of the treatment in the passed in data
#' @param xformla a one-sided formula for additional covariates (assumed
#'  not to be measured with error)
#' @param tau a vector containing particular quantiles that have been estimated ??
#' @param data a data.frame containing the data used for estimation
#' @param xdf If you want conditional distributions to be returned, pass in the value of the distribution here;
#'  otherwise the default behavior is to return a single distribution that averages over all values of X in the dataset
#' @param copula which type of copula to use (default is "gaussian")
#' @param Qyx quantile regression estimates (can be adjusted for measurement
#'  error) of Y on X
#' @param Qtx quantile regression estimates (can be adjusted for measurement
#'  error) of T on X
#' @param retFytxlist whether or not to return the conditional distribution
#'  for every value of x in xdf
#'  (default is FALSE because this can take up a lot of room in memory)
#' @inheritParams nlme
#' @inheritParams qrme
qr2me <- function(yname, tname, xformla, tau, data, xdf=NULL, tvals=NULL,
                  copula="gaussian",
                  Qyx, Qtx, retFytxlist=FALSE,
                  ndraws=100, messages=TRUE) {

  if (messages) {
    cat("\nqr2me method...\n")
    cat("----------------------")
    cat("\nCitation: Callaway, Brantly, Tong Li, and Irina Murtazashvili, Quantile Treatment Effects with Two-Sided Measurement Error, Working Paper, 2021....\n")
    cat("----------------------")
    cat("\n")
  }
  
  x <- model.matrix(xformla, data)
  n <-  nrow(data)

  tau_grid <- seq(0,1,length.out=100)
  Qyx_interpolated <- lapply(1:nrow(x), function(i) {
    xb <- as.numeric(t(as.matrix(x[i,]))%*%as.matrix(coef(Qyx)))
    sapply(tau_grid, function(tt_grid_val) {
      interpolateC(tau, xb, tt_grid_val, TRUE)
    })
  })
  Fyx <- lapply(Qyx_interpolated, function(Qyx) BMisc::makeDist(Qyx, tau_grid, sorted=TRUE))


  #Fyx1 <- predict(Qyx, newdata=as.data.frame(x)) ## this gives nxL matrix
  #Fyx <- lapply(1:nrow(Fyx1), function(i) BMisc::makeDist(Fyx1[i,], Fx=tau, sorted=TRUE))

  Qtx_interpolated <- lapply(1:nrow(x), function(i) {
    xb <- as.numeric(t(as.matrix(x[i,]))%*%as.matrix(coef(Qtx)))
    sapply(tau_grid, function(tt_grid_val) {
      interpolateC(tau, xb, tt_grid_val, TRUE)
    })
  })
  Ftx <- lapply(Qtx_interpolated, function(Qtx) BMisc::makeDist(Qtx, tau_grid, sorted=TRUE))
  
  #Ftx1 <- predict(Qtx, newdata=as.data.frame(x))
  #Ftx <- lapply(1:nrow(Ftx1), function(i) BMisc::makeDist(Ftx1[i,], Fx=tau, sorted=TRUE))

  ## Take passed in quantiles and create their conditional distribution
  ## Fyx <- predict(Qyx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)
  ## Ftx <- predict(Qtx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)

  ## ## Rearrangement so they are actually a distribution function
  ## Fyx <- lapply(Fyx, quantreg::rearrange)
  ## Ftx <- lapply(Ftx, quantreg::rearrange)

  ## ## Also, get their density (might want to change how we do this)
  ## fyx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")
  ## ftx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")

  if (messages) cat("Step 1 of 3: Converting QR to conditional density estimates...\n")
  ##fyx1 <- pbsapply(unique(data[,yname]), fy.x, betmat=t(coef(Qyx)), XX=x, tau=tau)
  fyx1 <- fYXmatC(Y=unique(data[,yname]), betmat=t(coef(Qyx)), X=x, tau=tau)
  eps <- 1e-300 ##.Machine$double.eps
  fyx <- apply(fyx1, 1, FUN=function(y) approxfun(x=unique(data[,yname]), y=y, yleft=eps, yright=eps))
  ##ftx1 <- pbsapply(unique(data[,tname]), fy.x, betmat=t(coef(Qtx)), XX=x, tau=tau)
  ftx1 <- fYXmatC(Y=unique(data[,tname]), betmat=t(coef(Qtx)), X=x, tau=tau)
  ftx <- apply(ftx1, 1, FUN=function(y) approxfun(x=unique(data[,tname]), y=y, yleft=eps, yright=eps))

  
  ## make draws from the mixture distribution
  Usig <- Qyx$sig
  Upi <- Qyx$pi
  Umu <- Qyx$mu
  Vsig <- Qtx$sig
  Vpi <- Qtx$pi
  Vmu <- Qtx$mu
  ksig <- length(Usig)
  Ucomponents <- sample(1:length(Usig), ndraws, replace=TRUE, prob=Upi)
  Us <- rnorm(ndraws, Umu[Ucomponents], Usig[Ucomponents])
  Vcomponents <- sample(1:length(Vsig), ndraws, replace=TRUE, prob=Vpi)
  Vs <- rnorm(ndraws, Vmu[Vcomponents], Vsig[Vcomponents])


  ################################################################
  if (messages) cat("\nStep 2 of 3: Estimating copula parameter...\n")
  ################################################################

  ## this is not right, but perhaps can make draws (e.g. similar to mcmc algorithm above)
  ## to get the draws right and then estimate this way.
  ## create a new dataset with the measurement error draws in order
  ## to estimate the copula parameter
  ## newids <- unlist(lapply(1:n, function(i) rep(i, ndraws))) ## just replicates Y and X over and over
  ## Yvals <- data[,yname]
  ## Tvals <- data[,tname]
  ## newdta1 <- cbind(Y=(Yvals[newids]-rep(Us,n)), T=(Tvals[newids]-rep(Vs,n)))

  if (copula=="frank") {
    cop <- copula::frankCopula()
  } else if (copula=="gumbel") {
    cop <- copula::gumbelCopula()
  } else if (copula=="clayton") {
    cop <- copula::claytonCopula()
  } else if (copula=="gaussian") {
    cop <- copula::normalCopula()
  } else {
    stop(paste0("copula type:", copula, " is not supported"))
  }

  ##this is not right, but perhaps can make draws (e.g. similar to mcmc algorithm above)
  ##to get the draws right and then estimate this way.
  #newdta1 <- data.frame(Y=Qyx$Ystar, T=Qtx$Ystar)
  #ranks1 <- copula::pobs(newdta1)
  ##cop <- copula::fitCopula(cop, ranks1, method="irho") ## irho inverts spearman's rho; it
  ##is very fast though (I think) not all copulas (exception=(I think)Gumbel) have 1-1 relationship
  ##with Spearman's rho, but in practice they seem very similar.


  
  ##delt <- rep(attributes(cop)$estimate, nrow(x))

  # estimation with maximum likelihood 
  res <- optimize(ll, c(0,1), maximum=TRUE, 
                  y=data[,yname], t=data[,tname], x=x, copula=copula,
                  Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx,
                  Us=Us, Vs=Vs)

  delt <- rep(parms2coppar(res$maximum, copula=copula, x=1), nrow(x))


  # estimate copula-type parameters
  cop <- copula::setTheta(cop, delt[1])
  Ystar_Tstar_inner <- lapply(1:nrow(x), function(i) {
    cop_draws <- copula::rCopula(100, cop)
    eY <- cop_draws[,1]
    eT <- cop_draws[,2]
    Y_xb <- as.numeric(t(as.matrix(x[i,]))%*%as.matrix(coef(Qyx)))
    T_xb <- as.numeric(t(as.matrix(x[i,]))%*%as.matrix(coef(Qtx)))
    Ystar <- sapply(eY, function(ey_draw) {
      interpolateC(tau, Y_xb, ey_draw, TRUE)
    })
    Tstar <- sapply(eT, function(et_draw) {
      interpolateC(tau, T_xb, et_draw, TRUE)
    })
    cbind(Ystar,Tstar)
  })

  Ystar_Tstar <- do.call("rbind", Ystar_Tstar_inner)
  Ystar <- Ystar_Tstar[,1]
  Tstar <- Ystar_Tstar[,2]

  t_mat <- tmat(Ystar, Tstar)
  Ps <- cor(Ystar, Tstar, method="spearman")
  up_mob <- upMob(Ystar, Tstar)
  
  # compute conditional distributions if values of the treatment are specified
  if (!is.null(tvals)) {
    
    if (messages) cat("\nStep 3 of 3: Building conditional distributions...\n")
    ## If you don't set particular values of X to compute,
    ## just set it equal to the average values of X in the dataset
    ##if (is.null(xdf)) xdf <- as.data.frame(t(apply(x,2,mean))) ##x, for all data

    yvals <- quantile(data[,yname], seq(.01,.99,.01)) ## could also take all unique yvals or let user pass them all in
    yvals <- yvals[order(yvals)]

    
    if (is.null(xdf)) xdf <- x
    
    if (copula == "gaussian") {
      rho <- delt[1]
      FytXmat <- array(dim=c(nrow(xdf), length(yvals), length(tvals)))
      for (j in 1:length(tvals)) {
        tt <- tvals[j]
        this.Fytx <- lapply(1:nrow(xdf), function(i) {
          pnorm ( ( qnorm(Fyx[[i]](yvals)) - rho*qnorm(Ftx[[i]](tt)) ) / ( sqrt(1-rho^2) ) )
        })
        FytXmat[,,j] <- do.call("rbind", this.Fytx)
      }
    }

    if (copula=="clayton") {
      thet <- delt[1]
      FytXmat <- array(dim=c(nrow(xdf), length(yvals), length(tvals)))
      C2 <- function(u,v,thet) {
        (u^(-thet) + v^(-thet) - 1)^(-(1/thet)-1) * v^(-(thet-1))
      }
      for (j in 1:length(tvals)) {
        tt <- tvals[j]
        this.Fytx <- lapply(1:nrow(xdf), function(i) {
          C2(Fyx[[i]](yvals), Ftx[[i]](tt), thet)
        })
        FytXmat[,,j] <- do.call("rbind", this.Fytx)
      }
    }

    if (copula=="gumbel") {
      thet <- delt[1]
      FytXmat <- array(dim=c(nrow(xdf), length(yvals), length(tvals)))
      dgenerator_gumbel <- function(t, thet) {
        d <- (-thet/t)*(-log(t))^(thet-1)
        d <- sapply(d, function(dd) max(dd,-100000000))
        d
      }
      ## C2 <- function(u,v,thet) {
      ##   inside <- (-log(u))^thet + (-log(v))^thet
      ##   exp( - (inside^(1/thet)) ) * (-(1/thet) * inside^( (1/thet) - 1)) * (-thet/v) * ( (log(v))^(thet-1) )
      ## }
      for (j in 1:length(tvals)) {
        tt <- tvals[j]
        this.Fytx <- lapply(1:nrow(xdf), function(i) {
          num <- dgenerator_gumbel(Ftx[[i]](tt), thet)
          denom <- dgenerator_gumbel( gumbel::invphigumbel( gumbel::phigumbel(Ftx[[i]](tt), thet) + gumbel::phigumbel(Fyx[[i]](yvals),thet), thet) , thet )
          num/denom
        })
        FytXmat[,,j] <- do.call("rbind", this.Fytx)
      }
    }

    #-----------------------------------------------------------------------------
    # this is old way of computing conditional distributions given copula estimate
    # but now we just directly compute them instead of doing it numerically
    #-----------------------------------------------------------------------------
   
    ## ##xtdf <- cbind(tvals, xdf)
    ## ## todo, this gives the copula, now convert to conditional distribution
    ## ##delt <- parms2coppar(res$par, copula, xdf)
    ## ##tvals <- quantile(data[,tname], tau, type=1)
    ## yvals <- quantile(data[,yname], seq(.01,.99,.01)) ## could also take all unique yvals or let user pass them all in
    ## yvals <- yvals[order(yvals)]
    ## QQyx <- predict( Qyx, newdata=as.data.frame(rbind(xdf,x[1,])), stepfun=TRUE)  ## super hack:  but predict.rqs is throwing an error that I think it shouldn't, and this gets around it.
    ## QQyx <- QQyx[-length(QQyx)]
    ## QQyx <- lapply(QQyx, rearrange)
    ## QQyx2  <- predict(Qyx, newdata=as.data.frame(xdf))
    ## FFyx2  <- t(sapply(Fyx, function(fyx) fyx(yvals)))
    ## Ftx <- predict(Qtx, newdata=as.data.frame(rbind(xdf, x[1,])),
    ##                type="Fhat", stepfun=TRUE)
    ## Ftx <- Ftx[-length(Ftx)]
    ## Ftx <- lapply(Ftx, rearrange)
    ## FFtx <- t(sapply(Ftx, function(ftx) ftx(tvals)))

    ## U <- seq(0,1,length.out=ndraws)
    ## U <- tau

    ## if (copula=="frank") {
    ##   cop <- copula::frankCopula(as.numeric(delt[1])) ## all delts restricted to be the same so just choose first one
    ## } else if (copula=="gumbel") {
    ##   cop <- copula::gumbelCopula(as.numeric(delt[1]))
    ## } else if (copula=="clayton") {
    ##   cop <- copula::claytonCopula(as.numeric(delt[1]))
    ## } else if (copula=="gaussian") {
    ##   cop <- copula::normalCopula(as.numeric(delt[1]))
    ## } else {
    ##   stop(paste0("copula type:", copula, " is not supported"))
    ## }

    ## ## call C++ function to return matrix with distribution of Fytx
    ## FytXmat <- computeFytXC(yvals, tvals, QQyx2, FFtx, tau, "gumbel", delt[1])

    
    ## internal function for reordering arguments of BMisc::makeDist
    makeDist1 <- function(Fx, x, sorted = FALSE, rearrange=FALSE) {
      BMisc::makeDist(x, Fx, sorted, rearrange)
    }

    # note to self: can run into common support issues if not careful in simulations here.
   
    ## converts 3-dimensional matrix of Fytx into 2-dimensional matrix of distributio
    ## functions
    Fytx  <- apply(FytXmat, c(1,3), makeDist1, x=yvals)

    ## FytX contains a distribution function for every value of t and x
    ## this step averages over all the x's
    Fyt <- lapply(1:length(tvals), function(i) {
      BMisc::combineDfs(yvals, Fytx[,i], rearrange=TRUE)
    })
    

    ## old way, replaced this with calls to C++ functions
    ## Fytx <- pblapply(1:length(QQyx), function(i) {
    ##     qfun <- QQyx[[i]]
    ##     ffun <- Ftx[[i]]
    ##     lapply(tvals, function(tt) {
    ##         BMisc::makeDist(yvals, sapply(yvals, function(yy) {
    ##             mean(1*(qfun(U)<=yy)*dCopula(cbind(U, ffun(tt)), cop))
    ##         }))
    ##     }) ## might want to make this a function with makeRQS (modified) or makeDist eventually
    ## }, cl=cl)
    ## ##Qytx <- lapply(Fytx, function(FF) quantile(FF, tau, type=1))

    ## ## next we want to reverse the list
    ## reverseListIndex <- function(l) {
    ##     outlist <- list()
    ##     length(outlist) <- length(l[[1]])
    ##     outlist <- lapply(outlist, function(f) {
    ##         g <- list()
    ##         length(g) <- length(l)
    ##         g
    ##     })
    
    ##     for (i in 1:length(l)) {
    ##         for (j in 1:length(l[[1]])) {
    ##             outlist[[j]][[i]] <- l[[i]][[j]]
    ##         }
    ##     }
    ##     outlist
    ## }

    ## Fytx <- reverseListIndex(Fytx)
    

    ##Fyt <- lapply(Fytx, function(FFytx) {
    ##    combineDfs(yvals, FFytx)
    ##})

    ## out <- list(cop.param=parms2coppar(res$maximum, copula=copula, x=1),
    ##             copula=copula, Fytxlist=Fytx, Fyt=Fyt, tvals=tvals, x=xdf)

    if (!retFytxlist) { ## often want to drop this because it is huge
      Fytxlist <- NULL
    }
    out <- list(cop.param=delt[1], copula=copula, Fytxlist=Fytx, Fyt=Fyt, tvals=tvals, x=xdf,
                t_mat=t_mat, Ps=Ps, up_mob=up_mob)

    ### only do above if you want the results for a particular value of t and x;
    ### otherwise can just return all results 
  } else {
    ## out <- list(cop.param=parms2coppar(res$maximum, copula=copula, x=x),
    ##             copula=copula)
    out <- list(cop.param=delt[1], copula=copula, t_mat=t_mat, Ps=Ps, up_mob=up_mob)
  }

  out$Qyx <- Qyx
  out$Qtx <- Qtx

  class(out) <- "qr2meobj"

  return(out)
}

#########################################################
##
## Additional helper functions for working with and
## manipulating main results
##
#########################################################

## check if this works, I think it is for putting back together an
## unconditional distribution from a list of them, but need to step through
## it
## avgDist <- function(Fytxlist, yvals) {
  
##   Fyt <- lapply(Fytxlist, function(FFytx) {
##     combineDfs(yvals, FFytx)
##   })

##   Fyt
## }


#' qr2meobj
#'
#' class for qr2meobj
#'
#' @param cop.param copula parameter
#' @param copula type of copula
#' @param tvals values of the treatment that conditional distributions were
#'  estimated for
#' @param x matrix of covariates
#' @param Fytxlist list of conditional distributions
#' @param Qyx estimates of quantiles of Y conditional on X
#' @param Qtx estimates of quantiles of T conditional on X
#'
#' @export
qr2meobj <- function(cop.param, copula, tvals, x, Fytxlist, Qyx, Qtx) {
  out <- list()
  out$cop.param <- cop.param
  out$copula <- copula
  out$tvals <- tvals
  out$x <- x
  out$Fytxlist <- Fytxlist
  out$Qyx <- Qyx
  out$Qtx <- Qtx

  class(out) <- "qr2meobj"
}

#' print.qr2meobj
#'
#' @param x a qr2meobj
#' @param ... other arguments
#' 
#' @export
print.qr2meobj <- function(x, ...) {
  print(x$Qyx)
  print(x$Qtx)
  cat("\n\n")
  cat("Copula Type: ")
  cat(x$copula)
  cat("\n")
  cat("Copula Paramter: ")
  cat(x$cop.param)
  cat("\n\n")
}

#' print.merr
#'
#' @param x an merr object
#' @param ... other arguments
#' 
#' @export
print.merr <- function(x,...) {
  coef <- round(t(as.matrix(x$coefficients)),4)
  rownames(coef) <- x$tau
  colnames(coef) <- c("Intercept", attr(x$terms,"term.labels"))
  cat("Coefficients:\n")
  print(coef)
  cat("\n\n")
  cat("Measurement Error Distribution:\n")
  U <- round(cbind(x$pi, x$mu, x$sig^2),4)
  rownames(U) <- sapply(1:nrow(U), function(i) paste0("Comp. ",i))
  colnames(U) <- c("Prob.", "Mean", "Variance")
  print(U)    
}

## getListElement <- function(listolists, whichone=1) {
##   lapply(listolists, function(l) l[[whichone]])
## }

#' gg_qr2meobj
#'
#' plot a qr2meobj
#'
#' @param obj a qr2meobj
#' @param tau which quantiles to plot
#' @param ylim limits of y-axis
#' @param ylab label for y-axis
#' @param xlab label for x-axis
#'
#' @export
gg_qr2meobj <- function(obj, tau=c(.1,.5,.9), ylim=NULL,
                          ylab=NULL, xlab=NULL) {

  if (is.null(tau)) tau <- c(.1,.5,.9)
  qq <- t(sapply(obj$Fyt, function(FF) quantile(FF, probs=tau)))
  tvals <- obj$tvals
  
  cmat <- cbind.data.frame(tvals, qq)
  colnames(cmat) <- c("tvals", paste0("c",tau*100))
  cmat <- gather(cmat, quantile, value, -tvals)

  p <- ggplot(cmat, mapping=aes(x=tvals,y=value, group=quantile, color=quantile)) +
    geom_line() +
    geom_point()

  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits=ylim)
  }

  if (!is.null(ylab)) {
    p <- p + ggplot2::ylab(ylab)
  }

  if (!is.null(xlab)) {
    p <- p + ggplot2::xlab(xlab)
  }

  p
}

#' addplot
#'
#' @param obj new object to plot
#' @param p existing plot
#' @param tau which quantiles to plot
#'
#' @export
addplot <- function(obj, p, tau=c(.1,.5,.9)) {

  qq <- t(sapply(obj$Fyt, function(FF) quantile(FF, probs=tau)))

  cmat <- cbind.data.frame(tvals, qq)
  colnames(cmat) <- c("tvals", paste0("c",tau*100))
  cmat <- tidyr::gather(cmat, quantile, value, -tvals)
  p <- p + geom_line(data=cmat, mapping=aes(x=tvals, y=value, group=quantile,
                                            color=quantile),
                     linetype="dashed")
  p <- p + geom_point(data=cmat, mapping=aes(x=tvals, y=value, group=quantile,
                                             color=quantile))
  p
}



