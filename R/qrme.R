# =============================================================================
# Title: QRME Main Functions
# Description: compute.qrme() and qrme() implement quantile regression with
#              mixture-of-normals measurement error in the dependent variable
#              via the pseudo-EM algorithm in em.R. Also includes qr2me() for
#              two-sided measurement error and supporting S3 methods.
# Author: Brant Callaway
# Last update: 2026-05-14
# Date created: 2026-05-07
# =============================================================================

#' @title compute.qrme
#' @description does the heavy lifting on computing quantile regression with
#'  left hand side measurement error
#'
#' @inheritParams qrme
#'
#' @keywords internal
#' @export
compute.qrme <- function(formula, tau=0.5, data, me_dist="gaussian", nmix=1, start_beta=NULL, start_mu=NULL,
                         start_sigma=NULL, start_pi=NULL, mcmc_method="MH", tol=NULL, conv_crit="params",
                         ndraws_ll=100L, conv_patience=1L, mcmc_draws=200, mcmc_burnin=100, proposal_sd=NULL, ncores=1,
                         maxit=100, verbose=FALSE) {
  xformula <- formula
  xformula[[2]] <- NULL ## drop y variable
  x <- model.matrix(xformula, data)
  yname <- as.character(formula[[2]])
  y <- data[,yname]
  k <- ncol(x) ## number of x variables
  m <- nmix
  if (is.null(start_beta)) {
    ## this defaults the start values of the beta_0 to be
    ## the observed quantiles of the outcome and the other
    ## betas to be equal to 0
    betvals <- t(coef(quantreg::rq(formula, tau=tau, data=data)))

  } else {
    betvals <- start_beta
  }
  # Target ME mixture SD: half the outcome SD puts the ME on the right scale
  # without assuming it dominates the signal.
  target <- 0.5 * sd(y)

  if (is.null(start_mu)) {
    if (m == 1L) {
      muvals <- 0
    } else {
      # Spread component means evenly over [-target/2, target/2] and center.
      # Spanning ±target/2 keeps the means within the variance budget so that
      # sigma can always be backed out as a meaningful positive value.
      muvals <- seq(-target / 2, target / 2, length.out = m)
      muvals <- muvals - mean(muvals)
    }
  } else {
    muvals <- start_mu
  }
  if (is.null(start_sigma)) {
    # Back out sigma from the mixture variance identity with equal weights and
    # equal component SDs:
    #   Var(mixture) = sigma^2 + mean(muvals^2)
    # Setting Var(mixture) = target^2 gives sigma = sqrt(target^2 - mean(muvals^2)).
    # Floor at target/sqrt(m): the equal-share lower bound if means dominate.
    var_means <- mean(muvals^2)
    sigma <- sqrt(max(target^2 - var_means, (target / sqrt(m))^2))
    sigvals <- rep(sigma, m)
  } else {
    sigvals <- start_sigma
  }
  if (is.null(start_pi)) {
    pivals <- rep(1/m, m)
  } else {
    pivals <- start_pi
  }
  qrparams <- list(bet=betvals, m=m, pi=pivals, mu=muvals, sig=sigvals)


  ## Estimate QR model with measurement error using EM algorithm
  res <- em.algo(formula, data,
                 betmatguess=betvals, tau=tau,
                 me_dist=me_dist,
                 m=m, piguess=pivals, muguess=muvals,
                 sigguess=sigvals, mcmc_method=mcmc_method, tol=tol,
                 conv_crit=conv_crit, ndraws_ll=ndraws_ll,
                 conv_patience=conv_patience,
                 mcmc_draws=mcmc_draws, mcmc_burnin=mcmc_burnin, proposal_sd=proposal_sd, ncores=ncores,
                 maxit=maxit, verbose=verbose)


  out <- makeRQS(res, formula, data, tau=tau)

  class(out) <- c("merr", class(out))

  ## set final parameters for measurement error class
  out$bet <- res$bet
  out$pi <- res$pi
  out$mu <- res$mu
  out$sig <- res$sig
  out$x <- x
  out$y <- y
  out$me_dist <- me_dist
  out$n_iter <- res$n_iter
  out$tol <- res$tol
  out$conv_crit <- res$conv_crit
  out$conv_criteria <- res$conv_criteria
  out$conv_converged <- res$conv_converged
  out$mix_n_iter <- res$mix_n_iter
  out$mix_loglik <- res$mix_loglik
  out$mix_converged <- res$mix_converged
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
#' @param formula y ~ x
#' @param tau vector for which quantiles to compute quantile regression
#' @param data a data.frame that contains y and x
#' @param nmix The number of mixture components of the measurement error
#'  (default 1). Increase this for more flexible measurement-error
#'  distributions.
#' @param start_beta an LxK matrix of starting values for beta where
#'  L is the dimension of tau and K is the number of covariates (default is
#'  NULL and in this case, the starting values are set to be the QR
#'  coefficients coming from QR that ignores measurment error)
#' @param start_mu A vector of length nmix of starting values for the mean
#'  of each mixture component. When \code{NULL} (default), set to 0 for
#'  \code{nmix = 1} and to evenly-spaced values on
#'  \code{[-0.25 * sd(y), 0.25 * sd(y)]} for \code{nmix > 1}, mean-centered
#'  to satisfy the mean-zero ME constraint.
#' @param start_sigma A vector of length nmix of starting values for the
#'  standard deviation of each mixture component. When \code{NULL} (default),
#'  backed out from the mixture variance identity assuming equal weights and
#'  equal component SDs targeting a mixture SD of \code{0.5 * sd(y)}:
#'  \code{sigma = sqrt(max(target^2 - mean(start_mu^2), (target/sqrt(nmix))^2))}.
#'  This gives \code{0.5 * sd(y)} for \code{nmix = 1} and slightly larger
#'  values for \code{nmix > 1} (since the component means absorb some variance).
#' @param start_pi A vector of length nmix of starting values for the fraction
#'  of observations in each component of the mixture of normals distribution
#'  for the measurement error (default is NULL and in this case, the starting
#'  values are all set to be 1/nmix)
#' @param mcmc_method The type of simulation step to use in the EM algorithm.
#'  Currently only \code{"MH"} for Metropolis-Hastings is supported.
#' @param tol Convergence tolerance.  When \code{NULL} (default), a value is
#'  chosen automatically based on \code{conv_crit}.  See \code{\link{em.algo}}
#'  for details.
#' @param conv_crit Convergence criterion: \code{"params"} (default) uses the
#'  Euclidean norm of parameter changes; \code{"loglik"} uses the relative
#'  change in the observed-data log-likelihood.  See \code{\link{em.algo}}.
#' @param ndraws_ll Number of Monte Carlo draws for the log-likelihood
#'  convergence check when \code{conv_crit = "loglik"} (default 100).
#'  Ignored when \code{conv_crit = "params"}.
#' @param conv_patience Integer. Consecutive iterations that must satisfy the
#'  convergence criterion before stopping (default 1). Setting to 2 guards
#'  against false convergence from Monte Carlo noise. See \code{\link{em.algo}}.
#' @param mcmc_draws Total number of MCMC draws per EM step (default 200)
#' @param mcmc_burnin Number of MCMC draws to discard as burnin (default 100)
#' @param proposal_sd Standard deviation of the Metropolis-Hastings proposal
#'  (random-walk step size). When \code{NULL} (default), initialised to the
#'  SD of the ME mixture from the starting parameters and updated after each
#'  M-step to track the current ME scale. Pass a positive numeric to fix the
#'  proposal at that value for all iterations.
#' @param ncores Number of cores for parallel bootstrap computation (default 1)
#' @param maxit Maximum number of EM outer iterations. If convergence is not
#'  reached, the estimates from the final iteration are returned (default is 100)
#' @param se Whether or not to compute standard errors using the bootstrap
#'  (default is FALSE)
#' @param n_boot Number of bootstrap iterations for standard errors (default 100)
#' @param verbose Logical or nonnegative integer. If \code{FALSE} (default),
#'  suppresses progress output. If \code{TRUE} or \code{1}, reports major
#'  computational stages, EM iteration numbers, and convergence diagnostics.
#'  If \code{2}, also reports qrme-native details such as \code{pi},
#'  \code{mu}, \code{sigma}, tolerance, and finite-mixture summaries. If
#'  \code{3} or larger, also reports raw \code{normalmixEM()} output, which
#'  uses \code{lambda} for mixture probabilities.
#' 
#' @return an object of class "merr". Supports \code{logLik()}, \code{AIC()},
#'   and \code{BIC()} for comparing fits across different starting values.
#'
#' @export
qrme <- function(formula, tau=0.5, data, me_dist="gaussian", nmix=1, start_beta=NULL, start_mu=NULL,
                 start_sigma=NULL, start_pi=NULL, mcmc_method="MH", tol=NULL, conv_crit="params",
                 ndraws_ll=100L, conv_patience=1L, mcmc_draws=200, mcmc_burnin=100, proposal_sd=NULL, ncores=1,
                 maxit=100, se=FALSE, n_boot=100, verbose=FALSE) {



  res <- compute.qrme(formula=formula,
                      tau=tau,
                      data=data,
                      me_dist=me_dist,
                      nmix=nmix,
                      start_beta=start_beta,
                      start_mu=start_mu,
                      start_sigma=start_sigma,
                      start_pi=start_pi,
                      mcmc_method=mcmc_method,
                      tol=tol,
                      conv_crit=conv_crit,
                      ndraws_ll=ndraws_ll,
                      conv_patience=conv_patience,
                      mcmc_draws=mcmc_draws,
                      mcmc_burnin=mcmc_burnin,
                      proposal_sd=proposal_sd,
                      ncores=ncores,
                      maxit=maxit,
                      verbose=verbose)

  if (se) {
    eachIter <- pbapply::pblapply(1:n_boot, function(b) {
      n <- nrow(data)
      brows <- sample(1:n, size=n, replace=TRUE)
      bdata <- data[brows,]
      tryCatch({
        out <- compute.qrme(formula=formula,
                            tau=tau,
                            data=bdata,
                            me_dist=me_dist,
                            nmix=nmix,
                            start_beta=res$bet,
                            start_mu=res$mu,
                            start_sigma=res$sig,
                            start_pi=res$pi,
                            mcmc_method=mcmc_method,
                            tol=tol,
                            conv_crit=conv_crit,
                            ndraws_ll=ndraws_ll,
                            conv_patience=conv_patience,
                            mcmc_draws=mcmc_draws,
                              mcmc_burnin=mcmc_burnin,
                              proposal_sd=proposal_sd,
                              ncores=1,
                              maxit=maxit)
        out$Ystar=NULL ## just drop this because it takes up a lot of memory
        out
      }, error=function(cond) {
        return(NULL)
      })
    }, cl=ncores)

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
#' @param tvals a vector of values of the treatment variable at which to compute the conditional distribution 
#'  of Y given X and T
#' @param me_dist which type of measurement error distribution to use (default is "gaussian"), 
#'  "laplace" is also supported
#' @param copula which type of copula to use (default is "gaussian")
#' @param Qyx quantile regression estimates (can be adjusted for measurement
#'  error) of Y on X
#' @param Qtx quantile regression estimates (can be adjusted for measurement
#'  error) of T on X
#' @param retFytxlist whether or not to return the conditional distribution
#'  for every value of x in xdf
#'  (default is FALSE because this can take up a lot of room in memory)
#' @param copula_me_draws Number of measurement-error draws used in the
#'  second-stage copula likelihood (default 100). Increase this for a less
#'  noisy copula likelihood at the cost of speed.
#' @param mobility_copula_draws Number of copula draws per covariate row used
#'  to compute mobility summaries such as transition matrices and rank
#'  correlations (default 100).
#' @inheritParams tsme
#' @inheritParams qrme
qr2me <- function(yname, tname, xformla, tau, data, xdf=NULL, tvals=NULL,
                    me_dist="gaussian", copula="gaussian",
                    Qyx, Qtx, retFytxlist=FALSE,
                    copula_me_draws=100L, mobility_copula_draws=100L,
                    verbose=FALSE) {

  if (qrme_verbose_level(verbose) >= 1L) {
    message("qr2me method")
    message("----------------------")
    message(
      "Citation: Callaway, Brantly, Tong Li, Irina Murtazashvili, ",
      "and Emmanuel Tsyawo, Distributional Effects with Two-Sided ",
      "Measurement Error: An Application to Intergenerational Income ",
      "Mobility, Working Paper, 2025."
    )
    message("----------------------")
  }
  
  x <- model.matrix(xformla, data)
  n <-  nrow(data)

  # Internal numerical grid for turning sparse QR estimates into CDF objects.
  tau_grid <- seq(0, 1, length.out = 100)
  Qyx_interpolated <- lapply(1:nrow(x), function(i) {
    xb <- as.numeric(t(as.matrix(x[i,]))%*%as.matrix(coef(Qyx)))
    sapply(tau_grid, function(tt_grid_val) {
      interpolateC(tau, xb, tt_grid_val, TRUE)
    })
  })
  Fyx <- lapply(Qyx_interpolated, function(Qyx) BMisc::makeDist(Qyx, tau_grid, sorted=TRUE))


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

  qrme_progress(verbose, "Step 1 of 3: converting QR to conditional density estimates...")
  ##fyx1 <- pbsapply(unique(data[,yname]), fy.x, betmat=t(coef(Qyx)), XX=x, tau=tau)
  fyx1 <- fYXmatC(Y=unique(data[,yname]), betmat=t(coef(Qyx)), X=x, tau=tau)
  eps <- .Machine$double.eps
  fyx <- apply(fyx1, 1, FUN=function(y) approxfun(x=unique(data[,yname]), y=y, yleft=eps, yright=eps))
  ##ftx1 <- pbsapply(unique(data[,tname]), fy.x, betmat=t(coef(Qtx)), XX=x, tau=tau)
  ftx1 <- fYXmatC(Y=unique(data[,tname]), betmat=t(coef(Qtx)), X=x, tau=tau)
  ftx <- apply(ftx1, 1, FUN=function(y) approxfun(x=unique(data[,tname]), y=y, yleft=eps, yright=eps))

  
  ## make draws from the measurement error distribution
  if (me_dist == "laplace") {
    Us <- rlaplace(copula_me_draws, mu=0, sigma=Qyx$sig)
    Vs <- rlaplace(copula_me_draws, mu=0, sigma=Qtx$sig)
  } else { # gaussian
    Usig <- Qyx$sig
    Upi <- Qyx$pi
    Umu <- Qyx$mu
    Vsig <- Qtx$sig
    Vpi <- Qtx$pi
    Vmu <- Qtx$mu
    ksig <- length(Usig)
    Ucomponents <- sample(1:length(Usig), copula_me_draws, replace=TRUE, prob=Upi)
    Us <- rnorm(copula_me_draws, Umu[Ucomponents], Usig[Ucomponents])
    Vcomponents <- sample(1:length(Vsig), copula_me_draws, replace=TRUE, prob=Vpi)
    Vs <- rnorm(copula_me_draws, Vmu[Vcomponents], Vsig[Vcomponents])
  }

  ################################################################
  qrme_progress(verbose, "Step 2 of 3: estimating copula parameter...")
  ################################################################

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

  # estimation with maximum likelihood 
  res <- optimize(ll, c(0,1), maximum=TRUE, 
                  y=data[,yname], t=data[,tname], x=x, copula=copula,
                  Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx,
                  Us=Us, Vs=Vs)

  delt <- rep(parms2coppar(res$maximum, copula=copula, x=1), nrow(x))


  # estimate copula-type parameters
  cop <- copula::setTheta(cop, delt[1])
  Ystar_Tstar_inner <- lapply(1:nrow(x), function(i) {
    cop_draws <- copula::rCopula(mobility_copula_draws, cop)
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
    
    qrme_progress(verbose, "Step 3 of 3: building conditional distributions...")
    ## If you don't set particular values of X to compute,
    ## just set it equal to the average values of X in the dataset
    ##if (is.null(xdf)) xdf <- as.data.frame(t(apply(x,2,mean))) ##x, for all data

    yvals <- unique(data[,yname])#quantile(data[,yname], seq(.01,.99,.005)) ## could also take all unique yvals or let user pass them all in
    yvals <- unique(yvals[order(yvals)])


    # only re-calculate if xdf is set
    if (!is.null(xdf)) {
      
      Qyxdf_interpolated <- lapply(1:nrow(xdf), function(i) {
        xb <- as.numeric(as.matrix(xdf[i,])%*%as.matrix(coef(Qyx)))
        sapply(tau_grid, function(tt_grid_val) {
          interpolateC(tau, xb, tt_grid_val, TRUE)
        })
      })
      Fyx <- lapply(Qyxdf_interpolated, function(Qyx) BMisc::makeDist(Qyx, tau_grid, sorted=TRUE))


      Qtxdf_interpolated <- lapply(1:nrow(xdf), function(i) {
        xb <- as.numeric(as.matrix(xdf[i,])%*%as.matrix(coef(Qtx)))
        sapply(tau_grid, function(tt_grid_val) {
          interpolateC(tau, xb, tt_grid_val, TRUE)
        })
      })
      Ftx <- lapply(Qtxdf_interpolated, function(Qtx) BMisc::makeDist(Qtx, tau_grid, sorted=TRUE))
    }
    
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
        (u^(-thet) + v^(-thet) - 1)^(-(1/thet)-1) * v^(-thet-1)
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
