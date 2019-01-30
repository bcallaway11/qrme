#' @title cop.pdf
#' @description returns the pdf of several copula functions evaluated at
#'  u and v; the available copula functions are:  Gumbel (the default),
#'  Frank, (should implement Joe and Gaussian)
#'
#' @param u first copula argument (u and v can be vectors)
#' @param v second copula argument
#' @param type one of c("gumbel","frank")
#' @param delt the copula paramter
#'
#' @return vector of copula pdf values
#'
#' @export
cop.pdf <- function(u,v,type="gumbel",delt,eps=1e-300) {
    if ( !(all(0 < u) & all(u < 1)) | !(all(0 < v) & all(v < 1))) {
        stop("u and v must be between 0 and 1")
    }
    if (type == "gumbel") {
        ## lu <- -log(u)
        ## lv <- -log(v)
        ## out <- exp(-(lu^delt + lv^delt)^(1/delt)) * ( (lu^delt + lv^delt)^(1/delt) + delt - 1) * (lu^delt + lv^delt)^(1/delt - 2) * (lu*lv)^(delt-1) * (u*v)^(-1)
        gc <- copula::gumbelCopula(delt, use.indepC="TRUE")
    } else if (type=="frank") {
        gc <- copula::frankCopula(delt, use.indepC="TRUE")
    } else if (type=="clayton") {
        gc <- copula::claytonCopula(delt, use.indepC="TRUE")
    } else if (type=="gaussian") {
        gc <- copula::normalCopula(delt)
    } else {
        stop(paste0("copula: ", type, " not supported"))
    }
    out <- copula::dCopula(cbind(u,v), gc)
    out <- sapply(out, function(o) max(o,eps))
    return(out)
}

#' @title ll 
#' @description this is the likelihood function for estimating the copula parameter
#' @param params a value of the parameters to calculate the log likelihood for
#' @param y should be a vector of outcomes
#' @param t should be a vector of treatments
#' @param x a matrix of values of x
#' @param copula the type of copula, should be supported by cop.pdf
#' @param Fyx function that contains the cdf of y conditional on x
#' @param Ftx function that contains the cdf of t conditional on x
#' @param fyx function that contains the pdf of y conditional on x
#' @param ftx function that contains the pdf of t conditional on x
#' @param Upi vector of probabilities for measurement error mixture model for Y
#' @param Umu vector of means for measurement error mixture model for Y
#' @param Usig vector of std. deviations for measurement error mixture model for Y
#' @param Vpi vector of probabilities for  measurement error mixture model for T
#' @param Vmu vector of means for measurement error mixture model for T
#' @param Vsig vector of std. deviations for measurement error mixture model for T
#' @param ndraws the number of draws of U and V to make
#'
#' @return scalar negative value of log likelihood function
#'
#' @export
ll <- function(params, y, t, x, copula="gumbel", Fyx, Ftx, fyx, ftx, Us, Vs, ndraws=100, eps=1e-300) {
    k <- length(params)
    params <- as.matrix(params)
    x <- as.matrix(x)
    n <- nrow(x)
    
    intonly <- as.matrix(rep(1, nrow(x)))
    delt <- parms2coppar(params, copula, intonly)
    ## if (copula == "gumbel") {
    ##     delt <- 1 + exp(x%*%params)
    ## } else if (copula == "frank") {
    ##     delt <- x%*%params
    ## } else if (copula == "gaussian") {
    ##     delt <- 2 * exp(x%*%params) / (1 + exp(x%*%params)) - 1
    ## }

    lval <- sapply(1:n, function(i) {
        max(
            mean(
            cop.pdf(Fyx[[i]](y[i] - Us), Ftx[[i]](t[i] - Vs), type=copula,
                delt=delt[i])  *
            fyx[[i]](y[i] - Us) * ftx[[i]](t[i] - Vs)
            )
          , eps) ## here, eps
        ## avoids this being exactly equal to 0
        })
    sum(log(lval))
}




#' @title qrme
#'
#' @description Quantile Regression with measurement error in the dependent
#'  variable using an EM-algorithm
#'
#' @param formla y ~ x
#' @param tau vector for which quantiles to compute quantile regression
qrme <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL, startsig=NULL, startpi=NULL, method="BFGS", maxit=1000, betOnly=FALSE) {
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
        betvals <- t(coef(rq(formla, tau=tau, data=data)))
        
    } else {
        betvals <- startbet
    }
    if (is.null(startmu)) {
        muvals <- seq(-1, by=1, length.out=(m-1))
    } else {
        muvals <- startmu
    }
    if (is.null(startsig)) {
        sigvals <- rep(1,m)
    } else {
        sigvals <- startsig
    }
    if (is.null(startpi)) {
        pivals <- rep(1/m, (m-1))
    } else {
        pivals <- startpi
    }
    qrparams <- list(bet=betvals, m=m, pi=pivals, mu=muvals, sig=sigvals)


    ## TODO: fill in details for estimting this using EM
    
    out <- makeRQS(params, formla, data)

    class(out) <- c("merr", class(out))
    
    out$params <- res$par
    ## idx <- length(betvals)+1
    ## nidx <- length(betvals)+length(pivals)
    ## pi1 <- res$par[idx:nidx]
    ## out$pi <- c(pi1, 1-sum(pi1))
    ## idx <- nidx+1
    ## nidx <- nidx+length(muvals)
    ## mu1 <- res$par[idx:nidx]
    ## out$mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    ## idx <- nidx+1
    ## nidx <- nidx+length(sigvals)
    ## out$sig <- res$par[idx:nidx]
    out$pi <- params$pi
    out$mu <- params$mu
    out$sig <- params$sig

    out
}


#' @title em.algo
#'
#' @description A pseudo EM algorithm for quantile regression with measurement error.  The measurement error here follows a mixture of normals.
#' @param betmatguess Initial values for the beta parameters.  This should be
#'  an LxK matrix where L is the number of quantiles and K is the dimension
#'  of the covariates
#' @param k The number of components of the mixture distribution for the
#'  measurement error
#' @param piguess Starting value for the probabilities of each mixture distribution (should have length equal to k)
#' @param muguess Starting value for the mean of each mixture component (should
#'  have length equal to k)
#' @param sigguess Starting value for the standard deviation of each mixture
#'  component (should have length equal to k)
#' @param tol This is the convergence criteria.  When the change in the
#'  Euclidean distance between the new parameters (at each iteration) and
#'  the old parameters (from the previous iteration) is smaller than tol,
#'  the algorithm concludes.  In general, larger values for tol will result
#'  in a fewer number of iterations and smaller values will result in more
#'  accurate estimates.
#' @return QRME object
#'
#' @keywords internal
#' @export
em.algo <- function(betmatguess, k=1, piguess=1, muguess=0, sigguess=1, tol=.01) {
  
  stopIters <- 100
  counter <- 1
  
  while (counter <= stopIters) {
    newone <- em.algo.inner(betmatguess, k, piguess, muguess, sigguess)
    newbet <- newone$bet
    newpi <- newone$pi
    newmu <- newone$mu
    newsig <- newone$sig
    
    cat("\n\n\nIteration: ", counter, "\n bet: ", newbet, "\n pi: ", newpi, "\n mu: ", newmu, "\n sig: ", newsig, "\n\n")
    criteria <- sqrt(sum(c(newbet-betmatguess,newsig-sigguess,newpi-piguess,newmu-muguess)^2))
    cat(" convergence criteria: ", criteria, "\n\n")
    if ( criteria <= tol) { ## Euclidean norm
      cat("\n algorithm converged\n") 
      return(newone)
    }
    counter <- counter+1
    betmatguess <- newbet
    sigguess <- newsig
    muguess <- newmu
    piguess <- newpi
  }
  cat("\n algorithm failed to converge\n")
  return(newone)
}

em.algo.inner <- function(betmat, k=1, pi=1, mu=0, sig=1) {
  
  cat("\nSimulating measurement error...")
  newdta <- pbapply::pblapply(1:n, function(i) {
    e <- mh_mcmc(betmat=bmat, k=k, pi=pi, mu=mu, sig=sig, y=YY[i], x=XX[i,], u.seq=tau, iters=1000, burnin=500)
    cbind(YY[i]-e, matrix(XX[i,], nrow=length(e), ncol=ncol(XX), byrow=TRUE),e=e)
  }, cl=10)
  
  ## Note: this currently just works for one X; will need to update
  
  cat("\nEstimating QR including simulated measurement error...")
  newdta1 <- do.call(rbind.data.frame, newdta)
  colnames(newdta1) <- c("Y", "int", "X", "e")
  out <- rq(Y ~ X, tau=tau, data=newdta1, method="pfn")
  
  ## this is part I am not sure about, once you have a new beta then estimate a new sigma??
  ## also should probably restrict overall mean of error term to be equal to 0
  cat("\nEstimating finite mixture model...")
  if (k == 1) {
    nm <- list(k=1, lambda=1, mu=0, sigma=sd(newdta1$e))
  } else {
    nm <- normalmixEM(newdta1$e, k=k)
  }
  
  return(list(bet=t(coef(out)), k=k, pi=nm$lambda, mu=nm$mu, sig=nm$sigma))
}



#' @title betfun
#'
#' @description Turns matrices of QR parameter values into functions.
#'  The returned function will be a map from (0,1) to a vector of dimension
#'  K where K is the number of elements in X (also it is the number of columns
#'  in the passed in matrix of QR parameter values).  It works by linear
#'  interpolation of parameter values (for interior values of tau) and
#'  based on some assumptions to handle values in the tails.
#' @param betmat An L x K matrix of parameter values where 
#'  L is the number of tau and 
#'  K is the dimension of X
#' @param tau The particular quantiles that the parameters were estimated at.
#'  It should be an L-dimensional vector.
#' @return A function that can be called on any element from 0 to 1 and
#'  returns a vector of parameter values
#' @export
betfun <- function(betmat, tau) {
  betmat <- as.matrix(betmat)
  betfun.list <- lapply(1:ncol(betmat), function(j) {
    if (j==1) { ##then we are on the constant
      betfun.inner(betmat[,j], tau, isconst=TRUE)
    } else {
      betfun.inner(betmat[,j], tau)
    }
  })
  bet <- function(u) {
    unlist(lapply(betfun.list, function(b) b(u)))
  }
  return(bet)
}

#' @title betfun.inner
#'
#' @description Does the heavy lifting for betfun.  Basically, betfun is just
#'  a wrapper for this that can handle a matrix of values of parameters.
#'  This function does the work, but only for a single vector of betas.
#' @param betvec vector of parameter values
#' @param tau corresponding vector of quantiles where beta was estimated
#' @param isconst ?? -- what is this doing?
#' @return function that takes argument from (0,1)
#' @keywords internal
#' @export
betfun.inner <- function(betvec, tau, isconst=FALSE) {
  bet <- function(u) {
    ul.pos <- tail(which(tau <= u),1) ## position of ul (in text)
    uu.pos <- which(tau >= u)[1] ## position of uu (in text)
    ul <- tau[ul.pos]
    uu <- tau[uu.pos]
    lam1 <- 1 - min(tau)
    lam2 <- max(tau)
    isconst <- 1*isconst
    if (u >= 0 & u < min(tau)){  ## this is case with u between 0 and smallest element of tau
      return(betvec[uu.pos] + isconst*log(u/min(tau))/lam1)
    }
    if (u <= 1 & u > max(tau)) { ## case with u < 1 but bigger than greatest element in tau
      return(betvec[ul.pos] + isconst*log((1-u)/(1-max(tau)))/lam2)  
    }
    if (is.na(uu) | is.na(ul)) {
      stop("uu or ul not set, likely tried some value of u that is not between the minimum and maximum values in tau")
    }
    if (uu == ul) {
      betvec[ul.pos]
    } else {
      betvec[ul.pos] + (u-ul)*(betvec[uu.pos] - betvec[ul.pos])/(uu-ul) ## this is linear interpolation from the text
    }
  }
  return(bet)
}

## X should be n x k
## the density of y conditional on x
fy.x <- function(y, betmat, X, tau) {
  X <- as.matrix(X)
  
  fout <- apply(X, 1, FUN = function(x) {
    ## take a particular row of X
    x <- as.matrix(x)
    
    ## the index for the "X" part
    xtb <- t(x) %*% t(betmat)
    
    ## figure out if we are in an "inner case" (i.e. standard case)
    ## or in one of the tails
    ul.pos <- tail(which(xtb-y <= 0),1) ## position of ul (in text)
    uu.pos <- head(which(xtb-y > 0),1) ## position of uu (in text)
    
    ## code to handle tails uniquely
    lam1 <- 1 - min(tau)
    lam2 <- max(tau)
    if (length(uu.pos) > 0) { ## add extra check for some missing cases; don't
      ## need for other side because of the ordering
      if (uu.pos == 1) { ##we are way on left tail
        return(min(tau) * lam1 * exp(lam1*(y - t(x) %*%betmat[1,])))
      }
    }
    if (ul.pos == length(tau)) { ## we are way on the right tail
      return( (1-max(tau)) * lam2 * exp(-lam2 * (y - t(x)%*%betmat[length(tau), ])))  
    }
    ## standard case ("inner case")
    (tau[uu.pos] - tau[ul.pos]) / t(x)%*%(betmat[uu.pos,] - betmat[ul.pos,])
  })
  fout
}

fv.yx <- function(v, betmat, k, pi, mu, sig, Y, X, tau) {
  fy.xvals <- sapply(1:length(Y), function(i) {
    fy.x(y = (Y[i] - v), X=t(X[i,]), betmat=betmat, tau=tau)
  })
  fy.xvals * fv(v,k,pi,mu,sig)
}

fv <- function(v,ksig=1,pi=1,mu=0,sig=1) {
  ## when just normal distribution
  ##dnorm(v, sd=sig)
  
  ## mixture of normals
  sum(sapply(1:ksig, function(i) {
    pi[i]/sig[i] * dnorm( (v-mu[i])/ sig[i] )
  }))
}


## metropolis-hastings algorithm
## this may not be "fast" at the moment
mh_mcmc <- function(startval=0, iters=500, burnin=100, betmat, k, pi, mu, sig, y, x, tau) {
  x <- t(x)
  out <- rep(NA, iters)
  out[1] <- startval
  for (i in 2:iters) {
    trialval <- out[i-1] + rnorm(1, sd=sqrt(4))
    fvold <- fv.yx(out[i-1], betmat, k, pi, mu, sig, y, x, tau)
    fvnew <- fv.yx(trialval, betmat, k, pi, mu, sig, y, x, tau)
    if ( fvnew > fvold ) {
      out[i] <- trialval
    } else {
      if ( (fvnew / fvold) >= runif(1)) {
        out[i] <- trialval
      } else {
        out[i] <- out[i-1]
      }
    }
  }
  return(tail(out, iters-burnin))
}




#' @title getParams
#'
#' @description Helper function to take the result of the optimization routine
#' and converts back to parameters
#'
#' @param optout results of call to optim
#' @param ksig the dimension of the mixture model for the measurement error
#'
#' @return list of parameters
#'
#' @export
getParams <- function(optout, formla, data, tau, nmix) {
    xformla <- formla
    xformla[[2]] <- NULL ## drop y variable
    x <- model.matrix(xformla, data)
    kx <- ncol(x)*length(tau)
    if (class(optout) == "numeric") {
        par <- optout
    } else { ## output of optim
        par <- optout$par
    }
    bet <- par[1:kx]
    k <- kx/length(tau)
    n <- nrow(x)
    ktau <- length(tau)
    bet <- split(bet,ceiling(seq_along(bet)/k))
    kmu <- nmix-1
    if (nmix > 1) {
        pi1 <- par[(kx+1):(kx+kmu)]
        mu1 <- par[(kx+kmu+1):(kx+kmu+kmu)]
        pi <- c(pi1, 1-sum(pi1))
        mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    } else {
        pi <- 1
        mu <- 0
    }
    ksig <- nmix
    sig <- par[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]
    out <- list(bet=bet, pi=pi, mu=mu, sig=sig)
    class(out) <- "PARAM"
    out
}

#' @title makeRQS
#'
#' @description Take the results of the optimization and convert them
#'  into a quantile regression object, so we can use all the tools from
#'  the quantreg package (e.g. inverting the quantiles).  The key step
#'  here is rearrangement, because the objective function doesn't impose
#'  any ordering -- see the discussion in HLLP.  We follow HLLP's
#'  recommendation and order the QR parameters by what makes the quantiles
#'  to be increasing for the mean values of the x's.  This means that
#'  for any particular value of the x's, the quantile function may not
#'  necessarily be increasing in tau.  However, we can separately rearrange
#'  those as needed.  But this gives a way to report the QR parameters.
#'
#' @param optout results of call to optim
#'
#' @return rqs object
#'
#' @export
makeRQS <- function(params, formla, data) {
    xformla <- formla
    xformla[[2]] <- NULL ## drop y variable
    x <- model.matrix(xformla, data)
    optout <- list()
    class(optout) <- c("rqs")
    optout$terms <- terms(formla)
    optout$formula <- formla
    bet <- params$bet
    ##bet <- split(bet1,ceiling(seq_along(bet1)/k))
    ## rearrangement (as in HLLP though double check)
    barx <- apply(x, 2, mean)
    betorder <- order(sapply(bet, function(b) sum(b*barx)))
    bet <- simplify2array(bet)
    if (class(bet)=="numeric") { ## i.e. there is just one element in bet
        bet <- as.matrix(bet)
    } else {
        bet <- t(bet)
    }
    bet <- as.matrix(bet[betorder,])
    optout$coefficients <- t(bet)
    optout$tau <- tau
    optout
}

#' @title parms2coppar
#'
#' @description convert parameters from optimization routine into copula
#'  parameters, parameters are between 0 and 1
#'
#' @inheritParams ll
#' @param params a vector of parameters
#'
#' @return a vector of copula parameters
#'
#' @export
parms2coppar <- function(params, copula="gumbel",x) {
    x <- as.matrix(x)
    ## if (ncol(x)==1) {
    ##     x <- t(x)
    ## }
    params <- log(params/(1-params))
    if (copula=="gumbel") {
        deltx <- 1 + exp(x%*%params)
    } else if (copula=="frank") {
        deltx <- x%*%params
    } else if (copula=="clayton") {
        deltx <- x%*%params
    } else if (copula=="gaussian") {
        deltx <- 2 * exp(x%*%params) / (1 + exp(x%*%params)) - 1
    } else {
        stop( paste0("copula ", copula, " not supported") )
    }
    deltx
}

getCBounds <- function(copula) {
    if (copula=="gumbel") {
        lo <- 1
        hi <- 1000 ##Inf
    } else if (copula=="frank") {
        lo <- -1000 ##-Inf
        hi <- 1000 ##Inf
    } else if (copula=="clayton") {
        lo <- -1000 ##-Inf
        hi <- 1000 ##Inf
    } else if (copula=="gaussian") {
        lo <- -1
        hi <- 1
    } else {
        stop( paste0("copula ", copula, " not supported") )
    }
    return(c(lo,hi))
}


qr2me <- function(yname, tname, xformla, tau, data, xdf=NULL, tvals=NULL,
                  copula="gumbel", Qyx, Qtx, 
                  startparams=NULL, method="BFGS", maxit=1000, ndraws=100,
                  cl=1) {

    cat("\nqr2me method...\n")
    cat("----------------------")
    cat("\nCitation: Callaway, Brantly, Tong Li, and Irina Murtazashvili, Quantile Treatment Effects with Two-Sided Measurement Error, Working Paper, 2018....\n")
    cat("----------------------")
    cat("\n")
    x <- model.matrix(xformla, data)
    if (is.null(startparams)) {
        ##startparams <- rep(0, ncol(x)) ## do this if want params to be function of X
        startparams <- 0.5
    }

    Fyx <- predict(Qyx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)
    Ftx <- predict(Qtx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)

    Fyx <- lapply(Fyx, rearrange)
    Ftx <- lapply(Ftx, rearrange)

    fyx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")
    ftx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")
    
    

    
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


    ##ll(startparams, data[,yname], data[,tname], x=x, copula=copula, Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx, Us=Us, Vs=Vs)
    
    cat("\nStep 1 of 2: Estimating copula parameter...\n")

    res <- optimize(ll, c(0,1), maximum=TRUE, 
                 y=data[,yname], t=data[,tname], x=x, copula=copula,
                 Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx,
                 Us=Us, Vs=Vs)

    delt <- rep(parms2coppar(res$maximum, copula=copula, x=1), nrow(x))

    
    if (!is.null(tvals)) {

        if (is.null(xdf)) xdf <- as.data.frame(t(apply(x,2,mean))) ##x, for all data
 
        ##xtdf <- cbind(tvals, xdf)
        ## todo, this gives the copula, now convert to conditional distribution
         ##delt <- parms2coppar(res$par, copula, xdf)
        ##tvals <- quantile(data[,tname], tau, type=1)
        yvals <- quantile(data[,yname], seq(.01,.99,.01)) ## could also take all unique yvals or let user pass them all in
        yvals <- yvals[order(yvals)]
        QQyx <- predict( Qyx, newdata=as.data.frame(rbind(xdf,x[1,])), stepfun=TRUE)  ## super hack:  but predict.rqs is throwing an error that I think it shouldn't, and this gets around it.
        QQyx <- QQyx[-length(QQyx)]
        QQyx <- lapply(QQyx, rearrange)
        Ftx <- predict(Qtx, newdata=as.data.frame(rbind(xdf, x[1,])),
                        type="Fhat", stepfun=TRUE)
        Ftx <- Ftx[-length(Ftx)]
        Ftx <- lapply(Ftx, rearrange)

        U <- seq(0,1,length.out=ndraws)
        if (copula=="frank") {
            cop <- copula::frankCopula(as.numeric(delt[1])) ## all delts restricted to be the same so just choose first one
        } else if (copula=="gumbel") {
            cop <- copula::gumbelCopula(as.numeric(delt[1]))
        } else if (copula=="clayton") {
            cop <- copula::claytonCopula(as.numeric(delt[1]))
        } else if (copula=="gaussian") {
            cop <- copula::normalCopula(as.numeric(delt[1]))
        } else {
            stop(paste0("copula type:", copula, " is not supported"))
        }

        cat("\nStep 2 of 2: Building conditional distributions...\n")
        Fytx <- pblapply(1:length(QQyx), function(i) {
            qfun <- QQyx[[i]]
            ffun <- Ftx[[i]]
            lapply(tvals, function(tt) {
                BMisc::makeDist(yvals, sapply(yvals, function(yy) {
                    mean(1*(qfun(U)<=yy)*dCopula(cbind(U, ffun(tt)), cop))
                }))
            }) ## might want to make this a function with makeRQS (modified) or makeDist eventually
        }, cl=cl)
        ##Qytx <- lapply(Fytx, function(FF) quantile(FF, tau, type=1))

        ## next we want to reverse the list
        reverseListIndex <- function(l) {
            outlist <- list()
            length(outlist) <- length(l[[1]])
            outlist <- lapply(outlist, function(f) {
                g <- list()
                length(g) <- length(l)
                g
            })
        
            for (i in 1:length(l)) {
                for (j in 1:length(l[[1]])) {
                    outlist[[j]][[i]] <- l[[i]][[j]]
                }
            }
            outlist
        }

        Fytx <- reverseListIndex(Fytx)
            

        ##Fyt <- lapply(Fytx, function(FFytx) {
        ##    combineDfs(yvals, FFytx)
        ##})

        out <- list(cop.param=parms2coppar(res$maximum, copula=copula, x=1),
                    copula=copula, Fytxlist=Fytx, tvals=tvals, x=xdf)
        
        ##QytxOut <- list(Qytx=Qytx, x=xdf, t=tvals, tau=tau, delt=delt)

    ### only do above if you want the results for a particular value of t and x;
    ### otherwise can just return all results 
    } else {
        out <- list(cop.param=parms2coppar(res$maximum, copula=copula, x=x),
                    copula=copula)
    }

    out$Qyx <- Qyx
    out$Qtx <- Qtx

    class(out) <- "qr2meobj"

    return(out)
}

avgDist <- function(Fytxlist, yvals) {
    
    Fyt <- lapply(Fytx, function(FFytx) {
        combineDfs(yvals, FFytx)
    })

    Fyt
}
    

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

print.qr2meobj <- function(obj) {
    print(obj$Qyx)
    print(obj$Qtx)
    cat("\n\n")
    cat("Copula Type: ")
    cat(obj$copula)
    cat("\n")
    cat("Copula Paramter: ")
    cat(obj$cop.param)
    cat("\n\n")
}

print.merr <- function(obj) {
    coef <- round(t(as.matrix(obj$coefficients)),4)
    rownames(coef) <- obj$tau
    colnames(coef) <- c("Intercept", attr(obj$terms,"term.labels"))
    cat("Coefficients:\n")
    print(coef)
    cat("\n\n")
    cat("Measurement Error Distribution:\n")
    U <- round(cbind(obj$pi, obj$mu, obj$sig^2),4)
    rownames(U) <- sapply(1:nrow(U), function(i) paste0("Comp. ",i))
    colnames(U) <- c("Prob.", "Mean", "Variance")
    print(U)    
}

getListElement <- function(listolists, whichone=1) {
    lapply(listolists, function(l) l[[whichone]])
}

plot.qr2meobj <- function(obj, whichone=1, tau=c(.1,.5,.9), ylim=NULL,
                          ylab=NULL, xlab=NULL) {
    qq <- t(sapply(getListElement(obj$Fytxlist, whichone),
                   function(FF) quantile(FF, probs=tau)))
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

addplot <- function(obj, p, whichone=1, tau=c(.1,.5,.9)) {
    qq <- t(sapply(getListElement(obj$Fytxlist, 1),
                   function(FF) quantile(FF, probs=tau)))
    cmat <- cbind.data.frame(tvals, qq)
    colnames(cmat) <- c("tvals", paste0("c",tau*100))
    cmat <- gather(cmat, quantile, value, -tvals)
    p <- p + geom_line(data=cmat, mapping=aes(x=tvals, y=value, group=quantile,
                                              color=quantile),
                       linetype="dashed")
    p <- p + geom_point(data=cmat, mapping=aes(x=tvals, y=value, group=quantile,
                                               color=quantile))
    p
}


qrme.alt <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL, startsig=NULL, startpi=NULL, method="BFGS", maxit=1000, maxemiters=25, tol=1) {

        
    newparams <- em.alt(formla, tau, data, nmix, startbet, startmu, startsig, startpi, method, maxit)

    

    cval <- tol+1
    i <- 1

    while( (cval > tol) & (i <= maxemiters)) {
        paramvals <- newparams
        params <- getParams(newparams, formla, data, tau, nmix)
        sbet <- unlist(params$bet)
        ssig <- params$sig
        ksig <- length(ssig)
        smu <- params$mu[-ksig]
        spi <- params$pi[-ksig]
        newparams <- em.alt(formla, tau, data, nmix, sbet, smu, ssig, spi, method, maxit)
        cval <- sum((paramvals-newparams)^2)
        cat("iter: ", i, "\n")
        cat("change in objective function: ", round(cval,4),"\n")
        i <- i + 1
        ## if (cval < tol) {
        ##     break
        ## }
    }
    
    out <- makeRQS(getParams(newparams, formla, data, tau, nmix), formla, data)

    class(out) <- c("merr", class(out))
    
    out$params <- newparams
        
    ## idx <- length(betvals)+1
    ## nidx <- length(betvals)+length(pivals)
    ## pi1 <- res$par[idx:nidx]
    ## out$pi <- c(pi1, 1-sum(pi1))
    ## idx <- nidx+1
    ## nidx <- nidx+length(muvals)
    ## mu1 <- res$par[idx:nidx]
    ## out$mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    ## idx <- nidx+1
    ## nidx <- nidx+length(sigvals)
    ## out$sig <- res$par[idx:nidx]
    out$pi <- params$pi
    out$mu <- params$mu
    out$sig <- params$sig

    out
}

em.alt <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL, startsig=NULL, startpi=NULL, method="BFGS", maxit=1000) {
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
        ## betvals <- unlist(lapply(1:length(tau), function(i) c(quantile(y, probs=tau[i], type=1), rep(0, (k-1)))))
        betvals <- c(rq(formla, tau=tau, data=data)$coefficients)
        
    } else {
        betvals <- startbet
    }
    if (is.null(startmu)) {
        muvals <- seq(-1, by=1, length.out=(m-1))
    } else {
        muvals <- startmu
    }
    if (is.null(startsig)) {
        sigvals <- rep(1,m)
    } else {
        sigvals <- startsig
    }
    if (is.null(startpi)) {
        pivals <- rep(1/m, (m-1))
    } else {
        pivals <- startpi
    }
    qrparams <- c(betvals, pivals, muvals, sigvals)
    
    x <- as.matrix(x)
    kx <- ncol(x)*length(tau)
    ksig <- nmix
    kmu <- nmix-1
    params <- qrparams
    bet <- params[1:kx]
    k <- kx/length(tau)
    n <- nrow(x)
    y <- data[,yname]
    ktau <- length(tau)
    bet <- split(bet,ceiling(seq_along(bet)/k))
    if (kmu > 0) {
        pi1 <- params[(kx+1):(kx+kmu)]
        mu1 <- params[(kx+kmu+1):(kx+kmu+kmu)]
        pi <- c(pi1, 1-sum(pi1))
        mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    } else {
        pi <- 1
        mu <- 0
    }
    
    sig <- params[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]

    ndraws <- 100
    ##Ucomponents <- sample(1:ksig, ndraws, replace=TRUE, prob=Upi)
    ##Us <- rnorm(ndraws, Umu[Ucomponents], Usig[Ucomponents])

    out <- list()
    out$par <- params
    Qyx.rqs <- makeRQS(getParams(out, formla, data, tau, 1), formla, data)

    Qyx <- predict(Qyx.rqs, newdata=data, type="Qhat", stepfun=TRUE)
    Fyx <- predict(Qyx.rqs, newdata=data, type="Fhat", stepfun=TRUE)
    fyx <- predict(Qyx.rqs, newdata=data, type="fhat")

    Qyx <- rearrange(Qyx)
    Fyx <- rearrange(Fyx)

    Qyx1 <- lapply(Qyx, function(q) {
        approxfun(c(0,tau,1), c(min(y), q(tau), max(y)))
    })

    Fyx1 <- lapply(Qyx1, function(q) {
        y <- seq(min(y), max(y), length.out=101)
        tt <- seq(0,1,.01)
        approxfun(y, sapply(y, function(yy) sum(1*(q(tt) <= yy))/length(tt)),
                  yleft=0, yright=1)
    })
    
    qyxi <- function(u, i) {
        taul <- tau[max(which(tau<=u))]
        tauu <- tau[min(which(tau>u))]
        ql <- Qyx1[[i]](taul)
        qu <- Qyx1[[i]](tauu)
        if (u < min(tau)) {
            ql <- min(y)
            taul <- 0
        }
        if (u >= max(tau)) {
            qu <- max(y)
            tauu <- 1
        }
        max((qu-ql)/(tauu-taul), 1e-20)
        ##Qyx[[i]](taul) + (u-taul)*(Qyx[[i]](tauu) - Qyx[[i]](taul))/(tauu-taul) ## this just smooths
    }

    fyexi <- function(y, e, i) {
        1 / qyxi(Fyx1[[i]](y - e), i)
        ##ff <- fyx[[i]](y-e)
        ##if (is.na(ff)) ff <- 1e-10
        ##ff
    }
    
    u <- lapply(bet, function(b) {
        b <- as.matrix(b)
        y - x%*%b
    })

    fe <- function(e) {
        sum( pi * dnorm( ( e - mu ) / sig) )
    }
    
    fu <- sapply(u, function(uu) {
        apply(sapply(1:ksig, function(i) {
            pi[i]/sig[i] * dnorm( (uu - mu[i]) / sig[i] )
        }), 1, sum)
    }) ## this will contain a matrix with n rows and ktau columns

    fu <- apply(fu, 1, sum)
    
    
    ## this is the posterior distribution
    feyxi <- function(e, i) {
        fyexi(y[i], e, i) * fe(e) / fu[i]
    }

    ##metropolis hastings
    e0 <- 0
    evec <- c()
    
    elist <- list()
    set.seed(43952)
    rns <- rnorm(600*n, sd=sqrt(var(y)/2))
    unis <- runif(600*n)
    counter <- 1
    for (i in 1:n) {
        for (j in 1:600) {
            e1 <- e0 + rns[counter]
            fe1 <- feyxi(e1, i)
            fe0 <- feyxi(e0, i)
            if (fe1 > fe0) {
                evec <- c(evec, e1)
                e0 <- e1
            } else {
                u <- unis[counter]
                if (u <= fe1/fe0) {
                    evec <- c(evec, e1)
                    e0 <- e1
                } else {
                    evec <- c(evec, e0)
                }
            }
            counter <- counter+1
        }
        evec <- evec[501:600] ## can update these later, increase to around 1000 burn-in
        elist[[i]] <- evec
    }


    ## accept-reject simulator
    ## elist <- list()
    ## set.seed(43951)
    ## rns <- rnorm(5000*n, sd=sqrt(2))
    ## unis <- runif(5000*n)
    ## counter <- 1
    ## for (i in 1:n) {
    ##     evec <- c()
    ##     for (j in 1:5000) {
    ##         if (feyxi(rns[counter],i) > unis[counter]) {
    ##             evec <- c(evec, rns[counter])
    ##         }
    ##         counter <- counter + 1
    ##     }
    ##     elist[[i]] <- evec
    ## }

    ## elist <- lapply(elist, function(ee) ee[1:min(sapply(elist, length))])
    



    ddf <- lapply(1:n, function(i) {
        cbind(matrix(rep(c(y[i], x[i,]), length(elist[[i]])), nrow=length(elist[[i]]), byrow=TRUE), elist[[i]])
    })

    ddf <- do.call(rbind.data.frame, ddf)
    colnames(ddf) <- c(yname, colnames(x), "E")

    ddf$Ystar <- ddf[,yname] - ddf[,"E"]

    formla[[2]] <- as.name("Ystar")

    tau1 <- seq(.1, 0.9, .1)
    betp <- t(coef(rq(formla, tau=tau1, data=ddf)))    
    

    cbind(tau1, b0Y(tau1), b1Y(tau1), betp)

    


    llme.alt <- function(params, y, x, bet, tau, ksig, kmu=(ksig-1)) {
        x <- as.matrix(x)
        kx <- 0 ## this a bit hack cause just copied and pasted, but will work ##ncol(x)*length(tau)
        k <- kx/length(tau)
        n <- nrow(x)
        ktau <- length(tau)
        if (kmu > 0) {
            pi1 <- params[(kx+1):(kx+kmu)]
            mu1 <- params[(kx+kmu+1):(kx+kmu+kmu)]
            pi <- c(pi1, 1-sum(pi1))
            mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
        } else {
            pi <- 1
            mu <- 0
        }

        sig <- params[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]

        
        u <- lapply(bet, function(b) {
            b <- as.matrix(b)
            y - x%*%b
        })
        
        fu <- sapply(u, function(uu) {
            apply(sapply(1:ksig, function(i) {
                pi[i]/sig[i] * dnorm( (uu - mu[i]) / sig[i] )
            }), 1, sum)
        }) ## this will contain a matrix with n rows and ktau columns

        ll <- log(apply(fu, 1, sum)) ## p. 29 in hausman, liu, luo, palmer

        return(-sum(ll))

    }


    llme.gr.alt <- function(params, y, x, bet, tau, ksig, kmu=(ksig-1)) {
        x <- as.matrix(x)
        kx <- 0 ##ncol(x)*length(tau)
        k <- kx/length(tau)
        n <- nrow(x)
        ktau <- length(tau)
        if (kmu > 0) {
            pi1 <- params[(kx+1):(kx+kmu)]
            mu1 <- params[(kx+kmu+1):(kx+kmu+kmu)]
            pi <- c(pi1, 1-sum(pi1))
            mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
        } else {
            pi <- 1
            mu <- 0
        }

        sig <- params[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]
        ## u is a list with as many elements as the length of tau
        ## each element contains an nx1 vector containing y - x'beta(tau)
        u <- lapply(bet, function(b) {
            b <- as.matrix(b)
            y - x%*%b
        }) 
        ## 

        ## fu contains the density of u (the measurement error) evaluated
        ## at (y - x'beta(tau)) / sig, and given the values of the parameters
        ## for the mixture model
        ## this will contain a matrix with n rows and ktau columns
        fu <- sapply(u, function(uu) {
            apply(sapply(1:ksig, function(i) {
                pi[i]/sig[i] * dnorm( (uu - mu[i]) / sig[i] )
            }), 1, sum)
        }) 

        ## ifu contains the integrated values (over 0 to 1) of fu
        ## this object is useful throughout and is an nx1 vector
        ifu <- apply(fu, 1, sum) ## this is integration step

        ## the gradient with respect to pi
        ## this should contain a vector of length ksig x 1
        gr.pi <- function() {
            p1 <- 1/ifu

            p2a <- simplify2array(lapply(u, function(uu) {
                sapply(1:ksig, function(i) {
                    1/sig[i] * dnorm( (uu - mu[i]) / sig[i] )
                })
            }))
            
            p2 <- apply(p2a, c(1,2), sum) ## integration step, should be nxksig

            as.numeric(apply(p2*p1, 2, sum))
        }

        gr.mu <- function() {
            p1 <- 1/ifu

            p2a <- simplify2array(lapply(u, function(uu) {
                sapply(1:ksig, function(i) {
                    pi[i]/sig[i] * (uu - mu[i] / sig[i]^2) * dnorm( (uu - mu[i]) / sig[i] )
                })
            }))
            p2 <- apply(p2a, c(1,2), sum) ## integration step, should be nxksig

            as.numeric(apply(p2*p1, 2, sum))
        }

        gr.sig <- function() {
            p1 <- 1/ifu

            p2a <- simplify2array(lapply(u, function(uu) {
                sapply(1:ksig, function(i) {
                    -pi[i]/sig[i]^2 * dnorm( (uu - mu[i]) / sig[i] )
                })
            }))

            p2b <- simplify2array(lapply(u, function(uu) {
                sapply(1:ksig, function(i) {
                    1/sig[i] * (uu - mu[i])^2 / sig[i]^3 * dnorm( (uu - mu[i]) / sig[i] )
                })
            }))

            p2 <- apply(p2a+p2b, c(1,2), sum) ## integration step, should be nxksig
            
            as.numeric(apply(p2*p1, 2, sum))
        }

        cgr.pi <- function() {
            (gr.pi() - gr.pi()[ksig] - (mu*(1-sum(pi[-ksig]))+sum( (mu*pi)[-ksig]))*gr.mu()[ksig]/(1-sum(pi[-ksig]))^2)[-ksig]
        }

        cgr.mu <- function() {
            (gr.mu() - pi*gr.mu()[ksig]/(sum(pi[-ksig])))[-ksig]
        }
        
        c(-cgr.pi(), -cgr.mu(), -gr.sig())
    }
    
    eparams <- c(pivals, muvals, sigvals)

    cst <- NULL
    ust <- NULL
    ## check if this part is right, would notice if violating constraints
    if (nmix > 1) {
        cst <- c(rep(.01, length(pivals)),
                 -.99, ## thisone is so that 1-pi1-pi2>= .01
                 rep(.01, length(sigvals)))

        u1 <- rep(0, length(eparams))

        ust <- t(simplify2array(c(lapply(1:length(pivals), function(i) {
            u2 <- u1
            u2[i] <- 1
            u2
        }), lapply(length(pivals), function(l) {
            u2 <- u1
            for (i in 1:l) {
                u2[i] <- -1
            }
            u2
        }), lapply(1:length(sigvals), function(i) {
            u2 <- u1
            u2[length(pivals) + length(muvals) + i] <- 1
            u2
        }))))
    }


    cat("\n\n currently testing estimating bets, dont forget to uncomment code \n\n")
    
    ## if (nmix == 1) {
    ##     res <- optim(eparams, llme.alt, gr=llme.gr.alt, method="BFGS",
    ##                  y=y, x=x, tau=tau, ksig=m, bet=bet,
    ##                  control=list(maxit=maxit, trace=6))
    ## } else {
    ##     res <- constrOptim(eparams, llme.alt, llme.gr.alt, ui=ust, ci=cst,
    ##                        y=y, x=x, tau=tau, bet=bet, ksig=m,
    ##                        method="BFGS",
    ##                        control=list(maxit=maxit, trace=6))
    ## }


    ## just for testing in the case where distribution is known
    res$par <- eparams
    ###

    out <- c(t(betp), res$par)


    
    out
}
