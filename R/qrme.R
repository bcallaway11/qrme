#' @title cop.pdf
#' @description returns the pdf of several copula functions evaluated at
#'  u and v; the available copula functions are:  Gumbel (the default),
#'  Frank, (should implement Joe and Gaussian)
#'
#' @param u first copula argument (u and v can be vectors)
#' @param v second copula argument
#' @param type one of c("gumbel","frank")
#' @param delt the copula paramter
#' @param eps A small positive number to trim out zeros
#'
#' @return vector of copula pdf values
#'
#' @export
cop.pdf <- function(u,v,type="gumbel",delt,eps=1e-300) {
    if ( !(all(0 <= u) & all(u <= 1)) | !(all(0 <= v) & all(v <= 1))) {
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
#' @param Us draws for the measurement error in the outcome equation (these
#'  are made in the function that calls the likelihood function)
#' @param Vs draws for the measurement error in the treatment equation (these
#'  are made in the function that calls the likelihood function)
#' 
#' @return scalar negative value of log likelihood function 
#'
#' @keywords internal
#' 
#' @export
ll <- function(params, y, t, x, copula="gumbel", Fyx, Ftx, fyx, ftx,
               Us, Vs, ndraws=100, eps=1e-300) {
    
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
#' @param tol This is the convergence criteria.  When the change in the
#'  Euclidean distance between the new parameters (at each iteration) and
#'  the old parameters (from the previous iteration) is smaller than tol,
#'  the algorithm concludes.  In general, larger values for tol will result
#'  in a fewer number of iterations and smaller values will result in more
#'  accurate estimates.
#' @param cl The numbe of clusters to use for parallel computation (default
#'  is 1 so that computation is not done in parallel)
#' 
#' @return an object of class "merr"
#'
#' @export
qrme <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL,
                 startsig=NULL, startpi=NULL, tol=1, iters=400,
                 burnin=200, drawsd=4, cl=1) {
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
                   sigguess=sigvals, tol=tol,
                   iters=iters, burnin=burnin, drawsd=drawsd, cl=cl)


    out <- makeRQS(res, formla, data, tau=tau)

    class(out) <- c("merr", class(out))

    ## set final parameters for measurement error class
    out$bet <- res$bet
    out$pi <- res$pi
    out$mu <- res$mu
    out$sig <- res$sig
    out$Ystar <- res$Ystar

    out
}


#' @title em.algo
#'
#' @description A pseudo EM algorithm for quantile regression with measurement error.  The measurement error here follows a mixture of normals.
#' @param betmatguess Initial values for the beta parameters.  This should be
#'  an LxK matrix where L is the number of quantiles and K is the dimension
#'  of the covariates
#' @param tau An L-vector indicating which quantiles have been estimated
#' @param m The number of components of the mixture distribution for the
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
em.algo <- function(formla, data,
                    betmatguess, tau, m=1, piguess=1, muguess=0,
                    sigguess=1, tol=.01,
                    iters=400, burnin=200, drawsd=4, cl=1) {
    
    stopIters <- 100
    counter <- 1
    
    while (counter <= stopIters) {
        newone <- em.algo.inner(formla, data,
                                betmatguess, tau, m, piguess, muguess, sigguess, iters=iters, burnin=burnin, drawsd=drawsd, cl=cl)
        newbet <- newone$bet
        newpi <- newone$pi
        newmu <- newone$mu
        newsig <- newone$sig
        
        cat("\n\n\nIteration: ", counter, "\n pi: ", newpi, "\n mu: ", newmu, "\n sig: ", newsig, "\n\n")
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

#' @title Inner part of EM-algorithm for QR with measurement error
#'
#' @description Does the heavy-lifting of the EM-algorithm for QR with
#'  measurment error
#' @inheritParams em.algo
#' @inheritParams fv.yx
#'
#' @keywords internal
#'
#' @return A list of QR parameters and parameters for mixture of normals for
#'  the measurement error term
#' @export
em.algo.inner <- function(formla, data, 
                          betmat, tau, m=1, pi=1, mu=0, sig=1,
                          iters=400, burnin=200, drawsd=4, cl=1) {


    xformla <- formla
    xformla[[2]] <- NULL ## drop y variable
    X <- model.matrix(xformla, data)
    yname <- as.character(formla[[2]])
    Y <- data[,yname]
    k <- ncol(X) ## number of x variables
    n <- length(Y)
    cat("\nSimulating measurement error...")

    
    ## newdta <- pbapply::pblapply(1:n, function(i) {
    ##     ## old call to R code
    ##     ##e <- mh_mcmc(betmat=betmat, m=m, pi=pi, mu=mu, sig=sig, y=Y[i], x=X[i,], tau=tau, iters=1000, burnin=500)
    ##     ## new call to C++ code
    ##     e <- mh_mcmc_innerC(startval=0, iters=400, burnin=200, drawsd=4, betmat=betmat, m=m, pi=pi, mu=mu, sig=sig, y=Y[i], x=as.matrix(X[i,]), tau=tau)
    ##     cbind.data.frame(Y[i]-e, matrix(X[i,], nrow=length(e), ncol=ncol(X), byrow=TRUE),e=e)
    ## }, cl=cl)

    startval <- 0
    edraws <- mh_mcmcC(Y, X, startval=startval, iters=iters, burnin=burnin,
                       drawsd=drawsd, betmat=betmat, m=m,
                       pi=pi, mu=mu, sig=sig, tau=tau)

    newids <- unlist(lapply(1:n, function(i) rep(i, (iters-burnin)))) ## just replicates Y and X over and over

    newdta1 <- as.data.frame(cbind(Y=(Y[newids]-edraws), X=X[newids,], e=edraws))

    ## Note: this currently just works for one X; will need to update

    cat("\nEstimating QR including simulated measurement error...")
    ##newdta1 <- do.call(rbind.data.frame, newdta)
    colnames(newdta1) <- c(yname, colnames(X), "e")
    out <- quantreg::rq(formla, tau=tau, data=newdta1, method="pfn")
    
    ## this is part I am not sure about, once you have a new beta then estimate a new sigma??
    ## also should probably restrict overall mean of error term to be equal to 0
    cat("\nEstimating finite mixture model...")
    if (m == 1) {
        nm <- list(m=1, lambda=1, mu=0, sigma=sd(newdta1$e))
    } else {
        nm <- mixtools::normalmixEM(newdta1$e, k=m, epsilon=1e-03)
    }
    
    return(list(bet=t(coef(out)), m=m, pi=nm$lambda, mu=nm$mu, sig=nm$sigma, Ystar=newdta1[,yname]))
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
#' @title fy.x
#'
#' @description return the density of Y (for particular value y) conditional
#'  on X (which can include n observations) when Q(Y|X) has been estimated
#'  using QR.  This is used in our simulation approach.
#' @param y particular value of y to estimate f(y|x)
#' @param betmat LxK matrix of parameter values with L the number of quantiles
#'  and K the dimension of the covariates
#' @param XX An nxK matrix with n the number of observations of X
#' @param tau an L-vector containing the quantile at which Q(Y|X) was estimated
#'
#' @return An nx1 vector that contains f(y|X)
#'
#' @export
fy.x <- function(y, betmat, XX, tau) {
  X <- as.matrix(XX)
  
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


#' @title fv.yx
#'
#' @description Return the density of the measurement error conditional
#'  on y and x; this takes as given some QR parameters from Y* (the true
#'  outcome) conditional on X.  Here, we also presume that the distribution
#'  of the measurement error is a mixture of normal distributions
#' @param v A particular value of the measurement error to estimate f(v|y,x)
#' @inheritParams fy.x
#' @param m The dimension of the measurement error
#' @param pi The probability of each mixture component (should have length
#'  equal to m)
#' @param mu The mean of each mixture component (should have length equal
#'  to m)
#' @param sig The standard deviation of each mixture component (should have
#'  length equal to m)
#' @param Y An nx1 vector of outcomes
#' @param X An nxK matrix of covariates
#' @param tau an L-vector of all the quantiles where betas were estimated
#' @return n x 1 vector of f(v|Y,X)
#'
#' @keywords internal
#' @export
fv.yx <- function(v, betmat, m, pi, mu, sig, Y, X, tau) {
  fy.xvals <- sapply(1:length(Y), function(i) {
    fy.x(y = (Y[i] - v), X=t(X[i,]), betmat=betmat, tau=tau)
  })
  fy.xvals * fv(v,m,pi,mu,sig)
}


#' @title fv
#'
#' @description The distribution of the measurement error using a mixture of
#'  normal distributions
#'
#' @inheritParams fv.yx
#' @return scalar f(v)
#' @keywords internal
#' @export
fv <- function(v,m=1,pi=1,mu=0,sig=1) {
  ## mixture of normals
  sum(sapply(1:m, function(i) {
    pi[i]/sig[i] * dnorm( (v-mu[i])/ sig[i] )
  }))
}


#' @title mh_mcmc
#' @description A Metropolis-Hastings algorithm for drawing measurment errors.
#'
#' @param startval The first value in the markov chain
#' @param iters The total number of measurement error draws to make
#' @param burnin The number of draws to drop
#' @param drawsd Trial values are drawn from N(0, sd=drawsd), default is 4
#' @inheritParams fv.yx
#' @param y particular value of y
#' @param x particular value of x
#' @return vector of draws of measurement error
#' @export
mh_mcmc <- function(startval=0, iters=500, burnin=100, drawsd=sqrt(4), betmat, m, pi, mu, sig, y, x, tau) {
    x <- t(x)
    out <- rep(NA, iters)
    out[1] <- startval
    for (i in 2:iters) {
        trialval <- out[i-1] + rnorm(1, sd=drawsd)
        fvold <- fv.yx(out[i-1], betmat, m, pi, mu, sig, y, x, tau)
        fvnew <- fv.yx(trialval, betmat, m, pi, mu, sig, y, x, tau)
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



## seems to be unused since we switched to EM algorithm
## #' @title getParams
## #'
## #' @description Helper function to take the result of the optimization routine
## #' and converts back to parameters
## #'
## #' @param optout results of call to optim
## #' @param ksig the dimension of the mixture model for the measurement error
## #'
## #' @return list of parameters
## #'
## #' @export
## getParams <- function(optout, formla, data, tau, nmix) {
##     xformla <- formla
##     xformla[[2]] <- NULL ## drop y variable
##     x <- model.matrix(xformla, data)
##     kx <- ncol(x)*length(tau)
##     if (class(optout) == "numeric") {
##         par <- optout
##     } else { ## output of optim
##         par <- optout$par
##     }
##     bet <- par[1:kx]
##     k <- kx/length(tau)
##     n <- nrow(x)
##     ktau <- length(tau)
##     bet <- split(bet,ceiling(seq_along(bet)/k))
##     kmu <- nmix-1
##     if (nmix > 1) {
##         pi1 <- par[(kx+1):(kx+kmu)]
##         mu1 <- par[(kx+kmu+1):(kx+kmu+kmu)]
##         pi <- c(pi1, 1-sum(pi1))
##         mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
##     } else {
##         pi <- 1
##         mu <- 0
##     }
##     ksig <- nmix
##     sig <- par[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]
##     out <- list(bet=bet, pi=pi, mu=mu, sig=sig)
##     class(out) <- "PARAM"
##     out
## }

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
#' @param params an LxK matrix of QR parameters where L is the number of
#'  quantiles that parameters have been estimated at and K is the dimension
#'  of the covariates.
#' @param formla y ~ x, a formula for the outcome on the regressors
#' @param data a data.frame containing y and x
#' @param tau the vector of quantiles where QR was estimated
#'
#' @return rqs object
#'
#' @export
makeRQS <- function(params, formla, data, tau) {
    xformla <- formla
    xformla[[2]] <- NULL ## drop y variable
    x <- model.matrix(xformla, data)
    optout <- list()
    class(optout) <- c("rqs")
    optout$terms <- terms(formla)
    optout$formula <- formla
    bet <- params$bet
    #################
    ## OLD: this is for if you use ML as in HLLP
    ##bet <- split(bet1,ceiling(seq_along(bet1)/k))
    ## rearrangement (as in HLLP though double check)
    ## barx <- apply(x, 2, mean)
    ## betorder <- order(sapply(bet, function(b) sum(b*barx)))
    ## bet <- simplify2array(bet)
    ## if (class(bet)=="numeric") { ## i.e. there is just one element in bet
    ##     bet <- as.matrix(bet)
    ## } else {
    ##     bet <- t(bet)
    ## }
    ##bet <- as.matrix(bet[betorder,])
    ####################
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
qr2me <- function(yname, tname, xformla, tau, data, xdf=NULL, tvals=NULL,
                  copula="gumbel", Qyx, Qtx, 
                  startparams=NULL, method="BFGS", maxit=1000, ndraws=100) {

    cat("\nqr2me method...\n")
    cat("----------------------")
    cat("\nCitation: Callaway, Brantly, Tong Li, and Irina Murtazashvili, Quantile Treatment Effects with Two-Sided Measurement Error, Working Paper, 2018....\n")
    cat("----------------------")
    cat("\n")
    x <- model.matrix(xformla, data)
    n <-  nrow(data)
    if (is.null(startparams)) {
        ##startparams <- rep(0, ncol(x)) ## do this if want params to be function of X
        startparams <- 0.5
    }


    Fyx1 <- predict(Qyx, newdata=as.data.frame(x)) ## this gives nxL matrix
    Fyx <- lapply(1:nrow(Fyx1), function(i) BMisc::makeDist(Fyx1[i,], Fx=tau, sorted=TRUE))
    Ftx1 <- predict(Qtx, newdata=as.data.frame(x))
    Ftx <- lapply(1:nrow(Ftx1), function(i) BMisc::makeDist(Ftx1[i,], Fx=tau, sorted=TRUE))
    ## Take passed in quantiles and create their conditional distribution
    ## Fyx <- predict(Qyx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)
    ## Ftx <- predict(Qtx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)

    ## ## Rearrangement so they are actually a distribution function
    ## Fyx <- lapply(Fyx, quantreg::rearrange)
    ## Ftx <- lapply(Ftx, quantreg::rearrange)

    ## ## Also, get their density (might want to change how we do this)
    ## fyx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")
    ## ftx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")

    cat("Step 1 of 3: Converting QR to conditional density estimates...\n\n")
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
    cat("\nStep 2 of 3: Estimating copula parameter...\n")
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
    newdta1 <- data.frame(Y=Qyx$Ystar, T=Qtx$Ystar)
    ranks1 <- pobs(newdta1)
    cop <- copula::fitCopula(cop, ranks1, method="irho") ## irho inverts spearman's rho; it
    ##is very fast though (I think) not all copulas (exception=(I think)Gumbel) have 1-1 relationship
    ##with Spearman's rho, but in practice they seem very similar.


    browser()
    delt <- rep(attributes(cop)$estimate, nrow(x))

    ##estimation with maximum likelihood as in the original version of the paper
    res <- optimize(ll, c(0,1), maximum=TRUE, 
                 y=data[,yname], t=data[,tname], x=x, copula=copula,
                 Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx,
                 Us=Us, Vs=Vs)

    delt <- rep(parms2coppar(res$maximum, copula=copula, x=1), nrow(x))

    
    if (!is.null(tvals)) {

        cat("\nStep 3 of 3: Building conditional distributions...\n")
        ## If you don't set particular values of X to compute,
        ## just set it equal to the average values of X in the dataset
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
        QQyx2  <- predict(Qyx, newdata=as.data.frame(xdf))
        FFyx2  <- t(sapply(Fyx, function(fyx) fyx(yvals)))
        Ftx <- predict(Qtx, newdata=as.data.frame(rbind(xdf, x[1,])),
                        type="Fhat", stepfun=TRUE)
        Ftx <- Ftx[-length(Ftx)]
        Ftx <- lapply(Ftx, rearrange)
        FFtx <- t(sapply(Ftx, function(ftx) ftx(tvals)))

        U <- seq(0,1,length.out=ndraws)
        U <- tau

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

        ## call C++ function to return matrix with distribution of Fytx
        FytXmat <- computeFytXC(yvals, tvals, QQyx2, FFtx, tau, "gumbel", delt[1])

        
        ## internal function for reordering arguments of BMisc::makeDist
        makeDist1 <- function(Fx, x, sorted = FALSE, rearrange=FALSE) {
            BMisc::makeDist(x, Fx, sorted, rearrange)
        }

        ## converts 3-dimensional matrix of Fytx into 2-dimensional matrix of distributio
        ## functions
        Fytx  <- apply(FytXmat, c(1,3), makeDist1, x=yvals)

        ## FytX contains a distribution function for every value of t and x
        ## this step averages over all the x's
        Fyt <- lapply(1:length(tvals), function(i) {
            BMisc::combineDfs(yvals, Fytx[,i])
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
        
        out <- list(cop.param=delt[1], copula=copula, Fytxlist=Fytx, Fyt=Fyt, tvals=tvals, x=xdf)

    ### only do above if you want the results for a particular value of t and x;
    ### otherwise can just return all results 
    } else {
        ## out <- list(cop.param=parms2coppar(res$maximum, copula=copula, x=x),
        ##             copula=copula)
        out <- list(cop.param=delt[1], copula=copula)
    }

    out$Qyx <- Qyx
    out$Qtx <- Qtx

    class(out) <- "qr2meobj"

    return(out)
}

## check if this works, I think it is for putting back together an
## unconditional distribution from a list of them, but need to step through
## it
avgDist <- function(Fytxlist, yvals) {
    
    Fyt <- lapply(Fytxlist, function(FFytx) {
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
    ##qq <- t(sapply(getListElement(obj$Fyt, whichone),
    ##               function(FF) quantile(FF, probs=tau)))
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

addplot <- function(obj, p, whichone=1, tau=c(.1,.5,.9)) {

    ## qq <- t(sapply(getListElement(obj$Fytxlist, 1),
    ##                function(FF) quantile(FF, probs=tau)))
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


