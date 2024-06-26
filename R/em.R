#-----------------------------------------------------------------------------
# code for running EM algorithm to estimate distribution of outcome conditional on covariates
##-----------------------------------------------------------------------------


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
#' @param simstep Whether to use MH in EM algorithm or importance sampling
#'  in EM algorithm.  "MH" for MH, and "ImpSamp" for importance sampling.
#'  Default is MH.
#' @param tol This is the convergence criteria.  When the change in the
#'  Euclidean distance between the new parameters (at each iteration) and
#'  the old parameters (from the previous iteration) is smaller than tol,
#'  the algorithm concludes.  In general, larger values for tol will result
#'  in a fewer number of iterations and smaller values will result in more
#'  accurate estimates.
#' @inheritParams qrme
#' @return QRME object
#'
#' @export
em.algo <- function(formla, data,
                    betmatguess, tau, m=1, piguess=1, muguess=0,
                    sigguess=1, simstep="MH", tol=.01,
                    iters=400, burnin=200, drawsd=4, cl=1,
                    messages=FALSE) {
    
    stopIters <- 100
    counter <- 1
    
    while (counter <= stopIters) {
        newone <- em.algo.inner(formla, data,
                                betmatguess, tau, m, piguess, muguess, sigguess, simstep, iters=iters, burnin=burnin, drawsd=drawsd, cl=cl, messages=messages)
        newbet <- newone$bet
        newpi <- newone$pi
        newmu <- newone$mu
        newsig <- newone$sig

        if (messages) {
            cat("\n\n\nIteration: ", counter, "\n pi: ", newpi, "\n mu: ", newmu, "\n sig: ", newsig, "\n\n")
        }
        criteria <- sqrt(sum(c(newbet-betmatguess,newsig-sigguess,newpi-piguess,newmu-muguess)^2))
        if (messages) {
            cat(" convergence criteria: ", criteria, "\n\n")
        }
        if ( criteria <= tol) { ## Euclidean norm
            if (messages) cat("\n algorithm converged\n") 
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
                          simstep="MH",
                          iters=400, burnin=200, drawsd=4, cl=1, messages=FALSE) {
    

    xformla <- formla
    xformla[[2]] <- NULL # drop y variable
    X <- model.matrix(xformla, data)
    yname <- as.character(formla[[2]])
    Y <- data[,yname]
    k <- ncol(X) # number of x variables
    n <- length(Y)
    if (messages) {
        cat("\nSimulating measurement error...")
    }

    if (simstep=="MH") {
        startval <- 0

        edraws <- mh_mcmcC(Y, X, startval=startval, iters=iters, burnin=burnin,
                           drawsd=drawsd, betmat=betmat, m=m,
                           pi=pi, mu=mu, sig=sig, tau=tau)

        newids <- unlist(lapply(1:n, function(i) rep(i, (iters-burnin)))) # just replicates Y and X over and over

        newdta1 <- as.data.frame(cbind(Y=(Y[newids]-edraws), X=X[newids,], e=edraws))

        # some extra debugging code for acceptance ratio of mh algorithm
        # newdta1$id <- sapply(strsplit(rownames(newdta1), split="[.]"), function(cc) cc[1])
        # ib <- iters-burnin
        # acceptanceratio <- sapply(split(newdta1, f=newdta1$id), function(df) 1-mean(df$e[1:(ib-1)] == df$e[2:ib]))
        # hist(acceptanceratio)
        
        # Note: this currently just works for one X; will need to update
        if (messages) {
            cat("\nEstimating QR including simulated measurement error...")
        }
        #newdta1 <- do.call(rbind.data.frame, newdta)
        colnames(newdta1) <- c(yname, colnames(X), "e")

        #thetime <- Sys.time()
        #out  <- quantreg::rq(formla, tau=tau, data=newdta1, method="fn")
        #Sys.time() - thetime
        #out <- quantreg::rq(formla, tau=tau, data=newdta1, method="pfn", weights=rep(1,nrow(newdta1)))
        #out <- rq.fit.pfnb(x=model.matrix(formla, data=data), 
        #              y=model.response(model.frame(formla, data=data)), 
        #              tau=tau
        newdta1$w <- 1
        out <- quantreg::rq(formula=formla, 
                            tau=tau, 
                            weights=w,
                            data=newdta1, 
                            method="pfn")
        # this is part I am not sure about, once you have a new beta then estimate a new sigma??
        # also should probably restrict overall mean of measurement error term to be equal to 0
        if (messages) {
            cat("\nEstimating finite mixture model...")
        }
        if (m == 1) {
            nm <- list(m=1, lambda=1, mu=0, sigma=sd(newdta1$e))
        } else {
            nm <- mixtools::normalmixEM(newdta1$e, k=m, epsilon=1e-03)
        }
        nmorder <- order(nm$mu)# reorder results by mean of each component

        return(list(bet=t(coef(out)), m=m, pi=nm$lambda[nmorder], mu=nm$mu[nmorder], sig=nm$sigma[nmorder]))#, Ystar=newdta1[,yname]))
        
    } else if (simstep=="ImpSamp") {
        # importance sampling
        edraws <- rnorm((iters*n), 0, drawsd)

        newids <- unlist(lapply(1:n, function(i) rep(i, iters))) # just replicates Y and X over and over

        newdta1 <- as.data.frame(cbind(Y=(Y[newids]-edraws), X=X[newids,], e=edraws))  # prepopulate some fields in dataset

        # compute weights using importance sampling
        newdta1$w <- imp_sampC(Y=Y[newids], X=X[newids,], V=edraws, iters=iters, drawsd=drawsd,
                              betmat=betmat, m=m, pi=pi, mu=mu, sig=sig, tau=tau) # but use original versions of the data (not adjusted for measurement errors) to compute weights
        newdta1$w  <- sapply(1:length(edraws), function(i) max(1e-05,newdta1$w[i])) # drop negative weights (not many of these...)
        # run weighted quantile regression
        out <- quantreg::rq(formla, tau=tau, data=newdta1, method="fn", weights=newdta1$w)
        if (messages) {
            cat("\nEstimating finite mixture model...")
        }
        ##
        # need to make actual draws here
        ## 
        # set up newids again
        #this is going to have a different length from
        # newids above because we are going to use the length of tau rather than number of iterations
        newids <- unlist(lapply(1:n, function(i) rep(i, length(tau)))) 
        # draws of Y^* using X'\beta(U) and having U take all possible values of tau 
        Ystar  <- c(t(X%*%coef(out)))
        # finally, recover draws of measurement error
        U  <- Y[newids] - Ystar
        if (m == 1) {
            nm <- list(m=1, lambda=1, mu=0, sigma=sd(U))
        } else {
            nm <- mixtools::normalmixEM(U, k=m, epsilon=1e-03)
        }
        nmorder <- order(nm$mu)# reorder results by mean of each component
        return(list(bet=t(coef(out)), m=m, pi=nm$lambda[nmorder], mu=nm$mu[nmorder], sig=nm$sigma[nmorder]))#, Ystar=Ystar))
    } else {
        stop("provided simstep not supported")
    }
}


##############################################################
##
## Code for computing densities used in MH algorithm
##
## Note: most of this has been replaced by C++ code to
## improve performance, but keeping it here to have R
## versions of these.
##
##############################################################
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
#' @param X An nxK matrix with n the number of observations of X
#' @param tau an L-vector containing the quantile at which Q(Y|X) was estimated
#'
#' @return An nx1 vector that contains f(y|X)
#'
#' @export
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
