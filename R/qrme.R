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
        lu <- -log(u)
        lv <- -log(v)
        out <- exp(-(lu^delt + lv^delt)^(1/delt)) * ( (lu^delt + lv^delt)^(1/delt) + delt - 1) * (lu^delt + lv^delt)^(1/delt - 2) * (lu*lv)^(delt-1) * (u*v)^(-1)
    } else if (type=="frank") {
        fc <- copula::frankCopula(delt)
        out <- copula::dCopula(cbind(u,v), fc)
    } else {
        stop(paste0("copula: ", type, " not supported"))
    }
    out <- max(out, eps)
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
ll <- function(params, y, t, x, copula="gumbel", Fyx, Ftx, fyx, ftx, Upi, Umu, Usig, Vpi, Vmu, Vsig, ndraws=100, eps=1e-300) {
    k <- length(params)
    params <- as.matrix(params)
    x <- as.matrix(x)
    n <- nrow(x)
    if (copula == "gumbel") {
        delt <- 1 + exp(x%*%params)
    } else if (copula == "frank") {
        delt <- x%*%params
    }
    ## make draws from the mixture distribution
    ksig <- length(Usig)
    Ucomponents <- sample(1:ksig, ndraws, replace=TRUE, prob=Upi)
    Us <- rnorm(ndraws, Umu[Ucomponents], Usig[Ucomponents])
    Vcomponents <- sample(1:ksig, ndraws, replace=TRUE, prob=Vpi)
    Vs <- rnorm(ndraws, Vmu[Vcomponents], Vsig[Vcomponents])

    lval <- sapply(1:n, function(i) {
        max( mean(cop.pdf(Fyx[[i]](y[i] - Us), Ftx[[i]](t[i] - Vs), type=copula,
                delt=delt[i]) *
                fyx[[i]](y[i] - Us) * ftx[[i]](t[i] - Vs)), eps) ## here, eps
        ## avoids this being exactly equal to 0
        })
    -sum(log(lval))
}



#' @title llme
#'
#' @description Log likelihood function for estimating quantile regression model with measurement
#'  error using the approach in Hausman, Liu, Luo, and Palmer (2016).
#' @inheritParams ll
#' @param ksig the number of mixture components
#' @param kmu the number of means in the mixture model (this should almost always be ksig-1)
#'  as the last mean is pinned down by the other means under the condition that the mean
#'  of the measurement error is equal to 0.
#'
#' @return negative value of log likelihood function
#'
#' @export
llme <- function(params, y, x, tau, ksig, kmu=(ksig-1)) {
    x <- as.matrix(x)
    kx <- ncol(x)*length(tau)
    bet <- params[1:kx]
    k <- kx/length(tau)
    n <- nrow(x)
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

    return(-sum(ll)) ## for maximization 
}

#' @title qrme
#'
#' @description Quantile Regression with measurement error in the dependent
#'  variable using the approach in Hausman, Liu, Luo, and Palmer (2016)
#'
#' @param formla y ~ x
#' @param tau vector for which quantiles to compute quantile regression
qrme <- function(formla, tau=0.5, data, nmix=3, startbet=NULL, startmu=NULL, startsig=NULL, startpi=NULL, method="BFGS", maxit=1000) {
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
        betvals <- unlist(lapply(1:length(tau), function(i) c(quantile(dta$lincF, probs=tau[i], type=1), rep(0, (k-1)))))
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

    out <- optim(qrparams, llme, method=method, y=y, x=x, tau=tau, ksig=m, control=list(maxit=maxit))

    params <- getParams(out, formla, data, tau, nmix)

    makeRQS(params, formla, data)
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
    bet <- optout$par[1:kx]
    k <- kx/length(tau)
    n <- nrow(x)
    ktau <- length(tau)
    bet <- split(bet,ceiling(seq_along(bet)/k))
    kmu <- nmix-1
    if (nmix > 1) {
        pi1 <- optout$par[(kx+1):(kx+kmu)]
        mu1 <- optout$par[(kx+kmu+1):(kx+kmu+kmu)]
        pi <- c(pi1, 1-sum(pi1))
        mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    } else {
        pi <- 1
        mu <- 0
    }
    ksig <- nmix
    sig <- optout$par[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]
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
    class(optout) <- "rqs"
    optout$terms <- terms(formla)
    optout$formula <- formla
    bet <- params$bet
    ##bet <- split(bet1,ceiling(seq_along(bet1)/k))
    ## rearrangement (as in HLLP though double check)
    barx <- apply(x, 2, mean)
    betorder <- order(sapply(bet, function(b) sum(b*barx)))
    bet <- t(simplify2array(bet))
    bet <- bet[betorder,]
    optout$coefficients <- t(bet)
    optout$tau <- tau
    optout$params <- params
    optout
}

#' @title parms2coppar
#'
#' @description convert parameters from optimization routine into copula
#'  parameters
#'
#' @inheritParams ll
#' @param params a vector of parameters
#'
#' @return a vector of copula parameters
#'
#' @export
parms2coppar <- function(params, copula="gumbel",x) {
    x <- as.matrix(x)
    if (ncol(x)==1) {
        x <- t(x)
    }
    if (copula=="gumbel") {
        deltx <- 1 + exp(x%*%params)
    } else if (copula=="frank") {
        deltx <- x%*%params
    } else {
        stop( paste0("copula ", copula, " not supported") )
    }
    deltx
}


qr2me <- function(yname, tname, xformla, tau, data, xdf=NULL, tvals=NULL,
                  copula="gumbel", Qyx, Qtx, Upi, Umu, Usig, Vpi, Vmu, Vsig,
                  startparams=NULL, method="BFGS", maxit=1000, ndraws=1000) {

    cat("\nqr2me method...")
    cat("\nCite: Callaway, Brantly, Tong Li, and Irina Murtazashvili (2018)...")
    x <- model.matrix(xformla, data)
    if (is.null(startparams)) {
        startparams <- rep(0, ncol(x))
    }

    Fyx <- predict(Qyx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)
    Ftx <- predict(Qtx, newdata=as.data.frame(x), type="Fhat", stepfun=TRUE)

    Fyx <- lapply(Fyx, rearrange)
    Ftx <- lapply(Ftx, rearrange)

    fyx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")
    ftx <- predict(Qyx, newdata=as.data.frame(x), type="fhat")

    
    ##ll(startparams, data[,yname], data[,tname], x=x, copula=copula,
    ##   Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx, Upi=Upi, Umu=Umu, Usig=Usig,
    ##   Vpi=Vpi, Vmu=Vmu, Vsig=Vsig)


    cat("\nStep 1 of 2: Estimating copula parameter...")
    res <- optim(startparams, ll, method=method,
                 y=data[,yname], t=data[,tname], x=x, copula=copula,
                 Fyx=Fyx, Ftx=Ftx, fyx=fyx, ftx=ftx,
                 Upi=Upi, Umu=Umu, Usig=Usig,
                 Vpi=Upi, Vmu=Umu, Vsig=Usig,
                 control=list(maxit=maxit))

    if (!is.null(tvals)) {

        if (is.null(xdf)) xdf <- x
 
        ##xtdf <- cbind(tvals, xdf)
        ## todo, this gives the copula, now convert to conditional distribution
        U <- runif(ndraws)
        delt <- parms2coppar(res$par, copula, xdf)
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

        cat("\nStep 2 of 2: Building conditional distributions...\n")
        Fytx <- pblapply(1:length(QQyx), function(i) {
            qfun <- QQyx[[i]]
            ffun <- Ftx[[i]]
            lapply(tvals, function(tt) {
                BMisc::makeDist(yvals, sapply(yvals, function(yy) {
                    if (copula=="frank") {
                        cop <- copula::frankCopula(as.numeric(delt[i]))
                    } else if (copula=="gumbel") {
                        cop <- copula::gumbelCopula(as.numeric(delt[i]))
                    } else {
                        stop(paste0("copula type:", copula, " is not supported"))
                    }
                    mean(1*(qfun(U)<=yy)*dCopula(cbind(U, ffun(tt)), cop))
                }))
            }) ## might want to make this a function with makeRQS (modified) or makeDist eventually
        }, cl=4)
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
            

        Fyt <- lapply(Fytx, function(FFytx) {
            combineDfs(yvals, FFytx)
        })

        out <- list(parms=parms2coppar(res$par, copula=copula, x=x),
                    Fyt=Fyt, Fytxlist=Fytx)
        
        ##QytxOut <- list(Qytx=Qytx, x=xdf, t=tvals, tau=tau, delt=delt)

### only do above if you want the results for a particular value of t and x;
    ### otherwise can just return all results 
    } else {
        out <- parms2coppar(res$par, copula=copula, x=x) 
    }

    return(out)
}
