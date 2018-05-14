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
cop.pdf <- function(u,v,type="gumbel",delt) {
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
ll <- function(params, y, t, x, copula="gumbel", Fyx, Ftx, fyx, ftx, Upi, Umu, Usig, Vpi, Vmu, Vsig, ndraws=100) {
    k <- length(params)
    params <- as.matrix(params)
    x <- as.matrix(x)
    n <- nrow(x)
    if (copula == "gumbel") {
        delt <- 1 + exp(x%*%params)
    } else if (copula == "frank") {
        delt <- x%*%params
    }
    ksig <- length(Usig)
    Ucomponents <- sample(1:ksig, ndraws, replace=TRUE, prob=Upi)
    Us <- rnorm(ndraws, Umu[Ucomponents], Usig[Ucomponents])
    Vcomponents <- sample(1:ksig, ndraws, replace=TRUE, prob=Vpi)
    Vs <- rnorm(ndraws, Vmu[Vcomponents], Vsig[Vcomponents])

    lval <- sapply(1:n, function(i) {
        mean(cop.pdf(Fyx[[i]](y[i] - Us), Ftx[[i]](t[i] - Vs), type=copula,
                delt=delt[i]) *
            fyx[[i]](y[i] - Us) * ftx[[i]](t[i] - Vs))
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
getParams <- function(optout, ksig) {
    x <- as.matrix(x)
    kx <- ncol(x)*length(tau)
    bet <- optout$par[1:kx]
    k <- kx/length(tau)
    n <- nrow(x)
    ktau <- length(tau)
    bet <- split(bet,ceiling(seq_along(bet)/k))
    kmu <- ksig-1
    if (ksig > 1) {
        pi1 <- optout$par[(kx+1):(kx+kmu)]
        mu1 <- optout$par[(kx+kmu+1):(kx+kmu+kmu)]
        pi <- c(pi1, 1-sum(pi1))
        mu <- c(mu1, -sum(mu1*pi1)/(1-sum(pi1)))
    } else {
        pi <- 1
        mu <- 0
    }
    sig <- optout$par[(kx+kmu+kmu+1):(kx+kmu+kmu+ksig)]
    list(bet=bet, pi=pi, mu=mu, sig=sig)
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
makeRQS <- function(optout) {
    class(optout) <- "rqs"
    optout$terms <- qrfout$terms ## some of this needs to be "more general"
    optout$formula <- lincS ~ ageF + ageS ## here too
    bet1 <- optout$par[1:length(betvals)]
    bet <- split(bet1,ceiling(seq_along(bet1)/k))
    ## rearrangement (as in HLLP though double check)
    barx <- apply(x, 2, mean)
    betorder <- order(sapply(bet, function(b) sum(b*barx)))
    bet <- t(simplify2array(bet))
    bet <- bet[betorder,]
    optout$coefficients <- t(bet)
    optout$tau <- tau
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
parms2coppar <- function(params, copula="gumbel") {
    if (copula=="gumbel") {
        deltx <- 1 + exp(x%*%params)
    } else if (copula=="frank") {
        deltx <- x%*%params
    } else {
        stop( paste0("copula ", copula, " not supported") )
    }
    deltx
}
