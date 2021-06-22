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
