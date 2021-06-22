#-----------------------------------------------------------------------------
# Code related to copulas for \code{qmre} package
#-----------------------------------------------------------------------------

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

  # previous versions of code allowed for copula parameter
  # to vary by covariates, but our assumptions imply that
  # the conditional copula is the same across all covariates
  # this is a way to enforce that restriction without
  # re-writing old code.
  intonly <- as.matrix(rep(1, nrow(x)))

  # in order to estimate parameters across a variety of copulas,
  # we restrict the value of the parameters to be between 0 and
  # 1, this function converts the paramters into the correct scale
  # for the copula parameter itself
  delt <- parms2coppar(params, copula, intonly)

  lval <- sapply(1:n, function(i) {
    max(
      mean(
        cop.pdf(Fyx[[i]](y[i] - Us), Ftx[[i]](t[i] - Vs), type=copula,
                delt=delt[i])  *
          fyx[[i]](y[i] - Us) * ftx[[i]](t[i] - Vs)
      )
    , eps) # here, eps avoids this being exactly equal to 0
  })
  sum(log(lval))
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

#-----------------------------------------------------------------------------
# don't think this is called anywhere
# comment for now, delete later
#-----------------------------------------------------------------------------
## getCBounds <- function(copula) {
##   if (copula=="gumbel") {
##     lo <- 1
##     hi <- 1000 ##Inf
##   } else if (copula=="frank") {
##     lo <- -1000 ##-Inf
##     hi <- 1000 ##Inf
##   } else if (copula=="clayton") {
##     lo <- -1000 ##-Inf
##     hi <- 1000 ##Inf
##   } else if (copula=="gaussian") {
##     lo <- -1
##     hi <- 1
##   } else {
##     stop( paste0("copula ", copula, " not supported") )
##   }
##   return(c(lo,hi))
## }
