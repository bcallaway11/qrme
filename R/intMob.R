#-----------------------------------------------------------------------------
# code for computing intergenerational mobility parameters
#-----------------------------------------------------------------------------


#' @title tmat
#'
#' @description A function for computing transition matrices
#'
#' @param Y vector of outcomes
#' @param T vector of treatments
#' @param qcutoffs the quantiles where to draw the cutoff points in the transition matrix; default is at the quartiles
#'
#' @return A matrix containing values of the transition matrix.  The columns contain different groups by value of the treatment,
#'  the rows contain different values for the outcome (it is ordered from bottom to top though as is standard for a transition matrix)
tmat <- function(Y,T, qcutoffs=c(.25,.5,.75)) {
    ucutoffs <- c(qcutoffs,1)
    lcutoffs <- c(0,qcutoffs)

    gc <- copula::empCopula(copula::pobs(cbind(Y,T)))
    
    ucutoffs <- c(qcutoffs,1)
    lcutoffs <- c(0,qcutoffs)

    tmatcell <- function(r1, r2, s1, s2, cop) { ## r1, r2 are for Y, s1, s2 are for T
        c1 <- copula::pCopula(cbind(r2, s2), copula=cop)
        c2 <- copula::pCopula(cbind(r2, s1), copula=cop)
        c3 <- copula::pCopula(cbind(r1, s2), copula=cop)
        c4 <- copula::pCopula(cbind(r1, s1), copula=cop)
        denom <- s2 - s1
        (c1-c2-c3+c4)/denom
    }

    transmat <- matrix(nrow=length(ucutoffs), ncol=length(ucutoffs))
    for (i in 1:length(ucutoffs)) {
        for (j in 1:length(lcutoffs)) {
            transmat[i,j]  <- tmatcell(lcutoffs[i], ucutoffs[i], lcutoffs[j], ucutoffs[j],
                                       cop=gc) # rows should be for Y, cols for T
        }
    }
    apply(transmat, 2, rev) # transition matrix is going to be upside down, this fixes it.   
}


#' @title upMob
#'
#' @description Upward Mobility Paramters
#'
#' @inheritParams tmat
#' @param amount How much child's rank has to increase in order to count
#'  for upward mobility
#'
#' @return vector of amount of upward mobility by parents' income quartile
#'  (or other specified cutoffs)
#' @export
upMob <- function(Y, T, amount=0, qcutoffs=c(.25,.5,.75)) {
    cutoffs <- c(0, qcutoffs, 1)
    ranksY <- ecdf(Y)(Y)
    ranksT <- ecdf(T)(T)

    upMobcell <- function(amt, s1, s2) {
        theseObs <- which( (ranksT >= s1) & (ranksT <=s2))
        thisRanksY <- ranksY[theseObs]
        thisRanksT <- ranksT[theseObs]
        mean( 1*((thisRanksY - thisRanksT) > amt) )
    }

    out <- c()
    for (i in 1:(length(cutoffs)-1)) {
        out[i] <- upMobcell(amount, cutoffs[i], cutoffs[i+1])
    }
    out
}
