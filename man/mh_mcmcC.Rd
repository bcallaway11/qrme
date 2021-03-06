% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mh_mcmcC}
\alias{mh_mcmcC}
\title{mh_mcmcC}
\usage{
mh_mcmcC(Y, X, startval, iters, burnin, drawsd, betmat, m, pi, mu, sig, tau)
}
\arguments{
\item{Y}{vector of outcomes}

\item{X}{matrix of covariates}

\item{startval}{starting value for the markov chain}

\item{iters}{number of Monte Carlo iterations}

\item{burnin}{number of first MC iteration to drop}

\item{drawsd}{the standard deviation for the standard normal draws in the
MH algorithm}

\item{betmat}{matrix of QR parameters}

\item{m}{number of mixture components for measurement error}

\item{pi}{mixture probabilities}

\item{mu}{means of mixture components}

\item{sig}{standard deviations of mixture components}

\item{tau}{which values QR's have been estimated for}
}
\description{
mh_mcmcC
}
