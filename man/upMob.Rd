% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intMob.R
\name{upMob}
\alias{upMob}
\title{upMob}
\usage{
upMob(Y, T, amount = 0, qcutoffs = c(0.25, 0.5, 0.75))
}
\arguments{
\item{Y}{vector of outcomes}

\item{T}{vector of treatments}

\item{amount}{How much child's rank has to increase in order to count
for upward mobility}

\item{qcutoffs}{the quantiles where to draw the cutoff points in the transition matrix; default is at the quartiles}
}
\value{
vector of amount of upward mobility by parents' income quartile
 (or other specified cutoffs)
}
\description{
Upward Mobility Paramters
}
