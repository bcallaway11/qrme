% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{interpolateC}
\alias{interpolateC}
\title{interpolateC}
\usage{
interpolateC(x, y, xval, extrapolate)
}
\arguments{
\item{x}{vector of x's}

\item{y}{vector of y's}

\item{xval}{a particular value to interpolate for}

\item{extrapolate}{whether or not to linearly extrapolate beyond endpoints
of x}
}
\value{
interpolated value
}
\description{
Returns interpolated value at x from parallel arrays ( xData, yData )
 Assumes that xData has at least two elements, is sorted and
is strictly monotonic increasing
}
