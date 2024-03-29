% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unusablesensors.R
\name{unusablesensors}
\alias{unusablesensors}
\title{unusablesensors}
\usage{
unusablesensors(x)
}
\arguments{
\item{x}{input data, specified as a numeric matrix or vector. In case
of a vector it represents a single signal; in case of a matrix each column
is a signal. Alternatively, an object of class \code{\link{ctd}}.}
}
\value{
a named list:
  \describe{
    \item{n}{number of unusable channels}
    \item{unusable}{a named integer vector containing all unusable sensors}
    \item{sensors_NA}{a named integer vector containing sensors with
      missing values}
    \item{sensors_constant}{a named integer vector containing sensors
      with constant values}
  }
}
\description{
Detect unusable sensors containing NA or constant values
}
\details{
Constant values are defined to have a standard deviation or a median absolute
deviation smaller than 10e-10.
}
\examples{
bad <- unusablesensors(EEGdata)
if (bad$n > 0) good <- EEGdata[, -bad$unusable]
 
}
\author{
Original Matlab code by Tim Mullen (2011); ported to R and
  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com},
}
