% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removetrend.R
\name{removetrend}
\alias{removetrend}
\title{removetrend}
\usage{
removetrend(
  x,
  fs = 1,
  type = c("highpass", "linear", "mean", "constant", "none"),
  hpcf = 1,
  hptbw = 1,
  hpdev = 0.01
)
}
\arguments{
\item{x}{input time series, specified as a numeric matrix or vector. In case
of a vector it represents a single signal; in case of a matrix each column
is a signal. Alternatively, an object of class \code{\link{ctd}}}

\item{fs}{sampling frequency of \code{x} in Hz. Default: 1. Overruled if
\code{x} is a \code{ctd} object, in which case the sampling frequency is
\code{fs(x)}. Only used when \code{type = "highpass"}.}

\item{type}{character string indicating which type of trend removal
is performed for line noise removal bad channel detection; one of: 
\describe{
  \item{\code{highpass}}{use a high-pass FIR filter with cutoff frequency
   specified by the parameter \code{hpfc} (default)}
  \item{\code{linear}}{remove polynomial linear trend}
  \item{\code{mean|constant}}{remove mean}
  \item{\code{none}}{no detrending}
}}

\item{hpcf}{high-pass cutoff frequency in Hz, specified as a positive numeric
value. Default: 1 Hz}

\item{hptbw}{high-pass transition bandwith in Hz, specified as a positive
numeric value. For instance, if the high-pass cutoff frequency is 1 Hz, and
the transition bandwith is also 1 Hz, the transition band of 1 Hz is
located between 0.5 and 1.5 Hz. Default: 1 Hz.}

\item{hpdev}{deviation from desired stop- and passband of the high-pass filter,
specified as a positive numeric value}
}
\value{
a list containing consisting of:
  \describe{
    \item{\code{y}}{the detrended data}
    \item{\code{h}}{the filter coefficients }
    \item{\code{mean|constant}}{remove mean}
    \item{\code{none}}{no detrending}
  }
}
\description{
Perform detrending or high pass filtering to remove low frequencies
}
\details{
When high-pass filtering is requested, an optimal-length FIR filter is
designed, which has a linear phase. The filger delay is corrected for by the
function. Before filtering, reflected parts of the input data (equal to the
calculated filter length) are added to the beginning and end of the data to
avoid filter edge effects.
}
\examples{

data(EEGdata)
det <- removetrend(EEGdata, hpcf = 1, hptbw = 1, hpdev = 0.01)
plot(det$y)

}
\author{
Geert van Boxtel, \email{G.J.M.vanBoxtel@gmail.com}
}
