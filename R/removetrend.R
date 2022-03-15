# removetrend.R
# Copyright (C) 2022  Geert van Boxtel, <G.J.M.vanBoxtel@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Version history:
# 202201123  GvB           Initial setup v0.3-2
#---------------------------------------------------------------------------------------------------------------------

#' removetrend
#' 
#' Perform detrending or high pass filtering to remove low frequencies
#' 
#' When high-pass filtering is requested, an optimal-length FIR filter is
#' designed, which has a linear phase. The filger delay is corrected for by the
#' function. Before filtering, reflected parts of the input data (equal to the
#' calculated filter length) are added to the beginning and end of the data to
#' avoid filter edge effects.
#' 
#' @param x input time series, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}. Only used when \code{type = "highpass"}.
#' @param type character string indicating which type of trend removal
#'   is performed for line noise removal bad channel detection; one of: 
#'   \describe{
#'     \item{\code{highpass}}{use a high-pass FIR filter with cutoff frequency
#'      specified by the parameter \code{hpfc} (default)}
#'     \item{\code{linear}}{remove polynomial linear trend}
#'     \item{\code{mean|constant}}{remove mean}
#'     \item{\code{none}}{no detrending}
#'   }
#' @param hpcf high-pass cutoff frequency in Hz, specified as a positive numeric
#'   value. Default: 1 Hz
#' @param hptbw high-pass transition bandwith in Hz, specified as a positive
#'   numeric value. For instance, if the high-pass cutoff frequency is 1 Hz, and
#'   the transition bandwith is also 1 Hz, the transition band of 1 Hz is
#'   located between 0.5 and 1.5 Hz. Default: 1 Hz.
#' @param hpdev deviation from desired stop- and passband of the high-pass filter,
#'   specified as a positive numeric value
#' 
#' @return a list containing consisting of:
#'   \describe{
#'     \item{\code{y}}{the detrended data}
#'     \item{\code{h}}{the filter coefficients }
#'     \item{\code{mean|constant}}{remove mean}
#'     \item{\code{none}}{no detrending}
#'   }
#'   
#' @examples
#' 
#' data(EEGdata)
#' det <- removetrend(EEGdata, hpcf = 1, hptbw = 1, hpdev = 0.01)
#' plot(det$y)
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @export

removetrend <- function(x, fs = 1,
                        type = c("highpass", "linear", "mean", "constant", "none"), 
                        hpcf = 1, hptbw = 1, hpdev = 0.01) {
  
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
  } else {
    if (!(is.numeric(fs) && length(fs) == 1 && fs > 0)) {
      stop("Sampling frequency (fs) must be a positive numeric value")
    }
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  npts <- nrow(x)

  type <- match.arg(type)
  if (type == "highpass") {
    if (!(is.numeric(hpcf) && length(hpcf) == 1 && hpcf > 0)) {
      stop("high-pass cutoff frequency (hpcf) must be a positive numeric value")
    }
    if (!(is.numeric(hptbw) && length(hptbw) == 1 && hptbw > 0)) {
      stop("high-pass transition bandwith (hptbw) must be a positive numeric value")
    } 
    if (!(is.numeric(hpdev) && length(hpdev) == 1 && hptbw > 0)) {
      stop("high-pass filter deviation (hpdev) must be a positive numeric value")
    } 
  }
 
  h <- NULL 
  if (type == "highpass") {
    ford <- gsignal::kaiserord(c(hpcf - hptbw / 2, hpcf + hptbw / 2),
                               c(0, 1), dev = hpdev, fs = fs)
    h <- gsignal::fir1(ford$n, ford$Wc, ford$type, gsignal::kaiser(ford$n + 1, ford$beta))
    delay <- round(ford$n / 2)
    # add reflected parts to both ends
    ext <- rbind(apply(x[1:delay, ], 2, rev),
                 x,
                 apply(x[(npts - delay + 1):npts, ], 2, rev)
    )
    y <- gsignal::filter(h, ext)[(ford$n + 1):(npts + ford$n), ]
  } else if (type == "linear") {
    y <- gsignal::detrend(x, p = 1)
  } else if (type == "mean" || type == "constant") {
    y <- gsignal::detrend(x, p = 0)
  } else {
    y <- x
  }
  # reset t variable if present (this is faster then subsetting x)
  if ("t" %in% colnames(x)) {
    y[, "t"] <- x[, "t"]
  }
  attributes(y) <- attributes(x)
  list(y = y, h = h)
}