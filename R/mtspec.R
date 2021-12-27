# mtspec.R
# Copyright (C) 2021  Geert van Boxtel, <G.J.M.vanBoxtel@gmail.com>
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
# 20210103  GvB           Initial setup v0.3-2
#---------------------------------------------------------------------------------------------------------------------

#' Multitaper Spectrum
#' 
#' Compute the multitaper spectrum of a vector or matrix.
#' 
#' The multitaper method is a method for spectral density estimation, which
#' overcomes some of the limitations of conventional Fourier analysis. It is a
#' periodogram-based method that uses multiple tapers, or windows, to form
#' independent estimates of the spectral density to reduce variance of the
#' spectral density estimate. The method uses Slepian or discrete prolate
#' spheroidal sequences as tapers since these vectors are mutually orthogonal
#' and possess desirable spectral concentration properties. They average out
#' noise in the spectrum and reduce information loss at the edges of the window.
#' 
#' The package \code{multitaper} contains a more complete version of multitaper
#' spectral analysis (function \code{\link[multitaper]{spec.mtm}}), which
#' computes confidence intervals, F-tests , and more.
#' 
#' The present implementation uses sliding windows, i.e., overlapping segments,
#' tapering each segment, and then averaging over tapered segments. The tapers
#' are computed by the function \code{\link[multitaper]{dpss}}, with parameters
#' \code{n = w * fs} (number of samples in sliding window), \code{k = tbw * w -
#' 1} (number of tapers), and \code{nw = tbw * w / 2} (time bandwidth of the
#' tapers).
#' 
#' The function \code{mtspec} computes the multitaper spectrum for a vector or
#' matrix. It calls the function \code{mtchan}, which computes the multitaper
#' spectrum for a single channel (no parameter checking).
#' 
#' @param x input time series, specified as a numeric or complex vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}
#' @param detrend character string specifying detrending option; one of:
#'   \describe{
#'     \item{\code{long-mean}}{remove the mean from the data before
#'     splitting into segments (default)}
#'     \item{\code{short-mean}}{remove the mean value of each segment}
#'     \item{\code{long-linear}}{remove linear trend from the data before
#'     splitting into segments}
#'     \item{\code{short-linear}}{remove linear trend from each segment}
#'     \item{\code{none}}{no detrending}
#'  }
#' @param w length of sliding window (segments) in seconds. Default: 4
#' @param overlap proportion of overlap between the segments. Default: 0
#' @param tbw taper bandwidth, specified as a positive numeric value. Default: 2
#'   Hz
#' @param k number of tapers to use, specified as a positive numeric value.
#'   Default: \code{tbw * w - 1}
#'
#' @return An object of class \code{frd}, containing the spectra
#' 
#' @examples
#' data(EEGdata)
#' sp <- mtspec(EEGdata[, 1:28], fs = fs(EEGdata), detrend = "short-linear")
#' plot(sp, yscale = "dB", main = "EEGdata")
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}, based on Matlab
#'   code by Tim Mullen's (2011) adaptation of the Chronux 2 toolbox
#'   (\url{http://chronux.org/}) for the cleanline plugin in EEGlab.
#' 
#' @references \url{https://en.wikipedia.org/wiki/Multitaper}
#' 
#' @seealso \link[multitaper]{spec.mtm}, \link[multitaper]{dpss},
#'   \code{\link{frd}}
#' 
#' @rdname mtspec
#' @export


mtspec <- function(x, fs = 1,
                   detrend = c("short-linear", "short-mean",
                               "long-linear", "long-mean", "none"),
                   w = 4, overlap = 0,
                   tbw = 2, k = round(tbw * w - 1)) {
  
  # parameter checking
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
    if ("t" %in% colnames(x)) {
      x <- x[, -which(colnames(x) == "t")]
    }
  } else {
    if (!(is.numeric(fs) && length(fs) == 1 && fs > 0)) {
      stop("fs must be a positive numeric value")
    }
  }
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
  }
  detrend <- match.arg(detrend)
  if (!(is.numeric(w) && length(w) == 1 && w > 0)) {
    stop("w must be a positive numeric value")
  }
  if (!(is.numeric(overlap) && length(overlap) == 1 &&
        overlap >= 0 && overlap < 1)) {
    stop("overlap must be a numeric value between 0 and 1")
  }
  if (!(is.numeric(tbw) && length(tbw) == 1 && tbw > 0)) {
    stop("tbw must be a positive numeric value")
  }
  k <- round(k)
  if (!(is.numeric(k) && length(k) == 1 && k > 0)) {
    stop("k must be a positive numeric value")
  }

  # global trend removal
  if (detrend == "long-mean") {
    x <- gsignal::detrend(x, p = 0)
  } else if (detrend == "long-linear") {
    x <- gsignal::detrend(x, p = 1)
  }
  
  # compute and scale tapers 
  n <- round(w * fs)
  n <- 2^ceiling(log2(n))
  nw = tbw * w / 2
  tapers <- multitaper::dpss(n, k, nw)$v * sqrt(fs)
  
  y <- apply(x, 2, mtchan, fs, detrend, n, overlap, tapers)
  
  if (n %% 2 == 0) {    # one-sided, nfft is even
    psd_len <- n / 2 + 1
    y <- apply(y, 2, 
               function(x)
                 x[1:psd_len] + c(0, x[seq(n, psd_len + 1, -1)], 0))
  } else {                    # one-sided, nfft is odd
    psd_len <- (n + 1) / 2
    y <- apply(y, 2,
               function(x)
                 x[1:psd_len] + c(0, x[seq(n, psd_len + 1, -1)]))
  }
  
  f <- seq.int(0, psd_len - 1) * (fs / n)
  frd(y, f, fu = "Hz", fs = fs, type = "spectrum", yscale = "linear")
}

#' @param chan single channel
#' @param n segment length (nearest power of 2 for \code{w})
#' @param tapers tapers from dpss
#' @rdname mtspec
#' @export

mtchan <- function(chan, fs, detrend, n, overlap, tapers) {
  len <- length(chan)
  stepsz <- n - round(n * overlap)
  S <- rep(0, n)
  nsegs <- 0
  for (sseg in seq(1, len - n + 1, stepsz)) {
    eseg <- min(len, sseg + n - 1)
    #cat(sseg, eseg, "\n")
    if (detrend == "short-mean") {
      seg <- gsignal::detrend(chan[sseg:eseg], p = 0)
    } else if (detrend == "short-linear") {
      seg <- gsignal::detrend(chan[sseg:eseg], p = 1)
    } else {
      seg <- chan[sseg:eseg]
    }
    mtf <- mtfft(seg, fs, tapers)
    S <- S + apply(Re(mtf * Conj(mtf)), 1, mean) # average over tapers
    nsegs <- nsegs + 1
  }
  if (nsegs > 1) S <- S / nsegs
  S
}


#' The function \code{mtfft} returns the multitaper Fourier tranform of a vector, and
#' returns a complex vector with the same number of columns as tapers.
#' 
#' @param dat data
#' 
#' @rdname mtspec
#' @export

mtfft <- function(dat, fs, tapers) {
  l <- length(dat)
  k <- ncol(tapers)
  n <- nrow(tapers)
  if (l != n) stop("length of tapers must equal length of data")
  stats::mvfft(dat * tapers) / fs
}
