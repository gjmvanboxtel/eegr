# filterdata.R
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
# 20210102  GvB           Initial setup v0.3-0
#---------------------------------------------------------------------------------------------------------------------

#' Filter Data
#' 
#' Filter data in a \code{ctd} object.
#' 
#' Filtering is done using the \code{filtfilt} function in the \code{gsignal}
#' package, which minimizes filter startup transients by pre- and postpending
#' reflected pieces of the input signal, which are tapered to zero. Forward and
#' reverse filtering the signals corrects for phase distortion (not perfect in
#' practice).
#' 
#' @param ctd_obj object of class \code{'ctd'} containing continuous time domain
#'   data.
#' @param filt filter object recognized by the \code{gsignal} package (of class
#'   \code{Arma}, \code{Ma}, \code{Sos}, or \code{Zpg}, characterizing a FIR or
#'   a IIR filter.)
#'
#' @return A \code{ctd} object containing the filtered data.
#' 
#' @examples
#' 
#' data(EEGdata)
#' fs <- fs(EEGdata)
#' nyq <- fs / 2
#' 
#' ## 10 Hz low-pass (Butterworth; maximally flat)
#' but <- gsignal::butter(3, 10 / nyq)
#' lpd <- filterdata(EEGdata, but)
#' plot(lpd)
#' 
#' ## 8-12 Hz low-pass filter 
#' ## note that selecting a subset of signal
#' filt <- gsignal::fir1(50, c(8, 12) / nyq, "pass")
#' alpha <- filterdata(selectdata(EEGdata, signals = 1:28), filt)
#' plot(alpha, ylim=c(-10,10))
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @export


filterdata <- function(ctd_obj, filt) {
  
  # parameter checking
  if (!("ctd" %in% class(ctd_obj) || is.matrix(ctd_obj))) {
    stop("ctd_obj must be a continuous time domain data object (matrix)")
  }
  if (!class(filt) %in% c("Arma", "Ma", "Sos", "Zpg")) {
    stop("filt must be a filter object recognized by gsignal")
  }

  if ("t" %in% colnames(ctd_obj)) {
    t <- ctd_obj[, "t"]
    y <- gsignal::filtfilt(filt, selectdata(ctd_obj, signals = which(colnames(ctd_obj) != "t")))
    y <- cbind(y, t)
  } else {
    y <- gsignal::filtfilt(filt, ctd_obj)
  }
  attributes(y) <- attributes(ctd_obj)

  y
}
