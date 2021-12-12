# changefs.R
# Copyright (C) 2020  Geert van Boxtel, <G.J.M.vanBoxtel@gmail.com>
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
# 20201231  GvB           Initial setup v0.3-0
#---------------------------------------------------------------------------------------------------------------------

#' Change Sampling Rate
#' 
#' Change sampling rate of a \code{ctd} object.
#' 
#' Changing the sampling rate is done using the \code{resample} function in the
#' \code{gsignal} package.
#' 
#' @param ctd_obj object of class \code{'ctd'} containing continuous time domain
#'   data.
#' @param new_fs positive scalar indicating the new sampling frequency.
#'
#' @return A \code{ctd} object containing the resampled data.
#' 
#' @examples
#' data(EEGdata)
#' down <- changefs(EEGdata, 100)
#' up <- changefs(EEGdata, 512)
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @seealso \code{\link[gsignal]{resample}}
#' 
#' @export


changefs <- function(ctd_obj, new_fs) {
  
  # parameter checking
  if (!("ctd" %in% class(ctd_obj) || !is.matrix(ctd_obj))) {
    stop("ctd_obj must be a continuous time domain data object (matrix)")
  }
  
  if(!is.numeric(new_fs) || length(new_fs) != 1 || new_fs <= 0) {
    stop("new_fs must be a numeric value > 0")
  }
  
  # compute p and q
  old_fs <- eegr::fs(ctd_obj)
  pq <- numbers::contFrac(new_fs / old_fs)$rat

  # resample and make it a ctd object
  y <- gsignal::resample(ctd_obj, pq[1], pq[2])
  colnames(y) <- colnames(ctd_obj)
  class(y) <- c("ctd", "matrix")
  attr(y, "npts") <- as.numeric(nrow(y))
  ns <- ncol(y)
  if ("t" %in% colnames(y)) ns <- ns - 1
  attr(y, "ns") <- as.numeric(ns)
  attr(y, "fs") <- as.numeric(new_fs)
  
  # adapt t if necessary
  if ("t" %in% colnames(y)) {
    y[, "t"] <- (0:(nrow(y) - 1)) / rep(new_fs, nrow(y))
  }
  
  y
}
