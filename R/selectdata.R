# selectdata.R
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
# 20210101  GvB           Initial setup v0.3-0
# 20220511  GvB           use inherits() instead of class(res)
#---------------------------------------------------------------------------------------------------------------------

#' Select Data
#' 
#' Select data from a \code{ctd} object.
#' 
#' @param ctd_obj object of class \code{'ctd'} containing continuous time domain
#'   data.
#' @param seconds character string representing an expression to select the
#'   specified seconds from the data; can only used when \code{ctd_obj} contains
#'   a \code{'t'} variable alongside the signals.
#' @param samples samples to select from the data; must be between 1 and
#'   \code{npts(ctd_obj)}.
#' @param signals signals to select from the data. Either a vector of numeric
#'   values between 1 and \code{ns(ctd_obj)}, or a character vector with sensors
#'   names from \code{colnames(ctd_obj)}.
#' @param add.t logical indicating whether to add or renew the time variable t.
#'   If \code{TRUE} and \code{ctd_obj} did not contain a variable \code{t}, it
#'   is added. If it did contain a variabe \code{t}, it is updated starting at
#'   0. If \code{FALSE} (default), an existing variable \code{t} is left
#'   unchanged.
#'
#' @return A \code{ctd} object containing the selected data.
#' 
#' @examples
#' 
#' data(EEGdata)
#' sel <- selectdata(EEGdata, seconds = 'which(t < 3)')
#' sel <- selectdata(EEGdata, samples = 1:600)
#' sel <- selectdata(EEGdata, signals = 1:28)
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @export

selectdata <- function(ctd_obj, seconds = NULL, samples = NULL, signals = NULL,
                       add.t = FALSE) {
  
  # parameter checking
  if (!("ctd" %in% class(ctd_obj) || is.matrix(ctd_obj))) {
    stop("ctd_obj must be a continuous time domain data object (matrix)")
  }
  p <- which(!is.null(c(seconds, samples, signals)))
  if (length(p) <= 0 || length(p) > 2) {
    stop('selection can be signals and either seconds or samples')
  }
  if (!is.null(seconds) && !is.null(samples)) {
    stop('selection must be either seconds or samples')
  }
  
  # initialize
  smp <- 1:nrow(ctd_obj)
  sig <- 1:ncol(ctd_obj)
  tpresent <- ('t' %in% colnames(ctd_obj))

  # select seconds
  if (!is.null(seconds)) {
    if (!tpresent) {
      stop('ctd_obj must contain a variable named t')
    }
    if (!is.character(seconds) || length(seconds) != 1) {
      stop('seconds must be a character vector representing an expression')
    }
    expr <- gsub('t', 'ctd_obj[, "t"]', seconds)
    res <- try(smp <- eval(parse(text = expr)))
    if (inherits(res, "try-error")) {
      stop('seconds contains an invalid expression')
    }
  }
  
  # select samples (adapt to npts(ctd_obj))
  if(!is.null(samples)) {
    if (is.numeric(samples)) {
      smp <- samples[which(samples > 0 & samples <= npts(ctd_obj))]
    } else {
      stop('samples must be a numeric vector')
    }
  }
  
  # select signals (adapt to ns(ctd_obj) or colnames(ctd_object))
  if (!is.null(signals)) {
    if (is.numeric(signals)) {
      sig <- signals[which(signals > 0 & signals < ns(ctd_obj))]
    } else if (is.character(signals)) {
      sig <- signals[which(signals %in% colnames(ctd_obj))]
    } else {
      stop('signals must be a numeric or character vector')
    }
  }

  # select samples and signals
  # handle t variable should it not be selected
  y <- ctd_obj[smp, sig]
  if(is.character(sig)) {
    colnames(y) <- sig
  } else if (is.numeric(sig)) {
    colnames(y) <- colnames(ctd_obj)[sig]
  }
  if (tpresent && !('t' %in% colnames(y))) {
    y <- cbind(y, ctd_obj[smp, "t"])
    colnames(y)[ncol(y)] <- 't'
  }
  
  # set attributes
  class(y) <- c("ctd", "matrix")
  npts <- as.numeric(nrow(y))
  attr(y, "npts") <- npts
  ns <- ncol(y)
  if ("t" %in% colnames(y)) ns <- ns - 1
  attr(y, "ns") <- as.numeric(ns)
  fs <- as.numeric(fs(ctd_obj))
  attr(y, "fs") <- fs

  # add t if requested
  if (add.t) {
    if (tpresent) {
      y[, "t"] <- (0:(npts - 1)) / rep(fs, npts)
    } else {
      y <- cbind(y, (0:(npts - 1)) / rep(fs, npts))
      colnames(y)[ncol(y)] <- 't'
    }
  }
  
  y
}
