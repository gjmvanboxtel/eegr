# rerefdata.R
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
# 20210103  GvB           Initial setup v0.3-0
#---------------------------------------------------------------------------------------------------------------------

#' Rereference Data
#' 
#' Compute the average reference or a new common reference for data in a
#' \code{ctd} object.
#' 
#' @param ctd_obj object of class \code{'ctd'} containing continuous time domain
#'   data.
#' @param signals signals to compute the new reference on, either specified as a
#'   character vector of signal names (\code{colnames(ctd_obj)}, or as a vector
#'   of numerical values.). Alternatively, specifying \code{"all"} (default)
#'   computes the new reference of all signals, except the signal(s) specified
#'   in \code{ref}.
#' @param new_ref signals that form the new reference, either specified as a
#'   character vector of signal names (\code{colnames(ctd_obj)}), or as a vector
#'   of numerical values. If more than one signal is specified, then the average
#'   of these sinals will form the new reference. Alternatively, specifying
#'   \code{"average"} (default) computes the average reference over the signals
#'   specified by \code{signals}.
#' @param keep_ref logical indicating whether to keep the reference in the output.
#'   Default: \code{FALSE}.
#' @param old_ref list with the sensor locations, either ('label', 'theta',
#'   'phi'), or ('label, 'x', 'y', 'z'), specifying the location of the original
#'   reference. If specified, the old reference signal is reconstructed back
#'   into the data object.
#'
#' @return A \code{ctd} object containing the re-referenced data. Signals not
#'   included neither in \code{signals} nor in \code{ref} are copied to the
#'   output object.
#' 
#' @examples
#' 
#' 
#' #data(EEGdata)
#' #new <- rerefdata(EEGdata, 'all', 'average')
#' #plot(new)
#' #new <- rerefdata(EEGdata, 'all', c('M1', 'M2'))
#' 
#' @seealso \code{\link{sensorlocs}}
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @export


rerefdata <- function(ctd_obj, signals, new_ref, keep_ref = FALSE, old_ref) {

  stop("Function not implemented yet")
  
  # parameter checking
  if (!("ctd" %in% class(ctd_obj) || is.matrix(ctd_obj))) {
    stop("ctd_obj must be a continuous time domain data object (matrix)")
  }
  if (!(signals == "all" || 
       (is.numeric(signals) && max(signals) <= ns(ctd_obj)) ||
       (is.character(signals) && all(signals %in% colnames(ctd_obj))))) {
    stop('Invalid signals specified')
  }
  if (!(new_ref == "average" || 
        (is.numeric(new_ref) && max(new_ref) <= ns(ctd_obj)) ||
        (is.character(new_ref) && all(new_ref %in% colnames(ctd_obj))))) {
    stop('Invalid new reference specified')
  }
  if(!is.logical(keep_ref)) {
    stop('keep_ref must be TRUE or FALSE')
  }
  if(!missing(old_ref)) {
    #if(!)
  }

  #y
  ctd_obj
}
