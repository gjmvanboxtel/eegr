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
# 20210103  GvB           Initial setup v0.3-0
#---------------------------------------------------------------------------------------------------------------------

#' Detrend Data
#' 
#' Remove polynomial trend in a \code{ctd} object.
#' 
#' 
#' @examples
#' 
#' data(EEGdata)
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @export


dum <- function() {
  
  # # parameter checking
  # if (!("ctd" %in% class(ctd_obj) || is.matrix(ctd_obj))) {
  #   stop("ctd_obj must be a continuous time domain data object (matrix)")
  # }
  # 
  # if ("t" %in% colnames(ctd_obj)) {
  #   t <- ctd_obj[, "t"]
  #   y <- gsignal::detrend(selectdata(ctd_obj, signals = which(colnames(ctd_obj) != "t")), p)
  #   y <- cbind(y, t)
  # } else {
  #   y <- gsignal::detrend(ctd_obj, p)
  # }
  # attributes(y) <- attributes(ctd_obj)
  # 
  # y
}
