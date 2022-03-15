# unusablesensors.R
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
# 20220116  GvB           Initial setup v0.3-2
#---------------------------------------------------------------------------------------------------------------------

#' unusablesensors
#' 
#' Detect unusable sensors containing NA or constant values
#' 
#' Constant values are defined to have a standard deviation or a median absolute
#' deviation smaller than 10e-10.
#' 
#' @param x input data, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}.
#' 
#' @return a named list:
#'   \describe{
#'     \item{n}{number of unusable channels}
#'     \item{unusable}{a named integer vector containing all unusable sensors}
#'     \item{sensors_NA}{a named integer vector containing sensors with
#'       missing values}
#'     \item{sensors_constant}{a named integer vector containing sensors
#'       with constant values}
#'   }
#' 
#' @examples
#' bad <- unusablesensors(EEGdata)
#' if (bad$n > 0) good <- EEGdata[, -bad$unusable]
#'  
#' @author Original Matlab code by Tim Mullen (2011); ported to R and
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com},
#'   
#' @export

unusablesensors <- function(x) {
  
  if (!("ctd" %in% class(x) || is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex matrix or vector")
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  npts <- nrow(x)
  ns <- ncol(x)
  
  sensor_NA <- which(apply(x, 2, function(x) any(is.na(x))))
  sensor_constant <- 
    which(apply(x, 2, function(x) stats::mad(x) < 10e-10 ||
                                  stats::sd(x) < 10e-10))
  unusable <- union(sensor_NA, sensor_constant)
  n <- length(unusable)

  list(n = n, unusable = unusable, sensor_NA = sensor_NA, sensor_constant = sensor_constant)
}
