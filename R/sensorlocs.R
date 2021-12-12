# sensorlocs.R
# define object and methods for sensorlocs (sensor locations)
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
# Version history
# 20160622  GvB       Initial version
# 20190906  GvB       Setup for eegr 0.1.0
#
#---------------------------------------------------------------------------------------------------------------------

#' Sensor locations
#' 
#' Define the \code{sensorlocs} object and associated methods.
#' 
#' A \code{\link{sensorlocs}} object is a \code{\link{data.frame}} containing
#' sensor locations that can be used in further analyses or display. In order to
#' match R's plotting functions, the coordinate system used in \code{eegr}
#' models the head as a perfect sphere with a radius of 1 and the origin in the
#' middle of the head (0,0,0), between the left and right pre-auricular fossae,
#' and the line between teh nasion and inion (Towle et al, 1993). Sensor labels
#' are as per the 10-5 electrode positioning system (Oostenveld & Praamstra, )
#' 
#' Within this system, Cartesian coordinates are defined along the following axes:
#' \describe{
#'   \item{x}{the X-axis points from the left (-) to right (+) pre-auricular
#'   points, i.e., positive values of x refer to sensors over the right
#'   hemisphere, whereas negative values of x refer to sensors over the left
#'   hemisphere.}
#'   \item{y}{the Y-axis points from the back of the head (inion, -) to the
#'   front (nasion, +), i.e., positive values of y refer to sensors over the
#'   anterior part of the brain, whereas negative values of y refer to sensors
#'   over the posterior part of the brain.}
#'   \item{z}{the Z-axis point from the bottom of the head (-) to the top (+) of
#'   body, i.e., positive values of z refer to sensors above the center of the
#'   head, whereas negative values of z refer to sensors below the center of the
#'   head.}
#' }
#' Angular coordinates are defined as follows:
#' \describe{
#'   \item{theta}{Theta is the angle angle between the Z-axis and the X/Y-plane;
#'   positive values indicate sensors over the right hemisphere; negative values
#'   over the left.}
#'   \item{phi}{Phi is the angle between the X- and Y-axes, with positive values
#'   indicating counterclockwise rotation and negative values clockwise.}
#' }
#' 
#' @references Oostenveld, R., Praamstra, P., 2001. The five percent electrode
#'   system for high-resolution EEG and ERP measurements. Clin. Neurophysiol.
#'   112, 713-719.
#' @references Towle, V.L., et al. (1993). The spatial location of EEG
#'   electrodes: locating the best-fitting sphere relative to cortical anatomy.
#'   Electroencephal. Clin Neurophysiol. 86, 1-6.
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @param x The input data, either a data frame or a matrix (coerced to data
#'   frame). It should either contain the variables x, y, z or theta and phi.
#'   The coordinates should already be in the eegr coordinate system, i.e., no
#'   conversion is done.
#' 
#' @rdname sensorlocs
#' @export

sensorlocs <- function (x, ...) UseMethod ("sensorlocs")

#' @return a data frame of class \code{'sensorlocs'} containing at least the
#'   variables x, y, z, x2d, y2d, theta, phi as per the eegr internal format
#'   described above. Other variables in the input data frame or matrix are
#'   copied unchanged to the output data frame. If a \code{label} variable was
#'   absent from the input \code{x}, then labels will be added based on the
#'   minimum spherical distance with electrodes in the \code{EEGlocations}
#'   dataset.
#'   
#' The calculation of the distance is based on the 'great circle distance' often
#' used in Earth distance calculations. For further information see
#' \url{https://en.wikipedia.org/wiki/Great-circle_distance}. There are 3
#' methods often used to compute the great circle distance, the law of cosines,
#' the haversine formula, and the Vincenty formula. The latter is most accurate
#' for small distances and is used in this function.
#' 
#' @references  Vincenty, T. (1975). Direct and Inverse Solutions of Geodesics
#'   on the Ellipsoid with Application of Nested Equations. Survey Review.
#'   Kingston Road, Tolworth, Surrey: Directorate of Overseas Surveys. 23 (176):
#'   88–93. doi: \url{https://dx.doi.org/10.1179/sre.1975.23.176.88}.
#' 
#' @rdname sensorlocs
#' @export

sensorlocs.default <- function (x, ...) {
  
  # define function to unfactor theta/phi or xyz just in case factors are passed as arguments
  unfactor <- function (f) as.numeric(levels(f)[as.integer(f)])
  
  # conversion factor degrees to radians (* deg-> rad; / rad -> deg)
  deg2rad <- (2 * pi) / 360
  
  # coerce to data frame
  xx <- as.data.frame(x)
  
  # either theta and phi should be present, or x, y, z
  theta <- phi <- x <- y <- z <- x2d <- y2d <- label <- FALSE
  cn <- colnames(xx)
  for (i in seq_along(cn)) {
    assign((cn[i]), TRUE)
  }
  if (!((theta && phi) || (x && y && z))) {
    stop ('input data should either contain x,y,z, or theta,phi.');
  }
  
  # x y z given - calculate theta and phi if not present
  if (x && y && z) {
    if (is.factor(xx$x)) xx$x <- unfactor(xx$x)
    if (is.factor(xx$y)) xx$y <- unfactor(xx$y)
    if (is.factor(xx$z)) xx$z <- unfactor(xx$z)
    if (!(theta && phi)) {
      r = sqrt(xx$x^2 + xx$y^2 + xx$z^2)
      xx$theta <- acos(xx$z / r) / deg2rad
      xx$theta[which(xx$x < 0)] <- -xx$theta[which(xx$x < 0)]
      xx$phi <- atan(xx$y / xx$x) / deg2rad
      xx$phi[which(abs(xx$x) < 1e-100 & abs(xx$y) < 1e-100)] <- 0
      xx$phi[which(abs(xx$x) < 1e-100 & xx$y > 0)] <- 90
      xx$phi[which(abs(xx$x)< 1e-100 & xx$y < 0)] <- -90
    }
  }
  
  # theta and phi given - calculate x y z if not present
  if (theta && phi) {
    if (is.factor(xx$theta)) xx$theta <- unfactor(xx$theta)
    if (is.factor(xx$phi)) xx$phi <- unfactor(xx$phi)
    if (!(x && y && z)) {
      xx$x <- cos(deg2rad * xx$phi) * sin(deg2rad * xx$theta)
      xx$y <- sin(deg2rad * xx$phi) * sin(deg2rad * xx$theta)
      xx$z <- cos(deg2rad * xx$theta)
    }
  }
  
  # calculate 2D coordinates if not present
  if (!(x2d && y2d)) {
    xx$x2d <- (deg2rad * xx$theta) * cos(deg2rad * xx$phi)
    xx$y2d <- (deg2rad * xx$theta) * sin(deg2rad * xx$phi)
  }

  # add labels if not present - great circle distance based on Vincenty sphere method
  # lat=phi lon=theta; 1=xx, 2=EEGlocations
  if (!label) {
    xx$label=''
    for (i in 1:nrow(xx)) {
      lat1 <- xx[i,'phi'] * deg2rad
      lon1 <- xx[i, 'theta'] * deg2rad
      lat2 <- eegr::EEGlocations$phi * deg2rad
      lon2 <- eegr::EEGlocations$theta * deg2rad
      xl <- sqrt((cos(lat2) * sin(lon1 - lon2))^2
                 + (cos(lat1) * sin(lat2)
                    - sin(lat1) * cos(lat2) * cos(lon1 - lon2))^2)
      yl <- sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)
      d <- atan2(xl, yl)
      xx[i, 'label'] <- as.character(eegr::EEGlocations[which.min(d), 'label'])
    }
  }
  
  class(xx) <- c("sensorlocs", "data.frame")
  invisible(xx)
}

#' @rdname sensorlocs
#' @param ... other arguments passed to \code{plot} function
#' @export

plot.sensorlocs <- function (x, ...) {
  
  #set plotting parameters
  op <- graphics::par(pty="s", mai=c(0,0,0,0))
  on.exit(graphics::par(op))
  
  #draw electrodes
  graphics::plot(x$x2d, x$y2d, xlim = c(-2, 2), ylim = c(-2, 2), bty = "n", axes = FALSE,
                 xlab = "", ylab = "", pch=20, ...)
  
  #draw circle
  radius <- 90 * (2 * pi) / 360
  plotrix::draw.circle(0, 0, radius)
  
  # draw nose
  graphics::segments(-0.2, radius - 0.02, 0, radius + 0.15)
  graphics::segments(0, radius + 0.15, 0.2, radius - 0.02)
  
  # draw ears
  graphics::segments(-(radius - 0.02), 0.2, -(radius + 0.06), 0.3)
  graphics::segments(-(radius + 0.06), 0.3, -(radius + 0.06), -0.3)
  graphics::segments(-(radius + 0.06), -0.3, -(radius - 0.02), -0.2)
  graphics::segments(+(radius - 0.02), 0.2, +(radius + 0.06), 0.3)
  graphics::segments(+(radius + 0.06), 0.3, +(radius + 0.06), -0.3)
  graphics::segments(+(radius + 0.06), -0.3, +(radius - 0.02), -0.2)
  
  n <- nrow(x)
  if (n > 0) {
    scale <- 1/floor(log10(n))    #adjust text size
    if (n < 10) scale <- 1
    for (i in 1:n) {
      graphics::text(x = x$x2d[i], y = x$y2d[i] + (0.1*scale), labels = x$label[i], cex = scale)
    }
  }
}
