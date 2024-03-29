# EEGlocations.R
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Version history
# 20190918  GvB       Setup for eegr 0.1.0
# 20201224  GvB       Documentation adapted for v0.3-0
#---------------------------------------------------------------------------------------------------------------------

#' EEGlocations
#'
#' Data frame of class \code{\link{sensorlocs}} containing 346 sensor locations
#' corresponding to the 10-5 electrode system for EEG.
#' 
#' A \code{\link{sensorlocs}} object is a \code{\link{data.frame}} containing
#' sensor locations that can be used in further analyses or display. In order to
#' match R's plotting functions, the coordinate system used in \code{eegr}
#' models the head as a perfect sphere with a radius of 1 and the origin in the
#' middle of the head (0,0,0).
#' 
#' Within this system, Cartesian coordinates are defined along the following
#' axes:
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
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @examples
#' data(EEGlocations)
#' 
#' @format A data frame with 356 rows and 8 columns:
#' \describe{
#'   \item{label}{(\code{factor}) electrode label according to the 10-5 system,
#'   e.g., "AF1", "P3", etc.}
#'   \item{x}{(\code{numeric}) Cartesian x-axis (left-to-right) value}
#'   \item{y}{(\code{numeric}) Cartesian y-axis (back-to-front) value}
#'   \item{z}{(\code{numeric}) Cartesian z-axis (bottom-to-top) value}
#'   \item{x2d}{(\code{numeric}) x-axis value of orthogonal projection onto 2D
#'   plane}
#'   \item{y2d}{(\code{numeric}) y-axis value of orthogonal projection onto 2D
#'   plane}
#'   \item{theta}{(\code{numeric}) angle between 3D z-axis and x/y plane}
#'   \item{phi}{(\code{numeric}) angle between x- and y-axes}
#' }
"EEGlocations"