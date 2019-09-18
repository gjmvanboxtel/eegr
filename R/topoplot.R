# topoplot.R
# Draw a topographical plot of EEG data using spherical splines.
# Copyright (C) 2014, 2019  Geert van Boxtel,
# Tilburg University, G.J.M.vBoxtel@tilburguniversity.edu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# See also: http://www.gnu.org/licenses/gpl-2.0.txt
#
# Version history:
# 20140212    GvB           Initial setup
# 20190918    GvB           version for eegr v0.1.0
#---------------------------------------------------------------------------------------------------------------------

#' Topographic Plot
#' 
#' Draw a topographical plot of EEG data using spherical spline interpolation
#' 
#' The function \code{topoplot} draws a topographical map of EEG data \code{x} recorded at sensors \code{sens}.
#' Values in between sensor locations are interpolated using spherical splines, as proposed by Perrin et al. (1989,
#' see also the Corrigendum). The plot may optionally include contour lines, sensor positions, sensor labels, and a
#' legend. Scaling can be manual or automatic.
#' 
#' The function assumes the internal sensor coordinate system, which originates in the middle of the head (0,0,0). See
#' \code{\link{sensorlocs}}
#' 
#' Multiple topoplots can also be combined on a single page or display, using appropriate par(mfrow = c()), par(mfcol = c()),
#' or split.screen() commands. In such cases, you may want to draw a single legend for all graphs. To do so, do not specify plotting
#' a legend in the constituent plots, but add the legands after plotting, as in the example below.
#' 
#' @param x Numeric array (coerced). The data to be interpolated and plotted. The number of data points should
#' equal the number of sensors in the \code{sens} parameter.
#' @param sl \code{\link{sensorlocs}} object, or data frame with sensor coordinates. Should at least contain the variables
#' 'label', 'x2d' and 'y2d'.
#' @param res Numeric. Resolution of the output grid, Default 200.
#' @param scale Plot using this scale. Use "auto" (default) for automatic scaling; otherwise pass a pair of c(min, max),
#'  e.g., c(-40,40).
#' @param plot Character. Specifies how the plot will look like. Levels are always plotted, "sensorlocs" (default) plots sensor
#' locations as dots, "legend" plots a legend at the right side of the plot, "contour" specifies contour lines with values in the plot,
#' "labels" plots sensor labels, and "all" plots all of the above.
#' @param ... Additional parameters passed to the 'image' function, e.g. main=
#' 
#' @examples 
#' \dontrun{
#' data <- ReadEDF ("yourfile")
#' sensors <- getLocationsfromLabels(labels = data$head$label)
#' opar <- par(xpd=NA, mfrow=c(1,2));
#' ret <- topoplot (data$ctd[200,], sensors, scale = c(-10, 10), plot = c("contour","sensorlocs"), main="Plot1");
#' ret <- topoplot (data$ctd[400,], sensors, scale = c(-10, 10), plot = c("contour","sensorlocs"), main="Plot2");
#' usr <- par("usr")
#' plotrix::color.legend(xl = usr[2], yb = (usr[3] + 23), xr = (usr[2] + 23), yt = (usr[4] - 23)
#'                       rect.col = rev(rainbow(250,start=0,end=.7)),
#'                       gradient = "y", align = "rb",
#'                       legend = sprintf("%+2.0f", seq(-10, 10, length=5)))
#' par(opar)
#' # Note: use appropriate scaling min and max scaling limits
#' }
#' 
#' @return 
#' list with 4 components:
#' \describe{
#'   \item{x,y}{vectors of x- and y- coordinates of the output grid, each 'res' values ranging from -200 to + 200.}
#'   \item{z}{matrix of fitted z-values, dimensions 400 x 400. The fitted value at an sensor location with coordinates
#'   x2d and y2d can be found at z[(x2d * 360/pi), (y2d * 360/pi)].}
#'   \item{zlim}{range of z-values used in the interpolation}
#' }
#' 
#' @note Spline interpolation is implemented through the \code{sspline} package, which implements the methods proposed by Wahba (1981),
#' on which the work of Perrin et al. (1989) is also based. The \code{sspline} package works with latitudes and longitudes in degrees
#' instead of sensor locations as defined internally in the \code{eegr} package. This i s normally hidden from the user calling this
#' function, and only important for users who want to work with the underlying code.
#' 
#' @references Perrin F., Pernier J., Bertrand O., and Echallier JF. (1989). Spherical splines for scalp potential and current
#' density mapping. \emph{Electroencephalogr. Clin. Neurophysiol.}, \emph{72(2)}, 184-187.
#' \url{https://doi.org/10.1016/0013-4694(89)90180-6}
#' @references Perrin, F., Pernier, J., Bertrand, O. and Echallier, J.F. (1990) Corrigenda: EEG 02274.
#' \emph{Electroenceph. Clin. Neurophysiol.}, \emph{76}, 565. \url{https://doi.org/10.1016/0013-4694(90)90009-9}
#' @references Wahba, G., (1981). Spline interpolation and smoothing on the sphere. \emph{SIAM J. Sci. Stat. Comput.},
#' \emph{2(1)}, 5-16. \url{https://doi.org/10.1137/0902002}
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.

#' @importFrom graphics par image segments text image contour
#' @importFrom plotrix draw.circle
#' @importFrom sspline smooth.sspline predict.smooth.sspline
#' @export


topoplot <- function (x, sl, res = 200, scale = "auto", plot = c("sensorlocs", "legend"), ...) {
  
  x <- as.numeric(x)

  if (!"theta" %in% colnames(sl)) stop ('sensorlocs data frame must contain theta variable')
  if (!"phi" %in% colnames(sl)) stop ('sensorlocs data frame must contain phi variable')
  sl <- sensorlocs(sl)

  lat <- sl$theta
  lon <- sl$phi
  if (!(length(x) == length(lat) && length(x) == length(lon))) stop ('The number of data points must equal the number of sensors')
  if(is.na(res) || is.null(res) || !is.numeric(res)) res <- 200
  
  # define interpolated latitude and longitude grid
  # the output grid is of size res * res, varing from -200 to +200
  # (a little more than 180 degrees; the rest will be set to NA later)
  iplat <- seq(-200, 200, len = res)
  iplon <- seq(-200, 200, len = res)
  
  # spherical spline interolation using the 'sspline' package
  spl <- sspline::smooth.sspline(lat = lat, lon = rev(lon), x)
  pred <- sspline::predict.smooth.sspline(spl, lat = iplat, lon = iplon, grid = TRUE)
  
  # interpolation produces a square grid. Make it a circle
  # Sensor Locations range from from -pi/2 to +pi/2 radians, the interpolated data from to -180 to +180 degrees
  # So the conversion ratio is 360/pi
  rad2deg <- 360 / pi
  max_xy <- max((sl$x2d * rad2deg), (sl$y2d * rad2deg), 180)
  r <- sqrt(outer(iplon^2, iplat^2, "+"))
  pred[r > max_xy] <- NA
  
  # now define the scaling limits for the data
  # autoscaling: compute minimum and maximum of the prediction matrix
  if ("auto" %in% scale) {
    zlim <- c(floor(min(pred, na.rm=T)), ceiling(max(pred, na.rm=T)))
  } else {;
    zlim = scale
  }
  # clip values outside of limits
  pred[which(pred < zlim[1])] <- zlim[1]
  pred[which(pred > zlim[2])] <- zlim[2]
  
  #Now start plotting...
  plot <- match.arg(plot, several.ok = TRUE)
  
  #define plotting region (square, no margins)
  op <- graphics::par(pty="s", mar = c(2, 2, 2, 2))
  
  # square plot using 'image', no axes, rainbow colors
  graphics::image(iplat, iplon, pred, axes = FALSE, xlim=c(-200, 200), ylim = c(-200, 200), zlim=zlim, 
                  xlab = "", ylab = "", col = rev(rainbow(250, start=0, end=.7)), ...)

  #contour requested? compute range and define range/5 contour levels
  if ("contour" %in% plot || "all" %in% plot) {
    range <- sum(abs(zlim))
    graphics::contour(iplat, iplon, pred, add = TRUE, axes = FALSE, nlevels = 5)
  }
  
  #electrodes requested?
  if ("electrodes" %in% plot || "all" %in% plot) {
    graphics::points((sl$x2d * rad2deg), (sl$y2d * rad2deg), xlim = c(-200, 200), ylim=c(-200, 200), bty = "n", pch = 20)
  }
  
  #electrode labels requested?
  if ("labels" %in% plot || "all" %in% plot) {
    for (i in 1:nrow(sl)) {
      text(x = (sl$x2d[i] * rad2deg), y = ((sl$y2d[i] * rad2deg) + 12), labels = sl$label[i])
    }
  }
  
  #draw circle
  radius <- 180
  xy <- plotrix::draw.circle(0.0, 0, radius, lwd = 2)
  
  # draw nose
  graphics::segments(-22, (radius - 2), 0, (radius + 17), lwd = 2)
  graphics::segments(0, (radius + 17), 22, (radius - 2), lwd=2)
  
  # draw ears
  graphics::segments(-(radius - 2), 23, -(radius + 7), 34, lwd = 2)
  graphics::segments(-(radius + 7), 34, -(radius + 7), -34, lwd = 2)
  graphics::segments(-(radius + 7), -34, -(radius - 2), -23, lwd = 2)
  graphics::segments(+(radius - 2), 23, +(radius + 7), 34, lwd = 2)
  graphics::segments(+(radius + 7), 34, +(radius + 7), -34, lwd = 2)
  graphics::segments(+(radius + 7), -34, +(radius - 2), -23, lwd = 2)

  #legend requested? Set it in user coordinates at the left of the plot
  #electrode labels requested?
  if ("legend" %in% plot || "all" %in% plot) {
    usr <- par("usr")
    plotrix::color.legend(xl = usr[2], yb = (usr[3] + 23), xr = (usr[2] + 23), yt = (usr[4] - 23),
                          rect.col = rev(rainbow(250, start = 0, end = .7)),
                          gradient = "y", align = "rb",
                          legend = sprintf("%+2.0f", seq(zlim[1], zlim[2], length = 5))
                          )
    # # tick marks
    par (mar = c(0, 0, 0, 0))
    steps <- (usr[4] - 23) / 2
    for (i in 0:4) {
      graphics::segments(usr[2], (usr[3] + 23 + (i * steps)), (usr[2] + 23), (usr[3] + 23 + (i * steps)), lwd = 2)
    }
  }
  
  #reset plot parameters
  par(op)
  
  #ready
  invisible(list(x = iplat, y = iplon, z = pred, zlim = zlim))
}
