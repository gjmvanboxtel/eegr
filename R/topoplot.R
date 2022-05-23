# topoplot.R
# Draw a topographical plot of EEG data using spherical splines.
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
# 20140212    GvB           Initial setup
# 20190918    GvB           version for eegr v0.1.0
# 20210602    GvB           use Akima spherical splines
# 20220523    GvB           add fmt parameter for legend
#---------------------------------------------------------------------------------------------------------------------

#' Topographic Plot
#'
#' Draw a topographical plot of EEG data using spherical spline interpolation.
#'
#' The function \code{topoplot} draws a topographical map of EEG data \code{x}
#' recorded at sensors \code{sens}. Values in between sensor locations are
#' interpolated using spherical splines. The plot may optionally include contour
#' lines, sensor positions, sensor labels, and a legend. Scaling can be manual
#' or automatic.
#'
#' The function assumes the internal sensor coordinate system, which originates
#' in the middle of the head (0,0,0). See \code{\link{sensorlocs}}
#'
#' Multiple topoplots can also be combined on a single page or display, using
#' appropriate par(mfrow = c()), par(mfcol = c()), or split.screen() commands.
#' In such cases, you may want to draw a single legend for all graphs. To do so,
#' do not specify plotting a legend in the constituent plots, but add the
#' legends after plotting, as in the example below.
#'
#' @param x Numeric array (coerced). The data to be interpolated and plotted.
#'   The number of data points should equal the number of sensors in the
#'   \code{sl} parameter.
#' @param sl \code{\link{sensorlocs}} object, or data frame with sensor
#'   coordinates. Should at least contain the variables 'label', 'x2d' and
#'   'y2d'.
#' @param res Numeric. Resolution of the output grid, Default 200.
#' @param scale Plot using this scale. Use "auto" (default) for automatic
#'   scaling; otherwise pass a pair of c(min, max), e.g., c(-40,40).
#' @param plot Character. Specifies how the plot will look like. Levels are
#'   always plotted, "sensorlocs" plots sensor locations as dots, "legend" plots
#'   a legend at the right side of the plot, "contour" specifies contour lines
#'   with values in the plot, "labels" plots sensor labels, and "all" plots all
#'   of the above. Default: \code{c("sensorlocs", "Legend")}.
#' @param col color palette used for the map. Default: jet colors (as in
#'   Matlab), interpolated using \code{\link[grDevices]{colorRampPalette}}
#'   using 100 levels.
#' @param fmt format statement, specified as a character value, passed to sprint
#'   to format the labels of the legend. See \code{link{sprintf}} for further
#'   detail. Default: "\%+f".
#' @param xlab,ylab,main labels passed too plotting function. Default: "".
#' @param ... Additional parameters passed to the plotting function.
#'
#' @return A list with 4 components, returned invisibly: 
#' \describe{
#'   \item{x, y}{vectors of x- and y- coordinates of the output grid, each 'res'
#'   values ranging from -200 to + 200.}
#'   \item{z}{matrix of fitted z-values, dimensions 400 x 400. The fitted value
#'   at an sensor location with coordinates x2d and y2d can be found at
#'   \code{z[x2d, y2d]}.}
#'   \item{zlim}{range of z-values used in the interpolation}
#'  }
#'
#' @examples
#' data("EEGdata")
#' ncoi <- 28
#' coi <- colnames(EEGdata)[1:ncoi]
#' sensors <- getLocationsfromLabels(labels = coi)
#'
#' ## scalp distributuon of two time points, Matlab-style 'jet' colors
#' col <- grDevices::colorRampPalette(c("#00007F", "blue",
#'                                      "#007FFF", "cyan",
#'                                      "#7FFF7F", "yellow",
#'                                      "#FF7F00", "red",
#'                                      "#7F0000"))(100)
#' op <- par(mfrow = c(1, 2))
#' topoplot (EEGdata[400, coi], sensors, scale = c(-10, 10), main="Sample #400",
#'           col = col, plot = "sensorlocs")
#' topoplot (EEGdata[450, coi], sensors, scale = c(-10, 10), main="Sample #450",
#'           col = col, plot = "sensorlocs")
#' par(op)
#' title("Example displaying two scalp maps")
#' usr <- graphics::par("usr")
#' plotrix::color.legend(xl = usr[2] + 0.15, yb = usr[3] + 0.35,
#'       xr = usr[2] + 0.35, yt = usr[4] - 0.35,
#'       rect.col = col, gradient="y", align="rb",
#'       legend = sprintf("%+2.0f", seq(-10, 10, length = 5)))
#'
#' ## this color scheme might be more appropriate for scalp plots
#' topoplot(EEGdata[400, coi], sensors, scale = c(-20, 20),
#'          col = grDevices::colorRampPalette(c("blue", "white", "red"))(50))
#'
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#'
#' @seealso \code{\link[grDevices]{colorRampPalette}}
#'
#' @export

topoplot <- function (x, sl, res = 200, scale = "auto",
                      plot = c("sensorlocs", "legend"),
                      col = grDevices::colorRampPalette(c("#00007F", "blue",
                                                          "#007FFF", "cyan",
                                                          "#7FFF7F", "yellow",
                                                          "#FF7F00", "red",
                                                          "#7F0000"))(100),
                      fmt = "%+f", xlab = "", ylab = "", main = "", ...) {
  
  # parameter checking
  snames <- names(x)
  x <- as.numeric(x)
  if (anyNA(x)) {
    stop("'x' must not contain any missing values")
  }
  sl <- sensorlocs(sl)
  if (!("x2d" %in% colnames(sl)) || !("y2d" %in% colnames(sl))) {
    stop ('sensorlocs data frame must contain x2d and y2d variables')
  }
  if(!is.null(snames) && !identical(as.character(sl$label), snames)) {
    stop('names of variables in x must match labels in sensorlocs data frame')
  }
  if (is.na(res) || is.null(res) || !is.numeric(res)) {
    res <- 200
  }
  if (is.character(scale)) {
    if (scale != "auto") {
      stop('scale must be "auto" (or a numeric vector of length 2)')
    }
  } else if (is.numeric(scale)) {
    if (length(scale) != 2) {
      stop('scale must be a numeric vector of length 2 (or "auto")')
    }
  } else {
    stop('scale must be "auto" or a numeric vector of length 2')
  }
  if (!is.character(fmt) || substr(fmt, 1, 1) != "%") {
    stop("fmt must be a character string beginning with '%'")
  }
  plot <- match.arg(plot, 
                    c("sensorlocs", "legend", "contour", "labels", "all"),
                    several.ok = TRUE)
  # end of parameter checking

  # spherical spline interolation using the 'akima' package
  ip <- akima::interp(sl$x2d, sl$y2d, x, linear = FALSE, extrap = TRUE,
                      xo = seq(-2, 2, length.out = res),
                      yo = seq(-2, 2, length.out = res))
  
  # interpolation produces a square grid. Make it a circle - Sensor Locations
  # range from from -pi/2 to +pi/2 radians
  max_xy <- max(sl$x2d, sl$y2d, pi / 2)
  r <- sqrt(outer(ip$x^2, ip$y^2, "+"))
  ip$z[r > max_xy] <- NA 
  
  # now define the scaling limits for the data
  # autoscaling: compute minimum and maximum of the prediction matrix
  if ("auto" %in% scale) {
    zlim <- c(floor(min(ip$z, na.rm = TRUE)), ceiling(max(ip$z, na.rm = TRUE)))
  } else {
    zlim = scale
  }
  # clip values outside of limits
  ip$z[which(ip$z < zlim[1])] <- zlim[1]
  ip$z[which(ip$z > zlim[2])] <- zlim[2]
  
  #define plotting region (square, no margins)
  op <- graphics::par(pty = "s", mar = c(2.5, 2.5, 2.5, 2.5))
  on.exit(graphics::par(op))
  
  # square plot using 'image', no axes
  graphics::image(ip$x, ip$y, ip$z, axes = FALSE, col = col,
                  xlim = c(-2, 2), ylim = c(-2, 2), zlim = zlim, 
                  xlab = xlab, ylab = ylab, main = main, ...)

  # contour requested? compute range and define range/5 contour levels
  if ("contour" %in% plot || "all" %in% plot) {
    range <- sum(abs(zlim))
    graphics::contour(ip$x, ip$y, ip$z, add = TRUE, axes = FALSE, nlevels = 5)
  }
  
  # sensorlocs requested?
  if ("sensorlocs" %in% plot || "all" %in% plot) {
    graphics::points(sl$x2d, sl$y2d, xlim = c(-2, 2), ylim = c(-2, 2),
                     bty = "n", pch = 20)
  }
  
  # sensor labels requested?
  if ("labels" %in% plot || "all" %in% plot) {
    for (i in 1:nrow(sl)) {
      graphics::text(x = sl$x2d[i], y = (sl$y2d[i] + 0.1), labels = sl$label[i])
    }
  }
  
  #draw circle
  deg2rad <- (2 * pi) / 360
  radius <- 90 * deg2rad
  xy <- plotrix::draw.circle(0.0, 0, radius, lwd = 2)

  # draw nose
  graphics::segments(-0.2, radius - 0.02, 0.0, radius + 0.15, lwd = 2)
  graphics::segments(0.0, radius + 0.15, 0.2, radius - 0.02, lwd = 2)
  
  # draw ears
  graphics::segments(-(radius - 0.02), +0.2, -(radius + 0.06), +0.3, lwd = 2)
  graphics::segments(-(radius + 0.06), +0.3, -(radius + 0.06), -0.3, lwd = 2)
  graphics::segments(-(radius + 0.06), -0.3, -(radius - 0.02), -0.2, lwd = 2)
  graphics::segments(+(radius - 0.02), +0.2, +(radius + 0.06), +0.3, lwd = 2)
  graphics::segments(+(radius + 0.06), +0.3, +(radius + 0.06), -0.3, lwd = 2)
  graphics::segments(+(radius + 0.06), -0.3, +(radius - 0.02), -0.2, lwd = 2)
  
  #legend requested? Set it in user coordinates at the left of the plot
  #electrode labels requested?
  if ("legend" %in% plot || "all" %in% plot) {
    usr <- graphics::par("usr")
    plotrix::color.legend(xl = usr[2], yb = usr[3] + 0.2,
                          xr = usr[2] + 0.2, yt = usr[4] - 0.2,
                          rect.col = col, gradient="y", align="rb",
                          legend = sprintf(fmt, seq(zlim[1], zlim[2],
                                                         length = 5))
    )
    # # tick marks
    graphics::par(mar = c(0, 0, 0, 0))
    steps <- (usr[4] - 0.2) / 2
    for (i in 0:4) {
      graphics::segments(usr[2], (usr[3] + 0.2 + (i * steps)),
                         (usr[2] + 0.2), usr[3] + 0.2 + (i * steps),
                         col = "black", lwd = 2)
    }
  }
  
  #ready
  invisible(list(x = ip$x, y = ip$y, z = ip$z, zlim = zlim))
}

