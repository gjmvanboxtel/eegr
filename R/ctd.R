# ctd.R - define object and methods for 'ctd'
# Copyright (C) 2019  Geert van Boxtel
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
# 20190918  GvB       Setup for eegr 0.1.0 (using matrix instead of data frame)
#
#---------------------------------------------------------------------------------------------------------------------

#' Continuous Time Domain (ctd) data
#' 
#' Define the \code{ctd} object and associated methods
#' 
#' The \code{ctd} object can be used to store \strong{c}ontinuous \strong{t}ime \strong{d}omain data. 
#'  A \code{ctd} object is a numeric matrix with attributes\cr
#'  \code{fs} - the sampling frequency,\cr
#'  \code{ns} - the number of signals/channels, (columns of the matrix), and \cr
#'  \code{npts} - the number of data points/samples (rows of the matrix).\cr
#'  Methods exist to print, summarize, and plot the data.
#'
#' @param x name of a matrix that contains the data, or an object that can be coerced to a matrix,
#' consisting of \code{npts} data points as rows, and \code{ns} signals as columns.
#' @param ... other arguments
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @rdname ctd
#' @export ctd
#' 
ctd <- function (x, fs, add.t = TRUE, ...) UseMethod ("ctd", x)

#' @param fs the sampling frequency in Hertz.
#' @param add.t flag indicating whether to add a time variable to the matrix
#'
#' @return A \code{ctd} object, a numeric matrix with attributes \code{fs}, \code{ns}, and \code{npts}.
#' If the input matrix does not have column names, the columns will be named \code{S1 ... Snpts}. An additional
#' variable \code{t} containing time points will be generated if not present already and \code{add.t = TRUE} 
#'
#' @examples
#' # simulate some rather silly data that could look like EEG (2 channels)
#' # assume 2 seconds of data sampled at 100 Hz
#' C3 <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 199)+(10*rnorm(200))
#' C4 <- arima.sim(list(order = c(1,1,0), ar = 0.7), n = 199)+(10*rnorm(200))
#' eeg <- ctd(cbind(C3,C4), 100)
#' print(eeg)
#' summary(eeg)
#' 
#' @rdname ctd
#' @export

ctd.default <- function (x, fs, add.t = TRUE, ...) {
  x <- as.matrix (x)
  npts <- nrow(x)
  ns <- ncol(x)
  cn <- colnames(x)
  if (is.null(cn)) cn <- paste0("S",seq_len(ns))
  if (add.t) {
    if ('t' %in% cn) {
      ns <- ns - 1
    } else {
      x <- cbind(x, (0:(npts-1))/rep(fs, npts))
      colnames(x)[ncol(x)] <- 't'
    }
  }
  colnames(x)[1:length(cn)] <- cn
  class(x) <- c("ctd", "matrix")
  attr(x, "npts") <- npts
  attr(x, "ns") <- ns
  attr(x, "fs") <- fs
  x
}

#'
#' \code{npts} returns the \code{npts} atrribute (number of points, samples, rows) of a \code{ctd} object (or \code{NULL}).
#' \code{npts<-} sets the \code{npts} attribute
#' 
#' @param object Object, usually of class \code{ctd}
#' @param value the attribute \code{npts}, \code{ns}, or \code{fs} will be set to this value
#'
#' @rdname ctd
#' @export
npts <- function(object) attr(object, "npts")

#' @rdname ctd
#' @export
`npts<-` <- function(object, value) attr(object, "npts") <- value

#' \code{ns} returns the \code{ns} atrribute (number of signals, samples) of a \code{ctd} object (or \code{NULL}).
#' \code{ns<-} sets the \code{ns} attribute
#' 
#' @rdname ctd
#' @export
ns <- function(object) attr(object, "ns")

#' @rdname ctd
#' @export
`ns<-` <- function(object, value) attr(object, "ns") <- value

#' \code{fs} returns the \code{fs} atrribute (sampling frequency) of a \code{ctd} object (or \code{NULL}).
#' \code{fs<-} sets the \code{fs} attribute
#' 
#' @rdname ctd
#' @export
fs <- function(object) attr(object, "fs")

#' @rdname ctd
#' @export
`fs<-` <- function(object, value) attr(object, "fs") <- value

#' @rdname ctd
#' @export
print.ctd <- function (x, ...)  {
  cat("Number of data points:", npts(x), "\n")
  cat("Number of signals:", ns(x), "\n")
  cat(colnames(x))
  cat("\nSampling frequency:", fs(x), "Hz\n")
}

#' @rdname ctd
#' @importFrom stats median sd
#' @export
summary.ctd <- function (object, ...) {
  cn <- colnames(object)[1:ns(object)]
  m <- matrix(0, 5, ns(object), 
              dimnames = list(c("min","median","max","mean","sd"), cn))
  for (is in 1:ns(object)) {
    m[1, is] <- base::min(object[,is])
    m[2, is] <- stats::median(object[,is])
    m[3, is] <- base::max(object[,is])
    m[4, is] <- base::mean(object[,is])
    m[5, is] <- stats::sd(object[,is])
  }
  ans <- data.frame(m)
  colnames(ans) <- cn
  class(ans) <- c("summary.ctd", "data.frame")
  ans
}

#' @rdname ctd
#' @export
print.summary.ctd <- function (x, ...) {
  for (c in 1:ncol(x)){
    cat(paste("\nSignal:", colnames(x)[c]))
    for (r in 1:nrow(x)) {
      cat(paste("\n  ", rownames(x)[r], ":", x[r,c]))
    }
  }
}

#' @rdname ctd
#' @importFrom graphics par plot mtext axis
#' @param sensors Numeric. Sensors to plot (default \code{1:ns})
#' @param xlim Numeric. Lower and upper limits of plot X-axis (default: 0,10)
#' @param ylim Numeric. Lower and upper limits of plot Y-axis (default: -50,50)
#' @export

plot.ctd <- function (x, sensors=1:ns(x), xlim=c(0,10), ylim=c(-50,50), ...) {

  n <- sensors
  if (is.numeric(sensors)) {
    n <- colnames(x)[sensors]
  }
  xlim <- pmax(xlim, x[1,"t"])
  xlim <- pmin(xlim, x[npts(x), "t"])
  xx <- x[x[,'t'] >= xlim[1] & x[,'t'] <= xlim[2], c(n, 't')]
  ns <- length(n)
  yrange<- abs(ylim[1])+abs(ylim[2])

  op <- graphics::par(mar=c(1.1,1,0,0), las=1)
  on.exit(graphics::par(op))
  graphics::par(new = FALSE)
  for (i in 1:ns) {
    pos <- (ns-i+0.5)*yrange
    graphics::plot(xx[,'t'], (xx[ ,n[i]] - mean(xx[ ,n[i]]) + pos), type = "l", ylim = c(0, yrange * ns), axes = FALSE)
    graphics::mtext(paste(n[i], '-'), side = 2, at = pos, cex = 0.7, adj = 0.2)
    graphics::par(new = TRUE)
  }
  graphics::axis(side = 1, tick = TRUE, labels = TRUE, cex.axis = 0.7, padj = -1.9)
  graphics::par(new = FALSE)
}
