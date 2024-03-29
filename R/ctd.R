# ctd.R - define object and methods for 'ctd'
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
# 20190918  GvB       Setup for eegr 0.1.0 (using matrix instead of data frame)
# 20201224  GvB       v0.3-0: documentation; make all attributed explicitely numeric
#                     (bug 20200416)
# 20211215  GvB       v0.3-2: added .ctd extension to fs(), ns(), npts()
#                     remove the assignment forms fs<- ns<- npts<-
# 20220122  GvB       added ref parameterg
# 20220205  GvB       bugfix in plot.ctd xlim calculation
#---------------------------------------------------------------------------------------------------------------------

#' Continuous Time Domain Data
#' 
#' Define the \code{ctd} object and associated methods.
#' 
#' The \code{ctd} object can be used to store \strong{c}ontinuous \strong{t}ime
#' \strong{d}omain data.
#'  A \code{ctd} object is a numeric matrix with attributes:
#'  \describe{
#'    \item{fs}{sampling frequency}
#'    \item{ns}{number of signals/channels, (columns of the matrix)}
#'    \item{npts}{the number of data points/samples (rows of the matrix)}
#'    \item{ref}{names or numbers of reference sensors}
#' }
#' Methods exist to print, summarize, and plot the data.
#'
#' @param x name of a matrix that contains the data, or an object that can be
#'   coerced to a matrix, consisting of \code{npts} data points as rows, and
#'   \code{ns} signals as columns.
#' @param fs the sampling frequency in Hertz.
#' @param ref names or numbers of reference sensors, "avg" or "average" for
#'   average reference. Default: NA
#' @param add.t flag indicating whether to add a time variable to the matrix
#' @param ... other arguments
#'
#' @return A \code{ctd} object, a numeric matrix with attributes \code{fs},
#'   \code{ns}, and \code{npts}. If the input matrix does not have column names,
#'   the columns will be named \code{S1 ... Snpts}. An additional variable
#'   \code{t} containing time points will be generated if not present already
#'   and \code{add.t = TRUE}
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
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @rdname ctd
#' @export

ctd <- function (x, fs, ref = NA, add.t = TRUE, ...) UseMethod ("ctd", x)

#' @rdname ctd
#' @export

ctd.default <- function (x, fs, ref = NA, add.t = TRUE, ...) {
  x <- as.matrix (x)
  npts <- nrow(x)
  ns <- ncol(x)
  cn <- colnames(x)
  if (is.null(cn)) {
    cn <- paste0("S",seq_len(ns))
  }
  if (!chks(x, ref)) {
    stop("invalid reference sensors specified")
  }
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
  attr(x, "npts") <- as.numeric(npts)
  attr(x, "ns") <- as.numeric(ns)
  attr(x, "fs") <- as.numeric(fs)
  attr(x, "ref") <- ref
  x
}

#'
#' \code{npts} returns the \code{npts} attribute (number of points, samples,
#' rows) of a \code{ctd} object (or \code{NULL}).
#'
#' @param object Object

#' @rdname ctd
#' @export
npts <- function (object) UseMethod ("npts")

#' @rdname ctd
#' @export
npts.ctd <- function(object) as.numeric(attr(object, "npts"))

#' \code{ns} returns the \code{ns} atrribute (number of signals, samples) of a
#' \code{ctd} object (or \code{NULL}).
#' 
#' @rdname ctd
#' @export
ns <- function (object) UseMethod ("ns")

#' @rdname ctd
#' @export
ns.ctd <- function(object) as.numeric(attr(object, "ns"))

#' \code{fs} returns the \code{fs} attribute (sampling frequency) of a
#' \code{ctd} object (or \code{NULL}).
#' 
#' @rdname ctd
#' @export
fs <- function (object) UseMethod ("fs")

#' @rdname ctd
#' @export
fs.ctd <- function(object) as.numeric(attr(object, "fs"))

#' \code{ref} returns the \code{ref} attribute (reference) of a
#' \code{ctd} object (or \code{NULL}).
#' 
#' @rdname ctd
#' @export
ref <- function (object) UseMethod ("fs")

#' @rdname ctd
#' @export
ref.ctd <- function(object) attr(object, "ref")

#' @rdname ctd
#' @export
print.ctd <- function (x, ...)  {
  cat("Number of data points:", npts(x), "\n")
  cat("Number of signals:", ns(x), "\n")
  cat(colnames(x))
  cat("\nSampling frequency:", fs(x), "Hz\n")
}

#' @rdname ctd
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
#' @param sensors Numeric. Sensors to plot (default \code{1:ns})
#' @param xlim Numeric. Lower and upper limits of plot X-axis (default: 0,10)
#' @param ylim Numeric. Lower and upper limits of plot Y-axis (default: -50,50)
#' @export

plot.ctd <- function (x, sensors=1:ns(x), xlim=c(0,10), ylim=c(-50,50), ...) {

  if(!"ctd" %in% class(x)) stop("not a ctd object")
  
  n <- sensors
  if (is.numeric(sensors)) {
    n <- colnames(x)[sensors]
  }
  lim <- pmax(xlim, x[1,"t"])
  lim <- pmin(xlim, x[npts(x), "t"])
  if (diff(lim) <= 0) lim <- xlim - xlim[1]
  xx <- x[x[,'t'] >= lim[1] & x[,'t'] <= lim[2], c(n, 't')]
  ns <- length(n)
  yrange<- abs(ylim[1])+abs(ylim[2])

  op <- graphics::par(mar=c(1.1,1,0,0), las=1)
  on.exit(graphics::par(op))
  graphics::par(new = FALSE)
  for (i in 1:ns) {
    pos <- (ns-i+0.5)*yrange
    graphics::plot(xx[,'t'], (xx[ ,n[i]] - mean(xx[ ,n[i]]) + pos),
                   type = "l", ylim = c(0, yrange * ns), axes = FALSE)
    graphics::mtext(paste(n[i], '-'), side = 2, at = pos, cex = 0.7, adj = 0.2)
    graphics::par(new = TRUE)
  }
  graphics::axis(side = 1, tick = TRUE, labels = TRUE, cex.axis = 0.7, padj = -1.9)
  graphics::par(new = FALSE)
}
