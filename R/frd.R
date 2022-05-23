# frd.R - define object and methods for 'frd'
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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Version history
# 20211222  GvB       Setup for eegr 0.3-1
#---------------------------------------------------------------------------------------------------------------------

#' Frequency Domain Data
#' 
#' Define the \code{frd} object and associated methods.
#' 
#' The \code{frd} object can be used to store \strong{fr}equency \strong{d}omain
#' data.
#'  A \code{frd} object is a numeric matrix with at least a column named \code{f} containing
#'  the frequencies, one or more data columns (usually named), and the following attributes:
#'  \describe{
#'    \item{fs}{sampling frequency}
#'    \item{fu}{frequency units}
#'    \item{ns}{number of signals, (columns of the matrix)}
#'    \item{npts}{the number of data points/samples (rows of the matrix)}
#'    \item{type}{type of data: spectrum, cross-spectrum, phase, coherence}
#'    \item{yscale}{Y-scale of the data: linear, log, dB}
#' }
#' Methods exist to print, summarize, and plot the data.
#'
#' @param x name of a matrix that contains the data, or an object that can be
#'   coerced to a matrix, consisting of \code{npts} data points as rows, and
#'   \code{ns} signals as columns.
#' @param f vector of frequencies
#' @param fu frequency units, \code{Hz} (default), \code{rad/s}, or
#'   \code{normalized}
#' @param fs sampling frequency (default: 1)
#' @param type type of data; one of \code{'spectrum'}, \code{'cross-spectrum'},
#'   \code{'phase'}, \code{'coherence'}
#' @param yscale scale of the data, one of \code{'linear'}, \code{'log'},
#'   \code{'dB'}
#' 
#' @examples 
#' # coming up
#' f <- 1
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @rdname frd
#' @export

frd <- function (x, ...) UseMethod ("frd", x)

#' @return A \code{frd} object, a numeric matrix with attributes \code{fs},
#'   \code{ns}, \code{npts}, \code{type}, and \code{yscale}. If the input matrix
#'   does not have column names, the columns will be named \code{S1 ... Snpts}.
#'
#' @examples
#' # coming up
#' 
#' @rdname frd
#' @export

frd.default <- function (x, f, fu = c("Hz", "rad/s", "normalized"), fs = 1,
                         type = c("spectrum", "cross-spectrum", "phase", "coherence"),
                         yscale = c("linear", "log", "dB"), ...) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  if (is.vector(f)) {
    f <- matrix(f, ncol = 1)
  }
  if (nrow(x) != nrow(f)) {
    stop("Number of rows in x and f do not match")
  }
  fu <- match.arg(fu)
  type <- match.arg(type)
  yscale <- match.arg(yscale)
  npts <- nrow(x)
  ns <- ncol(x)
  cn <- colnames(x)
  if (is.null(cn)) {
    cn <- paste0("S", seq_len(ns))
  }

  y <- cbind(x, f)
  colnames(y) <- c(cn, "f")
  class(y) <- c("frd", "matrix")
  attr(y, "npts") <- as.numeric(npts)
  attr(y, "ns") <- as.numeric(ns)
  attr(y, "fs") <- as.numeric(fs)
  attr(y, "fu") <- fu
  attr(y, "type") <- type
  attr(y, "yscale") <-yscale
  y
}

#'
#' \code{npts} returns the \code{npts} attribute (number of points, samples,
#' rows) of a \code{frd} object (or \code{NULL}).
#' 
#' @param object object
#'
#' @rdname frd
#' @export

npts.frd <- function(object) as.numeric(attr(object, "npts"))

#' \code{ns} returns the \code{ns} attribute (number of signals, samples) of a
#' \code{frd} object (or \code{NULL}).
#' 
#' @rdname frd
#' @export

ns.frd <- function(object) as.numeric(attr(object, "ns"))

#' \code{fs} returns the \code{fs} attribute (sampling frequency) of a
#' \code{frd} object (or \code{NULL}).
#' 
#' @rdname frd
#' @export

fs.frd <- function(object) as.numeric(attr(object, "fs"))

#' \code{fu} returns the \code{fu} attribute (frequency units) of a \code{frd}
#' object (or \code{NULL}).

#' @rdname frd
#' @export
fu <- function (object) UseMethod ("fu")

#' @rdname frd
#' @export

fu.frd <- function(object) attr(object, "fu")

#' \code{type} returns the \code{type} attribute (spectrum, cross-spectrum,
#' phase, coherence) of a \code{frd} object (or \code{NULL}).

#' @rdname frd
#' @export
type <- function (object) UseMethod ("type")

#' @rdname frd
#' @export

type.frd <- function(object) attr(object, "type")

#' \code{yscale} returns the \code{yscale} attribute (linear, log, dB) of a
#' \code{frd} object (or \code{NULL}).
#'
#' @param x object

#' @rdname frd
#' @export
yscale <- function (object) UseMethod ("yscale")

#' @rdname frd
#' @export

yscale.frd <- function(object) attr(object, "yscale")

#' @rdname frd
#' @export

print.frd <- function (x, ...)  {
  cat("Number of frequency points:", npts(x), "\n")
  cat("Number of signals:", npts(x), "\n")
  cat(colnames(x))
  cat("\nSampling frequency:", fs(x), fu(x), "\n")
  cat("Frequency resolution:", fs(x) / npts(x), fu(x) )
  cat("Type of data:", type(x), "\n")
  cat("Y-Scale:", yscale(x), "\n")
}


#' @rdname frd
#' @param sensors Numeric. Sensors to plot (default \code{1:ns})
#' @param xlim Numeric. Lower and upper limits of plot X-axis (default: 0,10)
#' @param ylim Numeric. Lower and upper limits of plot Y-axis (default: -50,50)
#' @param ... additional arguments passed to functions

#' @export

plot.frd <- function (x, sensors = setdiff(colnames(x), "f"),
                      yscale = c("linear", "log", "dB"),
                      xlim = c(0, fs(x) / 2), ylim = NULL, ...) {

  if(!("frd" %in% class(x))) stop("not an frd object")
  
  if (is.null(sensors)) {
    stop("not an frd object")
  } else {
    n <- sensors
  }
  ns <- length(n)
  
  yscx <- yscale(x)
  if (yscx == "linear") {
    if (yscale == "log") {
      x[, n] <- log10(x[, n])
    } else if (yscale == "dB") {
      x[, n] <- 10 * log10(x[, n])
    }
  } else if (yscx == "log") {
    if (yscale == "linear") {
      x[, n] <- 10^x[, n]
    } else if (yscale == "dB") {
      x[, n] <- 10 * x[, n]
    }
  } else if (yscx == "dB") {
    if (yscale == "linear") {
      x[, n] <- 10^(x[, n] / n) 
    } else if (yscale == "log") {
      x[, n] <- x[, n] / 10
    }
  }
  xlim <- pmax(xlim, x[1, "f"])
  xlim <- pmin(xlim, x[npts(x), "f"])
  xx <- x[x[,'f'] >= xlim[1] & x[,'f'] <= xlim[2], n]
  freq <- x[, "f"]
  if (is.null(ylim)) {
    ylim <- range(xx)
  }

  graphics::matplot(freq, xx, type = "l", xlim = xlim, ylim = ylim,
                    xlab = paste0("Frequency (", fu(x), ")"), 
                    ylab = paste0(type(x), " (", yscale, ")"), ...)
}
