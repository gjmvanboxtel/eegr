# eeginterp.R
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
# 20220126  GvB           Initial setup v0.3-2
# 20220215  GvB           Bugfix bad sensor index, catch errors in
#                         spheric_spline()
#---------------------------------------------------------------------------------------------------------------------

#' eeginterp
#' 
#' Spherical spline interpolation of EEG data
#' 
#' The data at bad sensors is replaced by interpolated data from remaining
#' sensors using spherical spline interpolation (Perrin et al., 1989).
#' 
#' @param x input time series, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param sl sensor locations of \code{x}, specified as a data frame according
#'   to \code{\link{sensorlocs}} format. Alternatively, a matrix with named
#'   columns x, y, and z representing the sensor locations. Default: obtained
#'   from the sensor names given by the column names of \code{x}.
#' @param bad sensor locations of the bad sensors (those to interpolate),
#'   specified in one of the following forms:
#'   \itemize{
#'     \item a data frame in \code{\link{sensorlocs}} format
#'     \item a matrix with named columns x, y, and z representing the sensor
#'     locations
#'     \item a character vector containing sensor names
#'     \item a numeric vector representing columns of x
#'   }
#' @param method method used for interpolation, specified as a character string.
#'   Currently, only "spherical" is supported.
#' 
#' @return interpolated bad sensors, returned as an object of the same class and
#'   dimensions as the input
#'   
#' @examples 
#' \dontrun{
#' x <- EEGdata[, 1:28]
#' y <- eeginterp(x, bad = "T7, T8")
#' plot(ctd(x, fs = fs(EEGdata)))
#' plot(ctd(y, fs = fs(EEGdata)))
#' }
#' 
#' 
#' @author Matlab code for EEGlab, Copyright (C) Arnaud Delorme, CERCO, 2006,
#'   \email{arno@@salk.edu}. Ported to R and adapted by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}
#'   
#' @references Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1989).
#'   Spherical splines for scalp potential and current density mapping.
#'   Electroencephalography and clinical neurophysiology, 72(2), 184-187.
#'
#' @export

eeginterp <- function (x, sl = getLocationsfromLabels(colnames(x)),
                       bad,
                       method = "spherical") {
  
  # check input arguments : x
  if (!is.matrix(x)) {
    stop("x must be a matrix or an object of class ctd")
  }
  atr <- attributes(x)
  
  # method
  method <- match.arg(method)
  
  # sl: sensorlocs
  if ("data.frame" %in% class(sl) || "matrix" %in% class(sl)) {
    if (nrow(sl) != ncol(x)) {
      stop("invalid specification of sensor locations")
    }
    if (!all(c("x", "y", "z") %in% colnames(sl))) {
      stop("sl must contain columns named x, y, and z")
    }
  } else {
    stop("sensor locations must be specified as a data frame or a matrix")
  }
  xsens <- sl[, "x"]
  ysens <- sl[, "y"]
  zsens <- sl[, "z"]
  rad <- sqrt(xsens^2 + ysens^2 + zsens^2)
  xsens <- xsens / rad
  ysens <- ysens / rad
  zsens <- zsens / rad

  # bad sensors  
  if ("data.frame" %in% class(bad) || "matrix" %in% class(bad)) {
    # specified as sensorlocs or matrix
    if (!all(c("x", "y", "z") %in% colnames(bad))) {
      stop("if 'bad' is a matrix or data frame, it must contain columns named x, y, and z")
    }
    xbad <- bad[, "x"]
    ybad <- bad[, "y"]
    zbad <- bad[, "z"]
  } else if (is.character(bad)) {
    # specified as character vector
    if (length(bad) == 1) {
      # perhaps specified as single string: split it
      bad <- unlist(strsplit(bad, "[ ,;]"))
      bad <- bad[bad != ""]
    }
    if ("label" %in% colnames(sl) && all(bad %in% sl$label)) {
      # try to find labels in sl
      rows <- sl[sl$label %in% bad, ]
    } else if (all(bad %in% colnames(x)) && all(bad %in% EEGlocations$label)) {
      # or take them from EEGlocations
      rows <- EEGlocations[EEGlocations$label %in% bad, ]
    } else {
      stop("invalid character specification of bad sensors")
    }
    xbad <- rows[, "x"]
    ybad <- rows[, "y"]
    zbad <- rows[, "z"]
  } else if (is.numeric(bad)) {
    #specified as numbers
    if (all(bad < nrow(x))) {
      rows <- sl[bad, ]
      xbad <- rows[, "x"]
      ybad <- rows[, "y"]
      zbad <- rows[, "z"]
    } else {
      stop("invalid numeric specification of bad sensors")
    }
  }
  rad <- sqrt(xbad^2 + ybad^2 + zbad^2)
  xbad <- xbad / rad
  ybad <- ybad / rad
  zbad <- zbad / rad
  
  # use index numbers of bad sensors
  #badidx <- which(xsens %in% xbad & ysens %in% ybad & zsens %in% zbad)
  badidx <- which(colnames(x) %in% bad)
  if (length(badidx) <= 0) {
    return(x)
  }
  
  # perform the interpolation
  # if interpolation failed, loop until it does
  if (method == "spherical") {
    while (length(badidx) > 0) {
      err <- try(ss <- spheric_spline(
        xsens[-c(badidx)], ysens[-c(badidx)], zsens[-c(badidx)],
        xbad, ybad, zbad, t(x[, -c(badidx)])),
        silent = TRUE)
      if (!"try-error" %in% class(err)) break
      badidx <- badidx[-1]
    }
  }
  
  # put the interpolated data back in x
  if (nrow(ss) > 0) {
    for (i in seq_along(badidx)) {
      x[, badidx[i]] <- ss[i, ]}
  }
  
  attributes(x) <- atr
  x
  
}


spheric_spline <- function(xsens, ysens, zsens, xbad, ybad, zbad, values) {
  
  newsens <- length(xbad)
  numpoints <- ncol(values)

  Gsens <- computeg(xsens, ysens, zsens, xsens, ysens, zsens)
  Gsph  <- computeg(xbad, ybad, zbad, xsens, ysens, zsens)
  
  # compute solution for parameters C
  #----------------------------------
  meanvalues <- apply(values, 2, mean)
  # make mean zero
  values <- values - pracma::repmat(meanvalues, nrow(values), 1)
  values <- rbind(values, rep(0, numpoints))
  C <- pracma::pinv(rbind(Gsens, rep(1, max(nrow(Gsens), ncol(Gsens))))) %*% values

  res <- matrix(0, newsens, numpoints)
  # apply results
  # -------------
  for (j in seq_len(nrow(Gsph))) {
    res[j, ] <- apply(
      C * pracma::repmat(matrix(Gsph[j, ], ncol = 1), 1, ncol(C)),
      2,
      sum)
  }
  res <- res + pracma::repmat(meanvalues, nrow(res), 1)
  res
}

computeg <- function(x, y, z, xelec, yelec, zelec) {
  unitmat <- matrix(1, length(x), length(xelec))
  EI <- unitmat - sqrt(
    (pracma::repmat(matrix(x, ncol = 1), 1, length(xelec)) - 
       pracma::repmat(xelec, length(x), 1))^2 +
    (pracma::repmat(matrix(y, ncol = 1), 1, length(xelec)) -
       pracma::repmat(yelec, length(x), 1))^2 +
    (pracma::repmat(matrix(z, ncol = 1), 1, length(xelec)) -
       pracma::repmat(zelec, length(x), 1))^2
  )
  
  g <- matrix(0, length(x), length(xelec))
  
  m <- 4  # 3 is linear, 4 is best according to Perrin's curve
  for (n in seq_len(7)) {
    # pracma legendre function cannot process 2-D matrices
    for (icol in seq_len(ncol(EI))) {
      tmpL = pracma::legendre(n, EI[, icol])
      if (icol == 1) L = array(0, dim=c(nrow(tmpL), ncol(tmpL), ncol(EI)))
      L[, , icol] = tmpL
    }
    g <- g + ((2 * n + 1) / (n^m * (n + 1)^m)) * drop(L[1, , ])
  }
  g <- g / (4 * pi)    
  g
}

