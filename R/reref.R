# reref.R
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
# 20220208  GvB           Initial setup v0.3-2                  
# 20220210  GvB           In robust referencing, do not interpolate
#                         if number of bad sensors > 50%
#---------------------------------------------------------------------------------------------------------------------

#' Rereference Data
#'
#' Compute the average reference or a new common reference for data.
#' 
#' The function \code{reref} computes the mean of the new reference sensors and
#' subtracts that from the channels specified in \code{sensors}. If the data
#' were originally referenced with respect to a sensor that is not present in
#' the data, then this sensor will be added, and it can be used in the new
#' reference.
#' 
#' Robust statistics can be used to compute the new reference (median or Huber
#' mean instead of arithmetic mean). In addition, \code{reref} can use an
#' iterative procedure to remove bad channels before or after re-referencing,
#' which makes the re-referencing even more robust to outlier sensors. The
#' detection of noisy sensors is done using the \code{\link{noisysensors}}
#' function, using the default setup.
#' 
#' @param x input time series, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}.
#' @param sensors vector of sensor names or numbers that are used for the
#'   calculations (hence they must contain EEG data). Default: all.
#' @param oldref vector of sensor names or numbers that formed the original
#'   reference. If specified, the old reference signal is reconstructed back
#'   into the data object (added if it was not present in the data). Ignored if
#'   \code{x} is an object of class \code{ctd} and contains a \code{"ref"}
#'   attribute. Default: NULL.
#' @param newref signals that form the new reference, either specified as a
#'   character vector of signal names (\code{colnames(ctd_obj)}), or as a vector
#'   of numerical values. If more than one signal is specified, then the average
#'   of these signals will form the new reference. Alternatively, specifying
#'   \code{"robust"} (default) computes an average reference with iterative
#'   detection and interpolation of bad channels, and \code{"average"} removes
#'   the average of the reference channels with no interpolation.
#' @param saveref logical indication whether to save to new reference as a 
#    separate sensor channel named "Ref" . Default: FALSE
#' @param interp For \code{newref = "robust"}, specifies whether to perform
#'   final channel interpolation. If \code{"post-reference"} (Dafault), after a
#'   final robust reference is computed, the channels are re-interpolated and
#'   the reference is corrected. In \code{"pre-reference"}, reref incrementally
#'   adds to the bad channel list and interpolates before computing the
#'   reference. If the initial estimate of the reference is poor, this is not a
#'   good approach. In the \code{"none"} option, reref removes the reference but
#'   does not interpolate. Bad channels remain in the signal. You may choose to
#'   remove them during post-processing.
#' @param estmean For \code{newref = "robust"}, statistic used for estimation of
#'   the mean; one of \code{"median"} (default), \code{"huber"} or
#'   \code{"mean"}. Median and mean are the obvious statistics; \code{"huber"}
#'   uses the Huber M-estimator of location with median absolute deviation (MAD)
#'   scale, which is computationally expensive.
#' @param sl For \code{newref = "robust"}, sensor locations of \code{x},
#'   specified as a data frame according to \code{\link{sensorlocs}} format.
#'   Alternatively, a matrix with named columns x, y, and z representing the
#'   sensor locations. Default: obtained from the sensor names given by the
#'   column names of \code{x}. Only used in case \code{"when"} is either
#'   \code{"before"} or \code{"after"}.
#' @param maxIter maximum number of iterations used for detecting noisy channels
#'   and recomputing the reference. Default: 4. Only used in case \code{"newref"}
#'   is \code{"robust"}.
#' 
#' @return A list containing the following elements:
#' \describe{
#'   \item{y}{An object of the same class as the input, containing the
#'   re-referenced data. Signals not included in \code{sensors} are copied to
#'   the output object. If \code{saveref = TRUE} then the new reference is added
#'   as a separate data channel (column). If the old reference was specified
#'   through the \code{oldref} parameters, and it was not also present in the
#'   data, then it is also added.}
#'   \item{fs}{the sampling frequency of \code{y}}
#'   \item{sensors}{the sensors that \code{reref} operated on (on input)}
#'   \item{oldref}{the original reference sensor(s)}
#'   \item{newref}{the new reference sensor(s)}
#'   \item{saveref}{logical indicating whether the new reference is saved}
#'   \item{estmean}{for robust referencing, the initial estimate of the mean}
#'   \item{interp}{for robust referencing the instant of interpolating noisy
#'   sensors, either before or after re-referencing}
#'   \item{iterations}{for robust referencing, the number of iterations
#'   performed to detect noisy reference sensors}
#'   \item{bad}{a list of sensors that were marked as bad by different
#'   criteria (see \code{\link{noisysensors}})}
#' }
#'
#' @references Bigdely-Shamlo, N., Mullen, T.,, Kothe, C.,, Su, K.-M., and
#'   Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for
#'   large-scale EEG analysis. Frontiers in Neuroinformatics, 9, Article 16,
#'   \url{https://www.frontiersin.org/article/10.3389/fninf.2015.00016}.
#' @references Huber, P.J. (1981) Robust Statistics. Wiley.
#'
#' @examples
#'
#' avgref <- reref(EEGdata, sensors = 1:29, oldref = "M1",
#'                 newref = "average")
#' lmref <- reref(EEGdata, sensors = 1:29, oldref = "M1", 
#'                 newref = "M1 M2")
#' 
#' \dontrun{
#' rref <- reref(EEGdata, 1:29, "M1", "robust")
#' }
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}, based on Matlab
#'   code of the PREP pipeline by Nima Bigdely-Shamlo and colleagues.
#'
#' @seealso \code{\link{noisysensors}}, \code{\link[MASS]{huber}}
#' 
#' @export

reref <- function(x, fs = 1, sensors = setdiff(colnames(x), "t"), 
                  oldref = NULL, newref = "robust", saveref = FALSE,
                  interp = c("post-reference", "pre-reference", "none"),
                  estmean = c("median", "huber", "mean"),
                  sl = getLocationsfromLabels(colnames(x[, sensors])),
                  maxIter = 4
                  ) {

  # parameter checking
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
    npts <- npts(x)
    ns <- ns(x)
    if (!is.null(attr(x, "ref"))) {
      oldref <- attr(x, "ref")
    }
  } else if (!is.matrix(x)) {
    stop("x must be a ctd object or a matrix)")
  } else {
    npts <- nrow(x)
    ns <- ncol(x)
  }
  
  # set sensor names if not present to prevent confusion on output
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("S", seq_len(ns))
  }

  # all sensors in x? also separate used and unused columns of x
  if (is.numeric(sensors)) {
    if (max(sensors) > ns) {
      stop("invalid numeric sensors specified")
    } else {
      used <- colnames(x[, sensors])
      notused <- colnames(x[, setdiff(seq_len(ncol(x)), sensors)]) #use ncol not ns
    }
  } else if (is.character(sensors)) {
    sensors <- splitsens(sensors)
    if (!all(sensors %in% colnames(x))) {
      stop("invalid character sensors specified")
    } else {
      used <- sensors
      notused <- setdiff(colnames(x), sensors)
    }
  } else {
    stop("sensors must be a numeric or character vector")
  }
  
  # oldref not in x? then add new channel(x)
  extra <- NULL
  orig_used <- used
  if (!is.null(oldref)) {
    if (is.numeric(oldref)) {
      n <- length(which(oldref > ns))
      if (n > 0) {
        extra <- matrix(0, npts, n)
        xu <- seq(length(used), length(used) + n)
        colnames(extra) <- paste0("S", xu)
        used <- c(used, xu)
      }
    } else if (is.character(oldref)) {
      oldref <- splitsens(oldref)
      names <- setdiff(oldref, sensors)
      n <- length(names)
      if (n > 0 && !(n == 1 && (names == "average" || names == "avg"))) {
        extra <- matrix(0, npts, n)
        colnames(extra) <- names
        used <- c(used, names)
      }
    } else {
      stop("oldref must be numeric or character")
    }
  }
  data <- cbind(x[, orig_used], extra)
  
  if (is.numeric(newref)) {
    if (any(newref > ncol(data))) {
      stop("invalid numeric new reference sensors specified")
    }
  } else if (is.character(newref)) {
    if (length(newref) == 1) {
      if (newref == "robust") {
        robust <- TRUE
        newref <- used
        interp <- match.arg(interp)
        estmean <- match.arg(estmean)
      }
      else if (newref == "average" || newref == "avg") {
        robust <- FALSE
        newref <- used
      } else {
        robust <- FALSE
        newref <- splitsens(newref)
      }
    } else {
      newref <- splitsens(newref)
    }
  } else {
    stop("newref must be numeric or character")
  }
  if (!all(newref %in% used)) {
    stop("invalid character new reference sensors specified")
  }
  
  if (!is.logical(saveref)) {
    saveref <- FALSE
  }

  # do the rereferencing - always work with sensor names
  ######################################################
  if (robust) {
    
    # Determine unusable channels and remove them from the reference channels
    signal <- removetrend(data, fs)$y
    noisyOrig <- noisysensors(signal, fs, newref, ransac = FALSE)
    bad <- updatebad(NULL, noisyOrig)
    bad$BadFromCorrelation <- bad$BadFromDropout <- integer(0)
    bad$TotalBad <- unique(c(bad$BadFromNA, bad$BadFromConstant,
                             bad$BadFromDeviation, bad$BadFromHFnoise))
    bad$BadNames <- newref[bad$TotalBad]
    ref <- setdiff(newref, bad$BadNames)
    
    # if more than 50% of the sensors are bad, then compute
    # a normal averaged reference
    if (length(bad$BadNames) / length(sensors) > 0.5) {
      refdata <- apply(data[, newref], 1, mean)
      data <- sweep(data, 1, refdata)
      estmean <- interp <- NULL
      iterations <- 0
    } else {
    
      # Get initial estimate of the mean by the specified method
      refTemp <- apply(signal[, ref], 1, estmean)
      signalTemp <- sweep(signal, 1, refTemp)
    
      # Remove reference from signal iteratively interpolating bad channels
      iterations <- 0
      noisyOld = NULL
      while (TRUE) {    # Do at least 1 iteration
        noisy <- noisysensors(signalTemp, fs)
        bad <- updatebad(bad, noisy)
        if (iterations > 1 && (length(noisy$TotalBad) <= 0 ||
            (length(setdiff(noisy$TotalBad, noisyOld$TotalBad)) <= 0 &&
             length(setdiff(noisyOld$TotalBad, noisy$TotalBad)) <= 0)) ||
            iterations > maxIter) {
          break
        }
        noisyOld <- noisy
        srcSens <- setdiff(ref, noisy$BadNames)
        if (length(srcSens) < 2) {
          stop("Could not perform a robust reference - not enough good channels")
        }
        if (length(noisy$BadNames) > 0) {
          signalTemp <- eeginterp(signal, bad = noisy$BadNames)
        } else {
          signalTemp <- signal
        }
        refSignal <- apply(signalTemp[, ref], 1, mean, na.rm = TRUE)
        signalTemp <- sweep(signal, 1, refSignal)
        iterations <- iterations + 1
      }

      if (interp == "pre-reference") {
        # Use the bad channels accumulated from reference search to robust
        if (length(bad$BadNames) <= 0) {
          refdata <- apply(data[, newref], 1, mean, na.rm = TRUE)
        } else {
          data <- eeginterp(data, bad = noisy$BadNames)
          refdata <- apply(data[, newref], 1, mean, na.rm = TRUE)
        }
        data <- sweep(data, 1, refdata)
      } else if (interp == "post-reference") {
        # Robust reference with interpolation afterwards
        if (length(bad$BadNames) <= 0) {
          refdata <- apply(data[, newref], 1, mean, na.rm = TRUE)
        } else {
          dataX <- eeginterp(data, bad = bad$BadNames)
          refdata <- apply(dataX[, newref], 1, mean, na.rm = TRUE)
          rm(dataX)
        }
        data <- sweep(data, 1, refdata)
        noisy <- noisysensors(removetrend(data, fs)$y, fs, newref)
        bad <- updatebad(NULL, noisy)
        
        # Bring forward unusable channels from original data
        bad <- updatebad(bad, noisyOrig)
        data <- eeginterp(data, bad = bad$BadNames)
        refdata <- apply(data[, newref], 1, mean, na.rm = TRUE)
        data <- sweep(data, 1, refdata)

      }
    }
  } else {
    # Do some type of non-robust referencing
    refdata <- apply(data[, newref], 1, mean)
    data <- sweep(data, 1, refdata)
    estmean <- interp <- iterations <- bad <- NULL
  }
  
  # save the reference if requested
  if (saveref) {
     data <- cbind(data, Ref = refdata)
  }
  
  # put the unused channels back 
  data <- cbind(data, x[, notused])
  
  # set attributes if x was ctd object
  if ("ctd" %in% class(x)) {
    class(data) <- class(x)
    attr(data, "npts") <- attr(x, "npts")
    ns <- attr(x, "ns")
    if (!is.null(ncol(extra))) ns <- ns + ncol(extra)
    if (saveref) ns <- ns + 1
    attr(data, "ns") <- ns
    attr(data, "fs") <- attr(x, "fs")
    attr(data, "ref") <- paste(newref, collapse = " ")
  }
  
  # return values
  list(
    y = data,
    fs = fs,
    sensors = sensors,
    oldref = oldref,
    newref = newref,
    saveref = saveref,
    estmean = estmean,
    interp = interp,
    iterations = iterations,
    bad = bad
  )
}

# function to calculate Huber mean
# x: M x N matrix
# calls MASS:huber
huber <- function(x) {
  h <- MASS::huber(x)
  h$mu
}

# function to update bad channels
updatebad <- function(bad = NULL, noisy) {
  
  if (is.null(bad)) {
    b <- list(
      BadFromNA = noisy$BadFromNA,
      BadFromConstant = noisy$BadFromConstant,
      BadFromDeviation = noisy$BadFromDeviation,
      BadFromHFnoise = noisy$BadFromHFnoise,
      BadFromCorrelation = noisy$BadFromCorrelation,
      BadFromDropout = noisy$BadFromDropout,
      BadFromRansac = noisy$BadFromRansac,
      TotalBad = noisy$TotalBad,
      BadNames = noisy$BadNames
    )
  } else {
    b <- list(
      BadFromNA = union(bad$BadFromNA, noisy$BadFromNA),
      BadFromConstant = union(bad$BadFromConstant, noisy$BadFromConstant),
      BadFromDeviation = union(bad$BadFromDeviation, noisy$BadFromDeviation),
      BadFromHFnoise = union(bad$BadFromHFnoise, noisy$BadFromHFnoise),
      BadFromCorrelation = union(bad$BadFromCorrelation, noisy$BadFromCorrelation),
      BadFromDropout = union(bad$BadFromDropout, noisy$BadFromDropout),
      BadFromRansac = union(bad$BadFromRansac, noisy$BadFromRansac),
      TotalBad = union(bad$TotalBad, noisy$TotalBad),
      BadNames = union(bad$BadNames, noisy$BadNames)
    )
  }
  b
}
