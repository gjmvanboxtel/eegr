# noisysensors.R
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
# 20220204  GvB           Initial setup v0.3-2
# 20220315  GvB           changed round() to mround(), bugfix in
#                         randsample(), catch errors in RANSAC, 
#                         bugfix cor_dev, corMed, rcorr, ransac_frac
#---------------------------------------------------------------------------------------------------------------------

#' noisysensors
#' 
#' Detect noisy or outlier sensors
#' 
#' The algorithm uses four primary measures: extreme amplitudes (deviation
#' criterion), unusual high frequency noise (noisiness criterion), lack of
#' correlation with any other sensor (correlation criterion), and lack of
#' predictability by other sensors (predictability criterion), in this order.
#' Several of the criteria use a robust z score, replacing the mean by the
#' median and the standard deviation by the robust standard deviation (0.7413
#' times the interquartile range). The algorithm also detects sensors that
#' contain any NAs that have significant periods with constant values or very
#' small values.
#' 
#' The **deviation criterion** calculates the robust z-score of the
#' robust standard deviation for each sensor. Sensors designated as bad have
#' a robust z-score greater than \code{devthr} (default 5).
#' 
#' The **noisiness criterion** of signal quality uses a robust estimate of
#' the ratio of the power of the high frequency components to the power in the
#' low frequency components. A 50 Hz low pass FIR filter is used to separate the
#' low and high frequency components. Noisiness is defined as the ratio of the
#' median absolute deviation of the high frequency component over the low
#' frequency component for each sensor. A z-score relative to all of the sensors
#' is computed and sensors with a z-score greater than a threshold
#' (\code{HFnoisethr}; default: 5) are marked as bad.
#' 
#' The **correlation criterion** is based on the observation that the low
#' frequency portion of EEG signals is somewhat correlated (but not too
#' correlated) among channels. Using signals low- pass filtered at 50 Hz, the
#' algorithm calculates the correlation of each sensor with the other channels
#' in small, non-overlapping time windows (\code{corsecs} parameter; 1 s by
#' default). The maximum absolute correlation is calculated as the 98th
#' percentile of the absolute values of the correlations with the other sensors
#' in each window. If this maximum correlation is less than a threshold
#' (\code{corthr}; 0.4 by default) for a certain percentage of the windows
#' (\code{timethr}; 1% by default), then the sensor is classified as bad.
#' 
#' The **predictability criterion** also relies on correlations of the low
#' frequency portion of the signals. The RANSAC (random sample consensus) method
#' of Fischler and Bolles (1981) is used to select a random subset of (so far)
#' good sensors to predict the behavior of each sensor in small non-overlapping
#' time windows (\code{rwin}, 5 seconds by default). A random subset of
#' predictor sensors is chosen for each sensor (\code{sfrac}; 25% by default).
#' The RANSAC algorithm uses a method of spherical splines for estimating scalp
#' potentials. Sensors which have a correlation less than a threshold
#' (\code{ransacthr}; 0.75 by default) with their RANSAC-predicted time courses
#' on more than a certain fraction (\code{ubtime}; 0.4 by default) of the
#' windows are flagged as bad.
#' 
#' @param x input data, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}. The data 
#'   is assumed to be high-pass filtered.
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}. .
#' @param sensors sensors to detect, specified as positive integers indicating
#'   sensors numbers, or as sensor names (colnames of \code{x}). Default: all
#'   sensors names in \code{x} excluding a variable \code{t} (if present).
#' @param devthr robust deviation threshold, specified as a numeric value.
#'   Sensors with robust z-score greater than this value are marked as bad.
#'   Default: 5
#' @param HFnoisethr high-frequency noise threshold, specified as a numeric
#'   value. Sensors with a robust z-score calculated from the ratio between
#'   high- and low-frequency activity are marked as bad. Default: 5
#' @param corsecs segment length in seconds for computing correlations between
#'   sensors. The data is segmented in non-overlapping segments of this duration
#'   before calculating correlations. Default: 1 second
#' @param corthr correlation threshold, specified as a numeric value. Sensors
#'   which have a correlation with other channels lower than this value, for a
#'   proportion of the data specified by \code{timethr}, are flagged as
#'   bad. Default: 0.4
#' @param timethr time threshold, specified as a numeric value. Proportion of
#'   segments during which the correlation is allowed to be lower than
#'   \code{corthr} before the sensor is marked as bad. Default: 0.01
#' @param ransac logical indicating whether to perform RANSAC (random sample
#'   consensus). Default: TRUE
#' @param slocs sensor locations, needed for RANSAC, specified as a data frame
#'   according to \code{\link{sensorlocs}} format. Default: obtained from the
#'   sensor names given by the column names of \code{x}
#' @param sfrac fraction of sensors for robust reconstruction. Default: 0.25
#' @param ubtime RANSAC unbroken time - cutoff fraction of time a sensor can
#'   have poor RANSAC predictability. Default: 0.4
#' @param ransacthr RANSAC threshold. Sensors which with a correlation less
#'   this value with their RANSAC-predicted time courses on more than
#'   \code{ubtime} proportion of the windows are flagged as bad. Default: 0.75
#' @param rwin RANSAC window in seconds. Default: 5
#' @param rss RANSAC sample size. Default: 50
#' 
#' @return A list with the following fields:
#' \describe{
#'   \item{BadFromNA}{**Unusable**: sensors containing NA values}
#'   \item{BadFromConstant}{**Unusable**: sensors containing constant values
#'   (usually 0)}
#'   \item{devthr}{**Deviation criterion**: sensors having a robust z-score of
#'   the robust standard deviation greater than this value are classified as
#'   bad}
#'   \item{dev}{**Deviation criterion**: robust standard deviation for each
#'   sensor}
#'   \item{devSD}{**Deviation criterion**: robust SD of deviation scores}
#'   \item{devMed}{**Deviation criterion**: median of deviation scores}
#'   \item{Deviation}{**Deviation criterion**: robust z-scores of deviation
#'   scores}
#'   \item{BadFromDeviation}{**Deviation criterion**: sensors marked as bad by
#'   the deviation criterion}
#'   \item{HFnoisethr}{**Noisiness criterion**: sensors having a robust
#'   noisiness z-score greater than this value are classified as bad}
#'   \item{HFnoise}{**Noisiness criterion**: ratio of high over low frequency
#'   noise (per median absolute deviation) for each sensor}
#'   \item{medHFnoise}{**Noisiness criterion**: median of HFnoise values}
#'   \item{sdHFnoise}{**Noisiness criterion**: robust SD of HFnoise values}
#'   \item{zHFnoise}{**Noisiness criterion**: robust z-scores of HFnoise values}
#'   \item{BadFromHFnoise}{**Noisiness criterion**: sensors marked as bad by the
#'   noisiness criterion}
#'   \item{corsecs}{**Correlation criterion**: for the correlation criterion,
#'   each sensor is segmented into non-overlapping periods of this length
#'   (seconds)}
#'   \item{corthr}{**Correlation criterion**: sensors with a correlation less
#'   than this value are marked as bad}
#'   \item{timethr}{**Correlation criterion**: percentage of segments of a
#'   sensor for which the correlation is allowed to be less than the threshold}
#'   \item{cors}{**Correlation criterion**: matrix of size [number of segments,
#'   number of sensors] containing, for each segment, the maximum absolute
#'   correlation (98% percentile of the absolute value of a sensor's correlation
#'   with the other sensors)}
#'   \item{cor_noise}{**Correlation criterion**: matrix of size [number of
#'   segments, number of sensors] containing, for each segment and sensor, a
#'   value of the noise level, defined as the ratio of the difference between
#'   the original and the 50 Hz low-pass filtered data}
#'   \item{cor_dev}{**Correlation criterion**: robust deviation scores for
#'   cor_noise}
#'   \item{corMed}{**Correlation criterion**: median of correlations per sensor}
#'   \item{BadFromCorrelation}{**Correlation criterion**: sensors marked as bad
#'   by the correlation criterion}
#'   \item{BadFromDropout}{**Correlation criterion**: sensors marked as bad by
#'   the cor_noise noisiness criterion (dropouts)}
#'   \item{ransacPerformed}{**Predictability criterion**: logical indicating
#'   whether RANSAC was performed}
#'   \item{ransacSensors}{**Predictability criterion**: sensors on which RANSAC
#'   was performed}
#'   \item{sfrac}{**Predictability criterion**: fraction of sensors on which
#'   RANSAC prediction was based}
#'   \item{ubtime}{**Predictability criterion**: RANSAC unbroken time - cutoff
#'   fraction of time a sensor can have poor RANSAC predictability}
#'   \item{ransacthr}{**Predictability criterion**: RANSAC threshold. Sensors
#'   which with a correlation less this value with their RANSAC-predicted time
#'   courses on more than ubtime proportion of the windows are flagged as bad}
#'   \item{rwin}{**Predictability criterion**: RANSAC window length in seconds}
#'   \item{rss}{**Predictability criterion**: RANSAC sample size}
#'   \item{rcorr}{**Predictability criterion**: matrix of size [number of
#'   sensors, number of windows] containing, for each segment (window) and
#'   sensor, the RANSAC correlations}
#'   \item{ransac_frac}{**Predictability criterion**: fraction of bad RANSAC
#'   windows}
#'   \item{BadFromRansac}{**Predictability criterion**: sensors marked as bad
#'   from the predictability criterion}
#'   \item{TotalBad}{**Summary**: indices of all bad sensors}
#'   \item{BadNames}{**Summary**: names of bad sensors}
#' }
#' 
#' @examples
#' \dontrun{ 
#' noisy <- noisysensors(EEGdata)
#' }
#'   
#' @author Methods 1 and 4 are adapted from code by Christian Kothe and Methods
#'   2 and 3 are adapted from code by Nima Bigdely-Shamlo; ported to R and
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com},
#' 
#' @references Bigdely-Shamlo, N., Mullen, T.,, Kothe, C.,, Su, K.-M., and
#'   Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for
#'   large-scale EEG analysis. Frontiers in Neuroinformatics, 9, Article 16,
#'   \url{https://www.frontiersin.org/article/10.3389/fninf.2015.00016}.
#' @references Fischler, M.A., and Bolles, R.C. (1981). Random sample consensus:
#'   A paradigm for model fitting with apphcatlons to image analysis and
#'   automated cartography. Communications of the ACM, 24(6)), 381-395.
#' 
#' @export

noisysensors <- function(x, fs = 1, sensors = setdiff(colnames(x), "t"),
                         devthr = 5,
                         HFnoisethr = 5,
                         corsecs = 1, corthr = 0.4, timethr = 0.01,
                         ransac = TRUE, slocs = eegr::getLocationsfromLabels(colnames(x)[sensors]),
                         sfrac = 0.25, ubtime = 0.4, ransacthr = 0.75, rwin = 5, rss = 50) {
  
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
  } else {
    if (!(is.numeric(fs) && length(fs) == 1 && fs > 0)) {
      stop("Sampling frequency (fs) must be a positive numeric value")
    }
  }
  if (is.vector(x)) {
    cn <- names(x)
    x <- matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    cn <- colnames(x)
    vec <- FALSE
  }
  npts <- nrow(x)
  ns <- ncol(x)
  
  # work with sensor numbers rather than names
  if (is.character(sensors)) {
    sensors <- sort(which(sensors %in% colnames(x)))
  }
  if (!is.numeric(sensors)) {
    stop("sensors must be specified as names or numbers")
  }
  if (!is.numeric(devthr) || length(devthr) != 1L) {
    stop('deviation threshold must be a numeric scalar value')
  } else {
    devthr <- abs(devthr)
  }
  if (!is.numeric(HFnoisethr) || length(HFnoisethr) != 1L) {
    stop('high-frequency noise threshold must be a numeric scalar value')
  }
  if (!is.numeric(corsecs) || length(corsecs) != 1L || corsecs <= 0) {
    stop('high-frequency noise threshold must be a positive numeric scalar value')
  }
  if (!is.numeric(corthr) || length(corthr) != 1L) {
    stop('high-frequency noise threshold must be a numeric scalar value')
  }
  if (!is.numeric(timethr) || length(timethr) != 1L || timethr <= 0) {
    stop('time threshold must be a positive numeric scalar value')
  }
  if (!is.logical(ransac) || length(ransac) != 1L) {
    stop('ransac must be a logical scalar value')
  } 
  if (ransac) {
    if (!is.numeric(sfrac) || length(sfrac) != 1L || sfrac <= 0) {
      stop('time threshold must be a positive numeric scalar value')
    }
    if (!is.numeric(ubtime) || length(ubtime) != 1L || ubtime <= 0) {
      stop('time threshold must be a positive numeric scalar value')
    }
    if (!is.numeric(ransacthr) || length(ransacthr) != 1L) {
      stop('ransac threshold must be a numeric scalar value')
    }
    if (!is.numeric(rwin) || length(rwin) != 1L || rwin <= 0) {
      stop('ransac window must be a positive numeric scalar value')
    }
  }
  
  # Detect unusable sensors
  orig_sensors <- sensors
  bad <- unusablesensors(x[, sensors])
  if (bad$n > 0) {
    sensors[which(sensors %in% bad$unusable)] <- NA
  }
  good <- which(!is.na(sensors))
  
  #############################################################################
  # Method 1: Unusually high or low amplitude (using robust std)

  dev <- apply(x[, good], 2, function(x) 0.7413 * stats::IQR(x))  # Robust estimate of SD
  devSD <-  0.7413 * stats::IQR(dev)
  devMed <- stats::median(dev, na.rm = TRUE)
  Deviation <- rep(NA, ns)
  names(Deviation) <- cn
  Deviation[good] <- (dev - devMed) / devSD
  noisyDev <- which(!is.na(Deviation) & abs(Deviation) > devthr)

  #############################################################################
  # Method 2: Compute the SNR (based on Christian Kothe's clean_channels)
  # Note: RANSAC uses the filtered values X of the data

  HFnoise <- zHFnoise <- rep(NA, ns)
  names(zHFnoise) <- cn
  if (fs > 100) {
    # Remove signal content above 50Hz (assume HP filtering already done)
    nyquist <- fs / 2
    h <- gsignal::fir2(100, c(0, 45, 50, nyquist) / nyquist, c(1, 1, 0, 0))
    X <- apply(x, 2, function(x, h) gsignal::filtfilt(h, x), h)
    HFnoise <- apply(x[, good], 2, function(x, X, good) 
      stats::mad(x - X[, good]) / stats::mad(X[, good]), X, good)
    medHFnoise <- stats::median(HFnoise, na.rm = TRUE)
    sdHFnoise <- stats::mad(HFnoise) * 1.4826
    zHFnoise[good] <- (HFnoise - medHFnoise) / sdHFnoise
    noisyHF <- which(!is.na(zHFnoise) & zHFnoise > HFnoisethr)
  } else {
    X = x
    medHFnoise <- 0
    sdHFnoise <- 1
    zHFnoise[sensors] <- rep(0, length(sensors))
    names(zHFnoise) <- cn
    noisyHF <- integer(0)
  }                        

  #############################################################################
  # Method 3: Global correlation criteria (from Nima Bigdely-Shamlo)
  cor_seglen <- corsecs * fs
  cor_offsets <- seq(1, (npts - cor_seglen), cor_seglen)
  nwin <- length(cor_offsets)

  cors <- cor_noise <- cor_dev <- matrix(NA, nwin, ns)
  colnames(cors) <- colnames(cor_noise) <- colnames(cor_dev) <- cn
  corMed <- frac_cor <- frac_do <- rep(NA, ns)
  names(corMed) <- names(frac_cor) <- names(frac_do) <- cn

  for (w in seq_len(nwin)) {
    o <- cor_offsets[w]
    X_seg <- X[o:(o + cor_seglen - 1), good]          #low-pass filtered data
    x_seg <- x[o:(o + cor_seglen - 1), good]          #original data
    if (all(X_seg == x_seg)) {
      warning("HF and LF data are identical")
    } else {
      wincor <- suppressWarnings(stats::cor(X_seg))     #correlation matrix
      abscor <- abs(wincor - diag(diag(wincor)))        #diagonal always contains 1s
      cors[w, good] <- apply(abscor, 2, stats::quantile, 0.98, na.rm = TRUE)
      cor_noise[w, good] <- apply(x_seg, 2, function(x, X) 
        stats::mad(x - X) / stats::mad(X), X_seg)
      cor_dev[w, good] <- apply(x_seg, 2, function(x) 0.7413 * stats::IQR(x))
    }
  }
  dropouts <- is.na(cors[, good]) | is.na(cor_noise[, good])
  thr_cors <- cors[, good] < corthr
  frac_cor[good] <- apply(thr_cors, 2, mean)
  frac_do[good] <- apply(dropouts, 2, mean)
  
  noisyCor <- which(frac_cor > timethr)
  noisyDO <- which(frac_do > timethr)
  med <- apply(cors[, good], 2, stats::median)
  corMed[good] <- med
  
  #############################################################################
  # Bad so far by amplitude and correlation (take these out before doing ransac)
  noisy = union(noisyDev, union(noisyHF, union(noisyCor, noisyDO)))
  
  
  #############################################################################
  # Method 4: Ransac corelation (may not be performed)
  
  # Setup for ransac (if a 2-stage algorithm, remove other bad channels first)
  ransac_frac <- 0
  ransacPerformed <- TRUE
  if (!ransac) {
    ransacPerformed <- FALSE
  } else if (is.null(slocs) || nrow(slocs) <= 0) { 
    warning("ransac could not be computed because there were no sensor locations")
    ransacPerformed <- FALSE
  } else if (!any(c("x", "y", "z") %in% colnames(slocs))) {
    warning("slocs must contain valid x, y and z variables")
    ransacPerformed <- FALSE
  } else {
    rsens <- setdiff(good, noisy)
    X <- X[, rsens]
  
    ransacSS <- mround(sfrac * ns) #rounding is different in Matlab and R
    if (ubtime < 1) {
      ubframes <- npts * ubtime
    } else {
      ubframes = fs * ubtime
    }
  
    rFrames <- rwin * fs
    rWindow <- seq(0, rFrames - 1)
    rOffsets <- seq(1, max(1, (npts - rFrames)), rFrames)
    WRansac <- length(rOffsets)
    
    nlocs <- nrow(slocs)
    if (nlocs < 3 || length(rsens) < 2) {
      warning("too many channels have failed quality tests to perform ransac")
      ransacPerformed <- FALSE
    }
  }

  noisyRansac <- integer(0)
  if (ransacPerformed) {
    err <- try({
      locs <- slocs[rsens, c("x", "y", "z")]
      P <- calc_projector(locs, rss, ransacSS)
      corrT <- matrix(0, nrow(locs), WRansac)

      # Calculate each channel's correlation to its RANSAC reconstruction for each window
      n <- length(rWindow)
      m <- length(rsens)
      p <- rss
      Xwin <- X[1:(n*WRansac), ]
      dim(Xwin) <- c(m, n, WRansac)
      for (k in seq_len(WRansac)) {
        corrT[, k] <- calculateRansacWindow(t(Xwin[, , k]), P, n, m, p)
      }
      flagged <- corrT < ransacthr
      if (length(flagged) > 0) {
        noisyRansac = which(sum(flagged) * rFrames > ubframes)
        ransac_frac = apply(flagged, 1, sum) / length(flagged)
        names(ransac_frac) <- colnames(x)[rsens]
      } else {
        noisyRansac <- integer(0)
        ransac_frac <- NA
      }
    }, silent = TRUE)
    if ("try-error" %in% err) {
      ransacPerformed <- FALSE
    }
  } else {
    rsens <- sfrac <- ubtime <- ransacthr <- rwin <- rss <- corrT <- ransac_frac <- NA
    noisyRansac <- integer(0)
  }

  # Combine bad channels detected from all methods
  noisy = sort(union(bad$unusable, union(noisy, noisyRansac)))
  
  # Set up return values
  rv <- list(
    # unusable
    BadFromNA = bad$sensor_NA,
    BadFromConstant = bad$sensor_constant,
    # 1: deviation
    devthr = devthr,
    dev = dev,
    devSD = devSD,
    devMed = devMed,
    Deviation = Deviation[sensors],
    BadFromDeviation = noisyDev,
    # 2: High Frequency Noise
    HFnoisethr = HFnoisethr,
    HFnoise = HFnoise,
    medHFnoise = medHFnoise,
    sdHFnoise = sdHFnoise,
    zHFnoise = zHFnoise[sensors],
    BadFromHFnoise = noisyHF,
    # 3: correlation
    corsecs = corsecs,
    corthr = corthr,
    timethr = timethr,
    cors = cors[, sensors],
    cor_noise = cor_noise[, sensors],
    cor_dev = cor_dev[, sensors],
    corMed = corMed[sensors],
    BadFromCorrelation = noisyCor,
    BadFromDropout = noisyDO,
    # 4: RANSAC
    ransacPerformed = ransacPerformed,
    ransacSensors = rsens,
    sfrac = sfrac,
    ubtime = ubtime,
    ransacthr = ransacthr,
    rwin = rwin,
    rss = rss,
    rcorr = corrT,
    ransac_frac = ransac_frac,
    BadFromRansac = noisyRansac,
    # summary
    TotalBad = noisy,
    BadNames = cn[noisy]
  )
  
  rv

}

# helper functions for RANSAC

calc_projector <- function(locs, numberSamples, subsetSize) {
  # Calculate a bag of reconstruction matrices from random channel subsets
  
  ps <- getRandomSubsets(locs, subsetSize, numberSamples)
  permutedLocations <- ps$permutedLocations
  subsets <- ps$subsets
  randomSamples = vector("list", numberSamples)
  for (k in seq_len(numberSamples)) {
    tmp <- matrix(0, nrow(locs), nrow(locs))
    slice <- subsets[k, ]
    tmp[slice, ] <- Re(ssinterp(t(permutedLocations[, , k]), locs))
    randomSamples[[k]] <- tmp
  }
  P = do.call(cbind, randomSamples)
  P
}

getRandomSubsets <- function(locs, subsetSize, numberSamples) {
  #stream = RandStream('mt19937ar', 'Seed', 435656);
  numberChannels <- nrow(locs)
  permutedLocations <- array(0, dim = c(3, subsetSize, numberSamples))
  subsets = matrix(0, numberSamples, subsetSize)
  for (k in seq_len(numberSamples)) {
    subset <- randsample(seq_len(numberChannels), subsetSize)
    subsets[k, ] <- subset
    permutedLocations[1:3, 1:subsetSize,  k] <- t(locs[subset, ])
  }
  list(permutedLocations = permutedLocations, subsets = subsets)
}

randsample <- function (X, num) {
  Y <- rep(0, num)
  for (k in seq_len(num)) {
    # make sure pick > 0
    while (TRUE) {
      pick <- mround(1 + (length(X) - 1) * stats::runif(1))
      if (pick > 0) break
    }
    Y[k] <- X[pick]
    X <- X[-pick]
  }
  Y
}

calculateRansacWindow <- function(XX, P, n, m, p) {
  YY <- XX %*% P
  dim(YY) <- c(n, m, p)
  YY = aperm(apply(YY, 1:2, sort), c(2,3,1))
  YY = YY[, , mround(dim(YY)[3] / 2)]
  rX = apply(XX * YY, 2, sum) /
    (sqrt(apply(XX^2, 2, sum)) * sqrt(apply(YY^2, 2, sum)))
  rX
}
