# cleanline.R
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
# Version history:
# 20211222  GvB           Initial setup v0.3-3
#---------------------------------------------------------------------------------------------------------------------

#' Cleanline
#' 
#' Estimate and remove sinusoidal (e.g. line) noise from EEG  channels using
#' multi-tapering and a Thompson F-statistic.
#' 
#' Sinusoidal noise can be a prominent artifact in recorded electrophysiological
#' data. This can stem from AC power line fluctuations (e.g. 50/60 Hz line noise
#' + harmonics), power suppliers (e.g. in medical equipment), fluorescent
#' lights, etc. Notch filtering is generally undesirable due to creation of
#' band-holes, and significant distortion of frequencies around the notch
#' frequency (as well as phase distortion at other frequencies and Gibbs
#' rippling in the time-domain). Other approaches for sinusoidal ("line") noise
#' reduction include adaptive regressive filtering approaches (e.g. RLS, LMS),
#' but these typically require a reference signal (e.g. externally-recorded
#' noise reference), which is often unavailable. Blind-source separation
#' approaches such as ICA may help mitigate line noise, but often fail to
#' completely remove the noise due to spatiotemporal non-stationarities in the
#' noise.
#' 
#' CleanLine uses an approach for line noise removal advocated by Partha Mitra
#' and Hemant Bokil in "Observed Brain Dynamics" (2007), Chapter 7.3.4.
#' 
#' In brief, the data is traversed by a sliding window. Within each window, the
#' signal is transformed to the frequency domain using a multi-taper FFT. The
#' complex amplitude (amplitude and phase) is thus obtained for each frequency.
#' Under the assumption of a deterministic sinusoid embedded in white noise, we
#' can set up a regression of the multi-taper transform (spectrum) of this
#' sinusoidal signal onto the multitaper spectrum of the original data at a
#' given frequency. The regression coefficient is a complex number representing
#' the complex amplitude (phase and amplitude) of the deterministic sinusoid.
#' From this, a time-domain representation of the sinusoid may be constructed
#' and subtracted from the data to remove the line.
#' 
#' Typically, one does not know the exact line frequency. For instance, in the
#' U.S.A., power line noise is not guaranteed to be at exactly 60 Hz (or even to
#' have constant phase over a given period of time). To ameliorate this problem
#' a Thompson F-Test may be applied to determine the statistical significance of
#' a non-zero coefficient in the above regression (indicating a sinusoid with
#' significantly non-zero amplitude). We can then search within a narrow band
#' around the expected location of the line for the frequency which maximizes
#' this F-statistic above a significance threshold (e.g. p=0.05).
#' 
#' Overlapping short (e.g. 2-4 second) windows can be used to adaptively
#' estimate the frequency, phase, and amplitude of the sinusoidal noise
#' components (which typically change over the course of a recording session).
#' The discontinuity at the point of window overlap can be smoothed using a
#' sigmoidal function.
#' 
#' @param x input time series, specified as a numeric or complex vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}
#' @param freq frequencies to remove. Will be expended to include harmonics up
#'   to the Nyquist frequency. Default: 50
#' @param fbw frequency bandwidth. Bandwidth centered on each \code{freq} to
#'   scan for significant lines. Default: 2 Hz.
#' @param detrend character string specifying detrending option; one of:
#'   \describe{
#'     \item{\code{long-mean}}{remove the mean from the data before
#'     splitting into segments (default)}
#'     \item{\code{short-mean}}{remove the mean value of each segment}
#'     \item{\code{long-linear}}{remove linear trend from the data before
#'     splitting into segments}
#'     \item{\code{short-linear}}{remove linear trend from each segment}
#'     \item{\code{none}}{no detrending}
#'  }
#' @param w length of sliding window (segments) in seconds. Default: 4
#' @param overlap proportion of overlap between the segments. Default: 0
#' @param tbw taper bandwidth, specified as a positive numeric value. Default: 2
#'   Hz
#' @param k number of tapers to use, specified as a positive numeric value.
#'   Default: \code{tbw * w - 1}
#' @param p significance level cutoff for F-test. Default 0.01
#' @param tau Window overlap smoothing factor. Default: 100
#' @param maxIter Maximum times to iterate removal. Default: 10
#' 
#' @return An object of the same size and class as the input, containing the
#'   cleaned data.
#'   
#' @examples
#' data(EEGdata)
#' ospec <- mtspec(EEGdata[, 1:28], fs = fs(EEGdata), detrend = "short-linear")
#' cleaned <- cleanline(EEGdata[, 1:28], fs = fs(EEGdata))
#' cspec <- mtspec(cleaned, fs = fs(EEGdata), detrend = "short-linear")
#' plot(ospec[, "f"], 10 * log10(ospec[, "Pz"]), type = "l",
#'       main = "Multitaper Spectrum\nEEGdata - Pz",
#'       xlab = "Frequency (Hz)", ylab = "Power/Frequency (dB)")
#' lines(cspec[, "f"], 10 * log10(cspec[, "Pz"]), col = "red")
#' legend("topright", legend = c("Original", "Cleaned"), lty = 1,
#'        col = c("black", "red"))
#' @seealso \code{\link{mtspec}}
#' 
#' @author Original Matlab code by Tim Mullen (2011); ported to R and
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com},
#'
#' @rdname cleanline   
#' @export

cleanline <- function(x, fs = 1, freq = 50, fbw = 2,
                      detrend = c("short-linear", "short-mean",
                                  "long-linear", "long-mean", "none"),
                      w = 4, overlap = 0,
                      tbw = 2, k = round(tbw * w - 1),
                      p = 0.01, tau = 100, maxIter = 10) {
  
  # parameter checking
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
    if ("t" %in% colnames(x)) {
      x <- x[, -which(colnames(x) == "t")]
    }
  } else {
    if (!(is.numeric(fs) && length(fs) == 1 && fs > 0)) {
      stop("fs must be a positive numeric value")
    }
  }
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
    vec = TRUE
  } else {
    vec = FALSE
  }
  if (!(is.numeric(freq) && all(freq > 0))) { #allow multiple values
    stop("freq must be a positive numeric value")
  }
  # expand frequencies up to Nyquist
  Nyq <- fs / 2
  while (TRUE) {
    freq <- unique(c(freq, 2 * freq))
    if (any(freq > Nyq)) break;
  }
  rmv  <- which(freq > Nyq)
  freq <- freq[-rmv]
  if (is.null(freq)) stop("No frequencies to correct below the Nyquist frequency")
  
  if (!(is.numeric(fbw) && length(fbw) == 1 && fbw > 0)) {
    stop("fbw must be a positive numeric value")
  }
  detrend <- match.arg(detrend)
  if (!(is.numeric(w) && length(w) == 1 && w > 0)) {
    stop("w must be a positive numeric value")
  }
  if (!(is.numeric(overlap) && length(overlap) == 1 &&
        overlap >= 0 && overlap < 1)) {
    stop("overlap must be a numeric value between 0 and 1")
  }
  if (!(is.numeric(tbw) && length(tbw) == 1 && tbw > 0)) {
    stop("tbw must be a positive numeric value")
  }
  k <- round(k)
  if (!(is.numeric(k) && length(k) == 1 && k > 0)) {
    stop("k must be a positive numeric value")
  }
  if (!(is.numeric(p) && length(p) == 1 && tau > 0)) {
    stop("p must be a positive numeric value")
  }
  if (!(is.numeric(tau) && length(tau) == 1 && tau > 0)) {
    stop("tau must be a positive numeric value")
  }
  if (!(is.numeric(maxIter) && length(maxIter) == 1 && maxIter > 0)) {
    stop("k must be a positive numeric value")
  }
  
  # global trend removal
  if (detrend == "long-mean") {
    x <- gsignal::detrend(x, p = 0)
  } else if (detrend == "long-linear") {
    x <- gsignal::detrend(x, p = 1)
  }
  
  # compute and scale tapers 
  n <- round(w * fs)
  n <- 2^ceiling(log2(n))
  nw = tbw * w / 2
  tapers <- multitaper::dpss(n, k, nw)$v * sqrt(fs)
  
  y <- apply(x, 2, cl_chan, fs, freq, fbw, detrend, n,
             overlap, tapers, p, tau, maxIter)
  if (vec) {
    y <- as.vector(y)
  }
  y
}

# define function to do one channel
cl_chan <- function(chan, fs, freq, fbw, detrend, n,
                    overlap, tapers, p, tau, maxIter) {
  
  # initializations
  len <- length(chan)
  nstep <- n - round(n * overlap)
  noverlap <- n - nstep
  smooth = 1 / (1 + exp(-tau *(ifelse(noverlap > 0, 1:noverlap, 0) - noverlap / 2) / noverlap))
  
  sseg <- seq(1, len - n + 1, nstep)
  nw <- length(sseg)
  datafit <- rep(0, sseg[nw] + n - 1)
  
  # initial mt spectrum
  initial <- mtchan(chan, fs, detrend, n, overlap, tapers)
  psd_len <- n / 2 + 1  #n is always even
  initial <- 10 * log10(
    initial[1:psd_len] + c(0, initial[seq(n, psd_len + 1, -1)], 0))
  f <- seq.int(0, psd_len - 1) * (fs / n)
  
  # find indices of frequencies to remove
  fidx <- rep(0, length(freq))
  for (fk in seq_along(freq)) {
    fidx[fk] <- which.min(abs(f - freq[fk]))
  }
  
  # iterate
  for (it in seq_len(maxIter)) {
    freq_mask <- rep(FALSE, length(freq))
    for (iseg in seq_along(sseg)) {
      bseg <- sseg[iseg]
      eseg <- min(len, bseg + n - 1)
      seg <- chan[bseg:eseg]
      if (detrend == "short-mean") {
        seg <- gsignal::detrend(chan[bseg:eseg], p = 0)
      } else if (detrend == "short-linear") {
        seg <- gsignal::detrend(chan[bseg:eseg], p = 1)
      } else {
        seg <- chan[sseg:eseg]
      }
      ff <- fitfreqs(seg, fs, freq, fbw, n, tapers, p)
      freq_mask <- freq_mask | ff$fSig
      if (iseg > 1 && noverlap > 0) {
        fitted[1:noverlap] <- smooth * datafit[1:noverlap] +
          (1 - smooth) * datafit0[(n - noverlap + 1):n]
      } else {
        fitted <- ff$fitted
      }
      datafit[bseg:eseg] <- fitted
      datafit0 <- datafit 
    }
    
    chan[1:length(datafit)] <- chan[1:length(datafit)] - datafit
    if (sum(freq_mask) > 0) {
      # Now find the line frequencies that have converged
      cleaned <- mtchan(chan, fs, detrend, n, overlap, tapers)
      cleaned <- 10 * log10(
        cleaned[1:psd_len] + c(0, cleaned[seq(n, psd_len + 1, -1)], 0))
        
      dBReduction <- initial - cleaned
      idx <- dBReduction[fidx] < 0
      freq <- freq[-which(idx | !freq_mask)]
      fidx <- fidx[-which(idx | !freq_mask)]
      initial <- cleaned
    }
    if (length(freq) <= 0) {
        break;
    } 
  }
  chan
}

# fit frequencies to one segment
fitfreqs <- function(seg, fs, freq, fbw, n, tapers, p) {
  tf <- testfreqs(seg, fs, freq, fbw, n, tapers, p)
  fit <- rep(0, n)
  f <- seq(0, fs / 2, length.out = length(tf$f))
  f_mask <- rep(FALSE, length(f))
  freq_mask <- rep(FALSE, length(freq))
  if (length(fbw) > 0) {
    # scan fbw around freq for largest significant peak of Fval
    for (fi in seq_along(freq)) {
      idx1 <- which.min(abs(f - (freq[fi] - fbw/2)))
      idx2 <- which.min(abs(f - (freq[fi] + fbw/2)))
      Fvals <- tf$Fval[idx1:idx2]
      Fvals[Fvals < tf$sig] <- 0
      if (any(Fvals > 0)) {
        maxidx <- which.max(Fvals)
        indx <- idx1 + maxidx - 1
        f_mask[indx] <- TRUE
        freq_mask[fi] <- TRUE
      }
    }
  } else {
    # remove exact freq if significant
    for (fi in seq_along(freq)) {
      indx <- which.min(abs(f - freq[fi]))
      f_mask[indx] <- tf$Fval[indx] >= tf$sig
    }
  }
  
  # Estimate the contribution of any significant freq
  fSig = f[f_mask]
  aSig = tf$A[f_mask]
  FvalSig = tf$Fval[f_mask]
  if (length(fSig) > 0) {
    fit <- Re(exp(1i * 2 * pi * seq(0, n - 1) * fSig / fs) * aSig +
              exp(-1i * 2 * pi * seq(0, n - 1) * fSig / fs) * Conj(aSig))
  }
  list(fitted = fit, fSig = fSig)
}

# test frequencies in one segment
testfreqs <- function(seg, fs, freq, fbw, n, tapers, p) {
  k <- ncol(tapers)
  f <- seq(1, round(n / 2) + 1)
  nf <- length(f)
  kodd <- seq(1, k, 2)
  keven = seq(2, k, 2)
  J <- mtfft(seg, fs, tapers) # tapered fft of data
  Jo = J[f, kodd] # drop the even ffts and restrict frequencies
  H0 <- apply(tapers[, kodd], 2, sum)
  H0sq <- sum(H0 * H0)
  JoH0 = apply(Jo * H0, 1, sum)
  A <- JoH0 / H0sq# % amplitudes for all frequencies and channels
  Jhat <- A %*% t(H0)

  num <- (k - 1) * (abs(A)^2) * H0sq
  den <- apply(abs(Jo - Jhat)^2, 1, sum) + apply(abs(J[f, keven])^2, 1, sum)
  Fval <- num / den
  sig <- stats::qf(1 - p, 2, 2 * k - 2)
  A <- A * fs
  
  list(Fval = Fval, A = A, f = f, sig = sig)
}

