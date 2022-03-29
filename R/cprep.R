# cprep.R
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
# 20220315  GvB           Initial setup v0.3-2
#---------------------------------------------------------------------------------------------------------------------

#' cprep
#' 
#' PREP pipeline for continuous data
#' 
#' @param x input time series, specified as a numeric matrix or vector. In case
#'   of a vector it represents a single signal; in case of a matrix each column
#'   is a signal. Alternatively, an object of class \code{\link{ctd}}
#' @param fs sampling frequency of \code{x} in Hz. Default: 1. Overruled if
#'   \code{x} is a \code{ctd} object, in which case the sampling frequency is
#'   \code{fs(x)}
#' @param dtsens sensors to detrend, specified as positive integers indicating
#'   sensors numbers, or as sensor names (colnames of \code{x}). Default: all
#'   sensors names in \code{x} excluding a variable \code{t} (if present).
#' @param detrend character string indicating which type of trend removal
#'   is performed for line noise removal bad channel detection; one of: 
#'   \describe{
#'     \item{\code{highpass}}{use a high-pass FIR filter with cutoff frequency
#'      specified by the parameter \code{hpcf} (default)}
#'     \item{\code{linear}}{remove polynomial linear trend}
#'     \item{\code{mean|constant}}{remove mean}
#'     \item{\code{none}}{no detrending}
#'   }
#' @param hpcf high-pass cutoff frequency in Hz, specified as a positive numeric
#'   value. Default: 1 Hz
#' @param hptbw high-pass transition bandwith in Hz, specified as a positive
#'   numeric value. For instance, if the high-pass cutoff frequency is 1 Hz, and
#'   the transition bandwith is also 1 Hz, the transition band of 1 Hz is
#'   located between 0.5 and 1.5 Hz. Default: 1 Hz.
#' @param hpdev deviation from desired stop- and passband of the high-pass filter,
#'   specified as a positive numeric value. Default: 0.01
#' @param fbw frequency bandwidth. Bandwidth centered on each \code{lfreq} to
#'   scan for significant lines. Default: 2 Hz.
#' @param lnsens sensors to remove line noise from, specified as positive
#'   integers indicating sensors numbers, or as sensor names (colnames of
#'   \code{x}). Default: idetical to \code{dtsens}.
#' @param lfreq line frequencies to remove, either specified as a single numeric
#'   value, or as a character string. If a single numeric value is specified,
#'   then this value will be expanded to include multiples of that frequency
#'   (harmonics) up to the Nyquist frequency. If a character string is
#'   specified, then the string will be converted to numeric values without
#'   further expansion. Default: 50
#' @param w length of sliding window (segments) in seconds. Default: 4
#' @param overlap proportion of overlap between the segments. Default: 0
#' @param tbw taper bandwidth, specified as a positive numeric value. Default: 2
#'   Hz
#' @param k number of tapers to use, specified as a positive numeric value.
#'   Default: \code{tbw * w - 1}
#' @param p significance level cutoff for F-test. Default 0.01
#' @param tau Window overlap smoothing factor. Default: 100
#' @param MaxIterLn Maximum times to iterate line noise removal. Default: 10
#' @param rrsens vector of sensor names or numbers that are used for the
#'   referencing (hence they must contain EEG data). Default: lnsens.
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
#' @param MaxIterRR maximum number of iterations used for detecting noisy
#'   channels and recomputing the reference. Default: 4.  Only used in case
#'   \code{"newref"} is \code{"robust"}.
#' @param report name of PDF file to which a detailed processing report is
#'   printed, or NULL to suppress the report. Default: "cprep_report.pdf", saved
#'   to the current working directory
#'   
#' @return A list consisting of the following elements:
#'   \describe{
#'     \item{y}{a \code{\link{ctd}} object containing the robust referenced data}
#'     \item{noisyPre}{list containing the output from running
#'       \code{\link{noisysensors}}} on the input data
#'     \item{noisyProst}{list containing the output from running
#'       \code{\link{noisysensors}}} on the output data
#' }
#' 
#' @seealso \code{\link[gsignal]{detrend}}, \code{\link{removetrend}}
#'   \code{\link{cleanline}}, \code{link{noisysensors}},
#'   \code{\link{reref}}
#'   
#' @examples
#' \dontrun{
#'   rr <- cprep(EEGdata[, 1:29], fs = fs(EEGdata))
#' }
#'   
#' @author Original Matlab code by Tim Mullen (2011); ported to R and
#'   adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com},
#'   
#' @export

cprep <- function(x, fs = 1,
                  dtsens = setdiff(colnames(x), "t"),
                  detrend = c("highpass", "linear", "mean", "constant", "none"), 
                  hpcf = 1, hptbw = 1, hpdev = 0.01,
                  lnsens = dtsens,  lfreq = 50, 
                  fbw = 2, w = 4, overlap = 0,
                  tbw = 2, k = round(tbw * w - 1),
                  p = 0.01, tau = 100, MaxIterLn = 10,
                  rrsens = lnsens, oldref = NULL,
                  newref = "robust", saveref = FALSE,
                  interp = c("post-reference", "pre-reference", "none"),
                  estmean = c("median", "huber", "mean"),
                  sl = getLocationsfromLabels(colnames(x[, rrsens])),
                  MaxIterRR = 4,
                  report = paste(getwd(), "cprep_report.pdf", sep = "/")
                  ) {
  
  x_name <- deparse(substitute(x))
  x_name <- gsub("$", paste0("\\", "$"), x_name, fixed = TRUE)   #Latex will choke on $
  
  # parameter checking: x and fs
  if ("ctd" %in% class(x)) {
    fs <- fs(x)
    npts <- npts(x)
  } else {
    if (!(is.numeric(fs) && length(fs) == 1 && fs > 0)) {
      stop("Sampling frequency (fs) must be a positive numeric value")
    }
    if (is.vector(x)) {
      x <- matrix(x, ncol = 1)
      vec <- TRUE
    } else {
      x <- data.matrix(x)
      vec <- FALSE
    }
    npts <- nrow(x)
  }
  ns <- NCOL(x)
  
  #parameter checking: detrend
  detrend <- match.arg(detrend)
  if (detrend == "highpass") {
    if (!(is.numeric(hpcf) && length(hpcf) == 1 && hpcf > 0)) {
      stop("high-pass cutoff frequency (hpcf) must be a positive numeric value")
    }
    if (!(is.numeric(hptbw) && length(hptbw) == 1 && hptbw > 0)) {
      stop("high-pass transition bandwith (hptbw) must be a positive numeric value")
    } 
    if (!(is.numeric(hpdev) && length(hpdev) == 1 && hptbw > 0)) {
      stop("high-pass filter deviation (hpdev) must be a positive numeric value")
    } 
  }
  if (is.numeric(dtsens)) {
    if(any(dtsens <= 0) || max(dtsens) > ncol(x)) {
      stop(paste("Invalid sensor numbers for detrending:", dtsens))
    }
  } else if (is.character(dtsens)) {
    if (!all(dtsens %in% colnames(x))) {
      stop(paste("Invalid sensor names for detrending:", dtsens))
    }
  } else {
    stop("Sensors for detrending must be specified as names or numbers")
  }
  
  # parameter checking: line noise removal
  if (!(is.numeric(lfreq) || is.character(lfreq))) {
      stop("line frequency must be specified as numeric or character")
  }
  if (is.numeric(lnsens)) {
    if(any(lnsens <= 0) || max(lnsens) > ncol(x)) {
      stop(paste("Invalid sensor numbers for line noise removal:", lnsens))
    }
  } else if (is.character(lnsens)) {
    if (!all(lnsens %in% colnames(x))) {
      stop(paste("Invalid sensor names for line noise removal:", lnsens))
    }
  } else {
    stop("Sensors for line noise removal must be specified as names or numbers")
  }
  if (!(is.numeric(fbw) && length(fbw) == 1 && fbw > 0)) {
    stop("fbw must be a positive numeric value")
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
  if (!(is.numeric(MaxIterLn) && length(MaxIterLn) == 1 && MaxIterLn > 0)) {
    stop("maxIterLn must be a positive numeric value")
  }
  
  # parameter checking: referencing
  if (is.numeric(rrsens)) {
    if(any(rrsens <= 0) || max(rrsens) > ncol(x)) {
      stop(paste("Invalid sensor numbers for referencing:", rrsens))
    }
  } else if (is.character(rrsens)) {
    if (!all(rrsens %in% colnames(x))) {
      stop(paste("Invalid sensor names for referencing:", lnsens))
    }
  } else {
    stop("Sensors for referencing must be specified as names or numbers")
  }
  if (!is.null(oldref) && !is.na(oldref) &&
      !is.numeric(oldref) && !is.character(oldref)) {
    stop("Invalid old reference")
  }
  if (is.numeric(newref)) {
    if (any(newref > ncol(x))) {
      stop("invalid numeric new reference sensors specified")
    }
  } else if (is.character(newref)) {
    if (length(newref) == 1 && newref == "robust") {
      interp <- match.arg(interp)
      estmean <- match.arg(estmean)
    }
  } else {
    stop("newref must be numeric or character")
  }
  if (!is.logical(saveref)) {
    saveref <- FALSE
  }
  if (!(is.numeric(MaxIterRR) && length(MaxIterRR) == 1 && MaxIterRR > 0)) {
    stop("maxIterRR must be a positive numeric value")
  }

  # parameter checking: report
  if (is.null(report)) {
    rep <- FALSE
  } else {
    rep <- TRUE
    if (!is.character(report)) {
      stop("Invalid report name")
    }
  }
  
  #############################################################################
  # Set up the processing report
  if (rep) {
    #on.exit(grDevices::graphics.off())
    Rmd <- tempfile(fileext = ".Rmd")
    con <- file(Rmd, "w+")
    writeLines('---', con)
    writeLines('title: "Results from the cPREP pipeline"', con)
    descr <- system.file("DESCRIPTION", package = "eegr")
    writeLines(paste('author: cPREP version eegr',
                     desc::desc_get_version(descr)), con)
    writeLines(paste('date: "report generated `r date()`"'), con)
    writeLines('output:', con)
    writeLines('  pdf_document:', con)
    writeLines('    toc: true', con)
    writeLines('    latex_engine: xelatex', con)
    writeLines('header-includes:', con)
    writeLines('- \\usepackage{booktabs}', con)
    writeLines('- \\usepackage{longtable}', con)
    writeLines('- \\usepackage{array}', con)
    writeLines('- \\usepackage{multirow}', con)
    writeLines('- \\usepackage{wrapfig}', con)
    writeLines('- \\usepackage{float}', con)
    writeLines('- \\usepackage{colortbl}', con)
    writeLines('- \\usepackage{pdflscape}', con)
    writeLines('- \\usepackage{tabu}', con)
    writeLines('- \\usepackage{threeparttable}', con)
    writeLines('- \\usepackage{threeparttablex}', con)
    writeLines('- \\usepackage[normalem]{ulem}', con)
    writeLines('- \\usepackage{makecell}', con)
    writeLines('- \\usepackage{xcolor}', con)
    writeLines('- \\usepackage{adjustbox}', con)
    writeLines('- \\usepackage{calc}', con)
    writeLines('- \\usepackage{caption}', con)
    writeLines('- \\usepackage{fontspec}', con)
    writeLines('- \\usepackage{graphicx}', con)
    writeLines('- \\usepackage{hhline}', con)
    writeLines('- \\usepackage{hyperref}', con)
    writeLines('- \\usepackage{siunitx}', con)
    writeLines('- \\usepackage{tabularx}', con)
    writeLines('- \\usepackage{threeparttable}', con)
    writeLines('- \\usepackage{ulem}', con)
    writeLines('- \\usepackage{caption}', con)
    writeLines('---', con)
    writeLines('', con)

    writeLines('\\captionsetup[table]{labelformat=empty}', con)
    writeLines('', con)
    
    writeLines('```{r setup, include=FALSE}', con)
    writeLines('knitr::opts_chunk$set(echo = FALSE, results = \'asis\')', con)
    writeLines('library(knitr, quietly = TRUE, warn.conflicts = FALSE)', con)
    writeLines('library(huxtable, quietly = TRUE, warn.conflicts = FALSE)', con)
    writeLines('```', con)
    writeLines('', con)
    
    writeLines('\\ ', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
    
    writeLines('# General Information', con)
    writeLines('', con)
    
    writeLines('This document is a summary report of the steps executed by the cPREP pipeline. ', con)
    writeLines('', con)

    col1 <- c(
      "Sampling frequency (Hz):",
      "Number of data points:",
      "Time (seconds):",
      "Total number of signals:",
      "Signal names:"
    )
    col2 <- c(
      fs,
      npts,
      npts / fs,
      ncol(x),
      paste(colnames(x), collapse = ', ')
    )
    table <- huxtable::hux(col1, col2)
    tabke <- huxtable::set_valign(table, huxtable::final(1), huxtable::everywhere, "top")
    huxtable::width(table) <- 0.9
    huxtable::wrap(table) <- TRUE
    huxtable::caption(table) <- paste('Characteristics of', x_name)
    huxtable::top_border(table)[1, ] <- 1
    huxtable::bottom_border(table)[nrow(table), ] <- 1
    
    writeLines(huxtable::to_latex(table), con)
    writeLines('', con)

    col1 <- c(
      "x:", "fs:", "dtsens:", "detrend:", "hpcf:", "hptbw:", "hpdev:",
      "lnsens:",  "lfreq:", "fbw:", "w:", "overlap:", "tbw:", "k:" ,
      "p:", "tau:", "MaxIterLn:", "rrsens:", "oldref:", "newref:", "saveref:",
      "interp:", "estmean:", "sl:",  "MaxIterRR:", "report:"
    )
    col2 <- c(x_name, fs, paste(dtsens, collapse = ', '), detrend, hpcf, hptbw, hpdev,
              paste(lnsens, collapse = ', '),  lfreq, fbw, w, overlap, tbw, k,
              p, tau, MaxIterLn, paste(rrsens, collapse = ', '), 
              paste(oldref, collapse = ', '), paste(newref, collapse = ', '), saveref,
              interp, estmean, deparse(substitute(sl)), MaxIterRR, paste(report, collapse = ", ")
    )
    table <- huxtable::hux(col1, col2)
    table <- huxtable::set_valign(table, huxtable::final(1), huxtable::everywhere, "top")
    huxtable::width(table) <- 0.9
    huxtable::wrap(table) <- TRUE
    huxtable::caption(table) <- "Arguments to cprep() function"
    huxtable::top_border(table)[1, ] <- 1
    huxtable::bottom_border(table)[nrow(table), ] <- 1
    
    writeLines(huxtable::to_latex(table), con)
    writeLines('', con)
    
    writeLines('The following pages provide details about the actions performed in each step of', con)
    writeLines('the cPREP pipeline.  ', con)
    writeLines('', con)

    writeLines('\\ ', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
    
  }
  
  
  #############################################################################
  # Phase 1: detrend the data
  tic <- Sys.time()
  
  if (rep) {
    writeLines('# Phase 1: Initial detrending', con)
    writeLines('', con)
    writeLines(paste0('Detrending requested: ', detrend, '.'), con)
    if (detrend == "highpass") {
      writeLines(paste('The requested cutoff frequency is', hpcf, 'Hz, '), con)
      writeLines(paste('with a transition bandwidth of', hptbw, 'Hz, '), con)
      writeLines(paste0('and a deviation of ', hpdev, '(', format(hpdev/100, scientific = FALSE), '%)'), con)
      writeLines(paste('both in the stopband and the passband.'), con)
      writeLines('', con)
      writeLines(paste0('Detrend requested on sensors ', paste(dtsens, collapse = ', '), '.'), con)
    } else if (detrend == "none") {
      writeLines(('. Note that is is not recommended.'), con)
      writeLines(('. The procedures work better with initial detrending..'), con)
    } else {
      writeLines(('.'), con)
    }
    writeLines('', con)
  }
  
  result <- try({
    detr <- removetrend(x[, dtsens], fs, detrend, hpcf, hptbw, hpdev)
  }, silent = TRUE)
  
  if (!"try-err" %in% class(result)) {
    x_filt <- detr$y
    if (rep) {
      if (detrend == "highpass") {
        hp <- detr$h
        l <- length(hp)
        ord <- l - 1
        delay <- round(ord / 2)
        writeLines('Here is the frequency response of the high-pass filter', con)
        writeLines(paste0('with filter order ', ord, ' (length ', l, ').'), con)
        writeLines('The filter has a linear phase.', con)
        writeLines(paste0('Its delay is ', delay, ' samples, corresponding to ', delay/fs, ' seconds.'), con)
        writeLines('The procedure will correct for this this delay.', con)
        writeLines('', con)
        png_detr_freqz <- tempfile(pattern = "img", fileext = ".png")
        grDevices::png(png_detr_freqz)
        print(gsignal::freqz(hp, fs = fs))
        rc <- grDevices::dev.off()
        writeLines('```{r, fig.align=\'center\'}', con)
        writeLines(paste0('knitr::include_graphics("', png_detr_freqz, '")'), con)
        writeLines(paste0('knitr::include_graphics("', png_detr_freqz, '")'), con)
        writeLines('```', con)
        writeLines('', con)
    
        writeLines('\\ ', con)
        writeLines('\\clearpage', con)
        writeLines('', con)
      }
      
      png_detr_before <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_detr_before, width = 14, height = 10, units = "cm", res = 600)
      plot(gsignal::pwelch(x[, dtsens], window = 4 * fs, overlap = 0.75, fs = fs, detrend = "none"),
           yscale = "dB", main = paste("Welch spectrum before detrending\n",
                                       "object:", x_name))
      rc <- grDevices::dev.off()
      png_detr_after <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_detr_after, width = 14, height = 10, units = "cm", res = 600)
      plot(gsignal::pwelch(x_filt[, dtsens], window = 4 * fs, overlap = 0.75, fs = fs, detrend = "none"),
           yscale = "dB", main = paste("Welch spectrum after detrending\n",
                                       "object:", x_name))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_detr_before, '")'), con)
      writeLines(paste0('knitr::include_graphics("', png_detr_after, '")'), con)
      writeLines('```', con)
      writeLines('', con)
    }
  } else {
    warning("Filtering unsuccesful. Error message:")
    warning(result)
    if (rep) {
      writeLines('\\', con)
      writeLines('Filtering unsuccesful. Error message:', con)
      writeLines(paste(result), con)
      writeLines('', con)
    }
    x_filt <- x
  }

  toc <- Sys.time()
  tdiff <- toc - tic
  if (rep) {
    writeLines(paste('Elapsed time:', format(tdiff, digits = 4, scientific = FALSE)), con)
    writeLines('\\ ', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
  }
  
  #############################################################################
  # Phase 2: remove line noise
  tic <- Sys.time()
  
  if (rep) {
    writeLines('# Phase 2: Line noise removal', con)
    writeLines('', con)
    writeLines('The following table lists the parameters for the line noise detrending procedure.', con)
    writeLines('', con)
    
    col1 <- c(
      "Signal names:",
      "Line noise frequency (Hz):",
      "Line frequency bandwidth (Hz):",
      "Taper window width (seconds):",
      "Overlap (seconds):",
      "Taper bandwidth (Hz):",
      "Number of tapers:",
      "Significance level (p):",
      "Smoothing factor (tau):",
      "Maximum iterations:"
    )
    col2 <- c(
      paste(lnsens, collapse = ', '),
      lfreq, 
      fbw,
      w,
      overlap,
      tbw,
      k,
      p,
      tau,
      MaxIterLn
    )
    table <- huxtable::hux(col1, col2)
    huxtable::valign(table) <- "top"
    huxtable::width(table) <- 0.9
    huxtable::wrap(table) <- TRUE
    huxtable::caption(table) <- paste('Line noise removal parameters for', x_name)
    huxtable::top_border(table)[1, ] <- 1
    huxtable::bottom_border(table)[nrow(table), ] <- 1
    
    writeLines(huxtable::to_latex(table), con)
    writeLines('', con)
  }
  
  result <- try(
    suppressWarnings(
      cleaned <- cleanline(x_filt[, lnsens], fs, lfreq, fbw, detrend = "none",
                           w, overlap, tbw, k, p, tau, MaxIterLn)
    ), silent = TRUE)

  if (!"try-err" %in% class(result)) {
    if (rep) {
      l <- length(cleaned$lfreq)
      if (l == 1) {
        str <- paste(cleaned$lfreq)
      } else if (l == 2) {
        str <- paste(cleaned$lfreq, collapse = " and ")
      } else {
        str <- paste0(paste(cleaned$lfreq[1:(l - 1)], collapse = ", "), ", and ", cleaned$lfreq[l])
      }
      writeLines(paste('The line noise frequencies that the procedure operated on were:', str, "Hz."), con)
      writeLines('', con)

      writeLines('The following table shows the number of iterations performed for each sensor.', con)
      writeLines('', con)
      
      df <- data.frame(colnames(cleaned$y), cleaned$niter)
      table <- huxtable::hux(t(df))
      huxtable::valign(table) <- "top"
      huxtable::width(table) <- 0.9
      huxtable::wrap(table) <- TRUE
      huxtable::caption(table) <- 'Number of iterations'
      huxtable::bottom_border(table)[1, ] <- 1
      huxtable::bottom_border(table)[nrow(table), ] <- 1
      if (ns > 10) {
        ltable <- huxtable::split_down(table, seq(10, ns, 10))
      } else {
        ltable <- list(table)
      }
      for (i in 1:length(ltable)) {
        huxtable::latex_float(ltable[[i]]) <- "h!"
        if (i > 1) huxtable::caption(ltable[[i]]) <- NA
        writeLines(huxtable::to_latex(ltable[[i]]), con)
        writeLines('', con)
      }
      writeLines('', con)
      
      # Welch spectra pre and post
      bspec <- gsignal::pwelch(x_filt[, lnsens], window = 4 * fs, overlap = 0.75, fs = fs, detrend = "none")
      espec <- gsignal::pwelch(cleaned$y, window = 4 * fs, overlap = 0.75, fs = fs, detrend = "none")
      
      #get noise reduction for the requested frequencies
      nfrex <- length(cleaned$lfreq)
      fidx <- rep(0, nfrex)
      for (fk in seq_len(nfrex)) {
        fidx[fk] <- which.min(abs(espec$freq - cleaned$lfreq[fk]))
        writeLines('', con)
      }
      nr <- matrix(10 * log10(espec$spec[fidx, ] / bspec$spec[fidx, ]), nrow = nfrex)
      rownames(nr) <- paste(cleaned$lfreq, "Hz")

      writeLines('And the following table shows the noise reduction for each sensor at each line frequency.', con)
      writeLines('', con)

      ns <- NCOL(cleaned$y)
      df <- rbind(colnames(cleaned$y), nr)      
      table <- huxtable::hux(df)
      huxtable::number_format(table) <- 2
      huxtable::number_format(table)[1, ] <- NA
      table <- huxtable::add_rownames(table)
      huxtable::valign(table) <- "top"
      huxtable::width(table) <- 0.9
      huxtable::wrap(table) <- TRUE
      huxtable::caption(table) <- 'Line noise reduction in dB'
      huxtable::bottom_border(table)[1, ] <- 1
      huxtable::bottom_border(table)[nrow(table), ] <- 1
      if (ns > 10) {
        ltable <- huxtable::split_down(table, seq(10, ns, 10))
      } else {
        ltable <- list(table)
      }
      for (i in 1:length(ltable)) {
        huxtable::latex_float(ltable[[i]]) <- "h!"
        if (i > 1) huxtable::caption(ltable[[i]]) <- NA
        writeLines(huxtable::to_latex(ltable[[i]]), con)
        writeLines('', con)
      }

      writeLines('\\ ', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      png_lnrem_before <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_lnrem_before, width = 14, height = 10, units = "cm", res = 600)
      plot(bspec, yscale = "dB", 
           main = paste("Welch spectrum before line noise removal\n",
                        "object:", x_name))
      rc <- grDevices::dev.off()
      png_lnrem_after <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_lnrem_after, width = 14, height = 10, units = "cm", res = 600)
      plot(espec, yscale = "dB", 
           main = paste("Welch spectrum after line noise removal\n",
                        "object:", x_name))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_lnrem_before, '")'), con)
      writeLines(paste0('knitr::include_graphics("', png_lnrem_after, '")'), con)
      writeLines('```', con)
      writeLines('', con)
    }
  } else {
    warning("Line noise removal unsuccesful. Error message:")
    warning(result)
    if (rep) {
      writeLines('\\', con)
      writeLines('Line noise removal unsuccesful. Error message:', con)
      writeLines(paste(result), con)
    }
    x_lnrem <- x_filt
  }
  
  toc <- Sys.time()
  tdiff <- toc - tic
  if (rep) {
    writeLines('', con)
    writeLines(paste('Elapsed time:', format(tdiff, digits = 4, scientific = FALSE)), con)
    writeLines('\\ ', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
  }

  # undo the effect of filtering
  # invoke garbage collection (useful for very large files)
  
  x <- x[, lnsens] - x_filt + cleaned$y
  remove(x_filt, cleaned, result)
  rc <- gc()
  
  #############################################################################
  # Phase 3: robust reference
  
  tic <- Sys.time()
  
  if (rep) {
    writeLines('# Phase 3: Robust Reference', con)
    writeLines('', con)
    writeLines('The following table lists the parameters for referencing procedure.', con)
    writeLines('', con)
    
    col1 <- c(
      "Signal names:",
      "Old reference:",
      "New reference:",
      "Save reference:",
      "Interpolation order:",
      "Initial estimate for mean:",
      "Maximum iterations:"
    )
    col2 <- c(
      paste(rrsens, collapse = ', '),
      deparse(oldref),
      newref,
      saveref,
      interp,
      estmean,
      MaxIterRR
    )
    table <- huxtable::hux(col1, col2)
    huxtable::valign(table) <- "top"
    huxtable::width(table) <- 0.9
    huxtable::wrap(table) <- TRUE
    huxtable::caption(table) <- paste('Referencing parameters for', x_name)
    huxtable::top_border(table)[1, ] <- 1
    huxtable::bottom_border(table)[nrow(table), ] <- 1
    
    writeLines(huxtable::to_latex(table), con)
    writeLines('', con)
    
    writeLines('', con)
    writeLines('\\ ', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
    
  }
  
  result <- try(
    suppressWarnings({
      noisyPre <- noisysensors(
        removetrend(x, fs, detrend, hpcf, hptbw, hpdev)$y,
        fs = fs)
      rr <- reref(x, fs, rrsens, oldref, newref, saveref,
                  interp, estmean, sl, MaxIterRR)
      noisyPost <- noisysensors(
        removetrend(rr$y, fs, detrend, hpcf, hptbw, hpdev)$y,
        fs = fs)
    }), silent = TRUE)
  
  if (!"try-err" %in% class(result)) {
    # Function to write bad sensors to string (or "None")
    wrs <- function(s) {
      if (length(s) > 0) paste(s, collapse = ' ') else 'None'  
    }
    if (rep) {
      writeLines('## Bad sensors', con)
      writeLines('', con)
      writeLines('### Before interpolation and referencing', con)
      writeLines('', con)
      writeLines(paste(ifelse(length(noisyPre$TotalBad) == 1, '1 sensor was', 
                              paste(length(noisyPre$TotalBad), 'sensors were'))
                       ,'classified as bad:',
                       wrs(noisyPre$BadNames)), con)
      writeLines('', con)
      writeLines(paste('Bad because of missing values:',
                       wrs(colnames(x)[noisyPre$BadFromNA])), con)
      writeLines('', con)
      writeLines(paste('Bad because data is constant:',
                       wrs(colnames(x)[noisyPre$BadFromConstant])), con)
      writeLines('', con)
      writeLines(paste('Bad because of dropouts:', 
                       wrs(colnames(x)[noisyPre$BadFromDropout])), con)
      writeLines('', con)
      writeLines(paste('Bad because of large deviation:',
                       wrs(colnames(x)[noisyPre$BadFromDeviation])), con)
      writeLines('', con)
      writeLines(paste('Bad because of HF noise:',
                       wrs(colnames(x)[noisyPre$BadFromHFnoise])), con)
      writeLines('', con)
      writeLines(paste('Bad because of poor max correlation:',
                       wrs(colnames(x)[noisyPre$BadFromCorrelation])), con)
      writeLines('', con)
      writeLines(paste('Bad because of poor RANSAC predictability:',
                       wrs(colnames(x)[noisyPre$BadFromRansac])), con)
      writeLines('', con)
      
      writeLines('### After interpolation and referencing', con)
      writeLines('', con)
      writeLines(paste(ifelse(length(noisyPost$TotalBad) == 1, '1 sensor was', 
                              paste(length(noisyPost$TotalBad), 'sensors were'))
                       ,'classified as bad:',
                       wrs(noisyPost$BadNames)), con)
      writeLines('', con)
      writeLines(paste('Bad because of missing values:',
                       wrs(colnames(x)[noisyPost$BadFromNA])), con)
      writeLines('', con)
      writeLines(paste('Bad because data is constant:',
                       wrs(colnames(x)[noisyPost$BadFromConstant])), con)
      writeLines('', con)
      writeLines(paste('Bad because of dropouts:',
                       wrs(colnames(x)[noisyPost$BadFromDropout])), con)
      writeLines('', con)
      writeLines(paste('Bad because of large deviation:',
                       wrs(colnames(x)[noisyPost$BadFromDeviation])), con)
      writeLines('', con)
      writeLines(paste('Bad because of HF noise:',
                       wrs(colnames(x)[noisyPost$BadFromHFnoise])), con)
      writeLines('', con)
      writeLines(paste('Bad because of poor max correlation:',
                       wrs(colnames(x)[noisyPost$BadFromCorrelation])), con)
      writeLines('', con)
      writeLines(paste('Bad because of poor RANSAC predictability:',
                       wrs(colnames(x)[noisyPost$BadFromRansac])), con)
      writeLines('', con)
      
      # function to plot sensor labels padded with codes
      plot_sens <- function(sl, noisy) {
        ss <- sl[sl$label %in% noisy$BadNames, ]
        if (nrow(ss) > 0) {
          # graphics::points(ss$x2d, ss$y2d, xlim = c(-2, 2), ylim = c(-2, 2),
          #                bty = "n", pch = 20)
          ss$label <- as.character(ss$label)
          for (i in 1:nrow(ss)) {
            if (ss$label[i] %in% colnames(x)[noisy$BadFromNA])
              ss$label[i] <- paste0(ss$label[i], 'n')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromConstant])
              ss$label[i] <- paste0(ss$label[i], 'z')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromDropout])
              ss$label[i] <- paste0(ss$label[i], 'd')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromDeviation])
              ss$label[i] <- paste0(ss$label[i], '+')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromHFnoise])
              ss$label[i] <- paste0(ss$label[i], 'x')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromCorrelation])
              ss$label[i] <- paste0(ss$label[i], 'c')
            if (ss$label[i] %in% colnames(x)[noisy$BadFromRansac])
              ss$label[i] <- paste0(ss$label[i], '?')
            graphics::text(x = ss$x2d[i], y = (ss$y2d[i] + 0.1), labels = ss$label[i])
          }
        }
      }
      
      writeLines('\\ ', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      writeLines('## Deviation criterion: robust sensor deviation', con)
      
      png_dev_before <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_dev_before, width = 16, height = 14, units = "cm", res = 600)
      topoplot(noisyPre$Deviation, sl, plot = c("sensorlocs", "legend"), main = "Before referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPre)
      rc <- grDevices::dev.off()
      png_dev_after <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_dev_after, width = 16, height = 14, units = "cm", res = 600)
      topoplot(noisyPost$Deviation, sl, plot = c("sensorlocs", "legend"), main = "After referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPost)
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_dev_before, '")'), con)
      writeLines(paste0('knitr::include_graphics("', png_dev_after, '")'), con)
      writeLines('```', con)
      writeLines('', con)

      writeLines('\\ ', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      # Windowed deviation (1 sec windows)
      winPre <- with(noisyPre, (cor_dev - devMed) / devSD)
      nsPre <- apply(winPre > noisyPre$devthr, 1, sum)
      winPost <- with(noisyPost, (cor_dev - devMed) / devSD)
      nsPost <- apply(winPost > noisyPost$devthr, 1, sum)
      
      # Deviation statistics table      
      col1 <- c(
        "Deviation threshold:",
        "Bad sensors because of deviation above this threshold:",
        "Number of segments for which deviation is calculated:",
        paste("Mean number of sensors with deviation above",
              "threshold in 1-second segments:"),
        paste("Mean proportion of sensors with deviation above",
              "threshold in 1-second segments:"),
        "Maximum raw deviation level:",
        paste("Number of windows with more than 25% of",
              "sensors above threshold:"),
        paste("Number of windows with more than 50% of",
              "sensors above threshold:"),
        "median window deviations",
        "SD window deviations:"
      )
      col2 <- c(
        noisyPre$devthr,
        wrs(colnames(x)[noisyPre$BadFromDeviation]),
        nrow(noisyPre$cor_dev),
        mnsPre <- mean(nsPre),
        mnsPre / ns,
        max(noisyPre$dev),
        sum(nsPre > round(0.25 * ns)),
        sum(nsPre > round(0.50 * ns)),
        stats::median(noisyPre$cor_dev),
        stats::sd(noisyPre$cor_dev)
      )
      col3 <- c(
        noisyPre$devthr,
        wrs(colnames(x)[noisyPost$BadFromDeviation]),
        nrow(noisyPre$cor_dev),
        mnsPost <- mean(nsPost),
        mnsPost / ns,
        max(noisyPost$dev),
        sum(nsPost > round(0.25 * ns)),
        sum(nsPost > round(0.50 * ns)),
        stats::median(noisyPost$cor_dev),
        stats::sd(noisyPost$cor_dev)
      )
      df <- data.frame(col1, col2, col3)
      colnames(df) <- c(" ", "Before", "After")
      table <- huxtable::hux(df)
      huxtable::valign(table) <- "bottom"
      huxtable::width(table) <- 0.9
      huxtable::wrap(table) <- TRUE
      huxtable::caption(table) <- paste('Robust deviation window statistics', x_name)
      table <- huxtable::set_bottom_border(table, huxtable::brdr(1, "solid", "black"))
      
      writeLines(huxtable::to_latex(table), con)
      writeLines('', con)

      writeLines('\\ ', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      # Cumulative distribution function of deviation scores
      png_cdfdev <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_cdfdev, width = 14, height = 10, units = "cm", res = 600)
      plot(stats::ecdf(winPre), xlim = c(-5, 5), pch = NA,
           xlab = "Deviation Score", ylab = "Cumulative probability",
           main = "Deviation score distribution")
      graphics::lines(stats::ecdf(winPost), col = "red", pch = NA)
      graphics::legend("left", legend = c("Before", "After"), lty = 1, col = c("black", "red"))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_cdfdev, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      # sensors above deviation threshold, by interval
      png_devsens <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_devsens, width = 14, height = 10, units = "cm", res = 600)
      mx <- ceiling(max(nsPre, nsPost))
      plot(nsPre, pch = "+", ylim = c(0, mx), 
           xlab = "Intervals", ylab = paste0("# sensors (out of ", ns, ")"),
           main = paste("Number of sensors above deviation threshold\n",
                        "(by intervals of", noisyPre$corsecs, "s)"))
      graphics::points(nsPost, pch = "o", col = "red")
      graphics::legend("topleft", legend = c("Before", "After"), pch = c("+", "o"),
             col = c("black", "red"))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_devsens, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      writeLines('', con)
      writeLines('\\ ', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      ### Median max correlation
      writeLines('## Correlation criterion: median maximum absolute correlation', con)
      
      png_medcor_before <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_medcor_before, width = 16, height = 14, units = "cm", res = 600)
      topoplot(noisyPre$corMed, sl, plot = c("sensorlocs", "legend"), main = "Before referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, 1))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPre)
      rc <- grDevices::dev.off()
      png_medcor_after <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_medcor_after, width = 16, height = 14, units = "cm", res = 600)
      topoplot(noisyPost$corMed, sl, plot = c("sensorlocs", "legend"), main = "After referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, 1))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPost)
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_medcor_before, '")'), con)
      writeLines(paste0('knitr::include_graphics("', png_medcor_after, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      # Bad min max correlation fraction
      thrcorr_before <- noisyPre$cors < noisyPre$corthr
      mmcorr_before <- apply(thrcorr_before, 2, mean)
      thrcorr_after <- noisyPost$cors < noisyPost$corthr
      mmcorr_after <- apply(thrcorr_after, 2, mean)
      
      mx <- max(mmcorr_before, mmcorr_after)

      writeLines('## Correlation criterion: proportion maximum absolute correlation below threshold', con)
      
      png_mmcor_before <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_mmcor_before, width = 16, height = 16, units = "cm", res = 600)
      topoplot(mmcorr_before, sl, plot = c("sensorlocs", "legend"), main = "Before referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, mx))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPre)
      rc <- grDevices::dev.off()
      png_mmcor_after <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_mmcor_after, width = 16, height = 14, units = "cm", res = 600)
      topoplot(mmcorr_after, sl, plot = c("sensorlocs", "legend"), main = "After referencing",
               col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, mx))
      graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
      plot_sens(sl, noisyPost)
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_mmcor_before, '")'), con)
      writeLines(paste0('knitr::include_graphics("', png_mmcor_after, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      # Windowed correlation (1 sec windows)
      nsPre <- apply(thrcorr_before, 1, sum)
      nsPost <- apply(thrcorr_after, 1, sum)
      
      # Correlation statistics table      
      col1 <- c(
        "Correlation threshold:",
        "Bad sensors because of correlation below this threshold:",
        "Number of segments for which correlations are calculated:",
        paste("Mean number of sensors with correlation below",
              "threshold in 1-second segments:"),
        paste("Mean proportion of sensors with correlation above",
              "threshold in 1-second segments:"),
        "Minimum raw correlation:",
        "Maximum raw correlation:",
        paste("Number of windows with more than 25% of",
              "sensors above threshold:"),
        paste("Number of windows with more than 50% of",
              "sensors above threshold:")
      )
      col2 <- c(
        noisyPre$corthr,
        wrs(colnames(x)[noisyPre$BadFromCorrelation]),
        nrow(noisyPre$cors),
        mnsPre <- mean(thrcorr_before),
        mnsPre / ns,
        min(noisyPre$cors),
        max(noisyPre$cors),
        sum(nsPre > round(0.25 * ns)),
        sum(nsPre > round(0.50 * ns))
      )
      col3 <- c(
        noisyPost$corthr,
        wrs(colnames(x)[noisyPost$BadFromCorrelation]),
        nrow(noisyPost$cors),
        mnsPost <- mean(thrcorr_after),
        mnsPost / ns,
        min(noisyPost$cors),
        max(noisyPost$cors),
        sum(nsPost > round(0.25 * ns)),
        sum(nsPost > round(0.50 * ns))
      )
      df <- data.frame(col1, col2, col3)
      colnames(df) <- c(" ", "Before", "After")
      table <- huxtable::hux(df)
      huxtable::valign(table) <- "bottom"
      huxtable::width(table) <- 0.9
      huxtable::wrap(table) <- TRUE
      huxtable::caption(table) <- paste('Correlation window statistics', x_name)
      table <- huxtable::set_bottom_border(table, huxtable::brdr(1, "solid", "black"))
      
      writeLines(huxtable::to_latex(table), con)
      writeLines('', con)
      
      writeLines('\\', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      # Cumulative distribution function of maximum correlation scores
      png_cdfcorr <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_cdfcorr, width = 14, height = 10, units = "cm", res = 600)
      plot(stats::ecdf(noisyPre$cors), xlim = c(0, 1), pch = NA,
           xlab = "Correlation", ylab = "Cumulative probability",
           main = "Maxiumum correlation distribution")
      graphics::lines(stats::ecdf(noisyPost$cors), col = "red", pch = NA)
      graphics::legend("left", legend = c("Before", "After"), lty = 1, col = c("black", "red"))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_cdfcorr, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      # sensors with correlation below threshold, by interval
      png_corrsens <- tempfile(pattern = "img", fileext = ".png")
      grDevices::png(png_corrsens, width = 14, height = 10, units = "cm", res = 600)
      mx <- ceiling(max(nsPre, nsPost))
      plot(nsPre, pch = "+", ylim = c(0, mx), 
           xlab = "Intervals", ylab = paste0("# sensors (out of ", ns, ")"),
           main = paste("Number of sensors below correlation threshold\n",
                        "(by intervals of", noisyPre$corsecs, "s)"))
      graphics::points(nsPost, pch = "o", col = "red")
      graphics::legend("topleft", legend = c("Before", "After"), pch = c("+", "o"),
                       col = c("black", "red"))
      rc <- grDevices::dev.off()
      writeLines('```{r, fig.align=\'center\'}', con)
      writeLines(paste0('knitr::include_graphics("', png_corrsens, '")'), con)
      writeLines('```', con)
      writeLines('', con)
      
      writeLines('', con)
      writeLines('\\clearpage', con)
      writeLines('', con)
      
      # bad ransac fraction
      writeLines('## Predictability criterion: RANSAC sensor fraction', con)

      if (noisyPre$ransacPerformed && noisyPost$ransacPerformed) {
      
        png_rsf_before <- tempfile(pattern = "img", fileext = ".png")
        grDevices::png(png_rsf_before, width = 16, height = 14, units = "cm", res = 600)
        topoplot(noisyPre$ransac_frac, sl[noisyPre$ransacSensors, ],
                 plot = c("sensorlocs", "legend"), main = "Before referencing",
                 col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, 1))
        graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
        plot_sens(sl, noisyPre)
        rc <- grDevices::dev.off()
        png_rsf_after <- tempfile(pattern = "img", fileext = ".png")
        grDevices::png(png_rsf_after, width = 16, height = 14, units = "cm", res = 600)
        topoplot(noisyPost$ransac_frac, sl[noisyPost$ransacSensors, ],
                 plot = c("sensorlocs", "legend"), main = "After referencing",
                 col = grDevices::colorRampPalette(c("blue", "white", "red"))(50), scale = c(0, 1))
        graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
        plot_sens(sl, noisyPost)
        rc <- grDevices::dev.off()
        writeLines('```{r, fig.align=\'center\'}', con)
        writeLines(paste0('knitr::include_graphics("', png_rsf_before, '")'), con)
        writeLines(paste0('knitr::include_graphics("', png_rsf_after, '")'), con)
        writeLines('```', con)
        writeLines('', con)
      
        # RANSAC statistics table
        thr_before <- noisyPre$rcor < noisyPre$ransacthr
        thr_after <- noisyPost$rcor < noisyPost$ransacthr
        nsPre <- apply(thr_before, 1, sum)
        nsPost <- apply(thr_after, 1, sum)
        
        col1 <- c(
          "RANSAC performed on sensors:",
          "Fraction of sensors used in RANSAC:",
          "RANSAC threshold:",
          "Segment length (seconds):",
          "Number of segments:",
          "Proportion of time a sensor is allowed to have poor RANSAC predictability:",
          "Sample size:",
          "Bad sensors because of low RANSAC predictability:",
          paste("Mean number of sensors with RANSAC predictability",
                "below threshold in segments:"),
          paste("Mean proportion of sensors with RANSAC predictability",
                "below threshold in segments:"),
          "Minimum RANSAC predictability:",
          "Maximum RANSAC predictability:",
          paste("Number of windows with more than 25% of",
                "sensors below threshold:"),
          paste("Number of windows with more than 50% of",
                "sensors below threshold:")
        )
        col2 <- c(
          paste(colnames(x)[noisyPre$ransacSensors], collapse = ", "),
          noisyPre$sfrac,
          noisyPre$ransacthr,
          noisyPre$rwin,
          ncol(noisyPre$rcorr),
          noisyPre$ubtime,
          noisyPre$rss,
          wrs(colnames(x)[noisyPre$BadFromRansac]),
          mnsPre <- mean(thr_before),
          mnsPre / ns,
          min(noisyPre$rcorr),
          max(noisyPre$rcorr),
          sum(nsPre > round(0.25 * ns)),
          sum(nsPre > round(0.50 * ns))
        )
        col3 <- c(
          paste(colnames(x)[noisyPost$ransacSensors], collapse = ", "),
          noisyPost$sfrac,
          noisyPost$ransacthr,
          noisyPost$rwin,
          ncol(noisyPost$rcorr),
          noisyPost$ubtime,
          noisyPost$rss,
          wrs(colnames(x)[noisyPost$BadFromRansac]),
          mnsPost <- mean(thr_after),
          mnsPost / ns,
          min(noisyPost$rcorr),
          max(noisyPost$rcorr),
          sum(nsPost > round(0.25 * ns)),
          sum(nsPost > round(0.50 * ns))
        )
        df <- data.frame(col1, col2, col3)
        colnames(df) <- c(" ", "Before", "After")
        table <- huxtable::hux(df)
        huxtable::valign(table) <- "bottom"
        huxtable::width(table) <- 0.9
        huxtable::wrap(table) <- TRUE
        huxtable::caption(table) <- paste('RANSAC window statistics', x_name)
        table <- huxtable::set_bottom_border(table, huxtable::brdr(1, "solid", "black"))
      
        writeLines(huxtable::to_latex(table), con)
        writeLines('', con)
      
        writeLines('\\', con)
        writeLines('\\clearpage', con)
        writeLines('', con)
      
        # Cumulative distribution function of RANSAC correlation scores
        png_cdfransac <- tempfile(pattern = "img", fileext = ".png")
        grDevices::png(png_cdfransac, width = 14, height = 10, units = "cm", res = 600)
        plot(stats::ecdf(noisyPre$rcorr), xlim = c(0, 1), pch = NA,
             xlab = "Correlation", ylab = "Cumulative probability",
             main = "RANSAC correlation distribution")
        graphics::lines(stats::ecdf(noisyPost$rcorr), col = "red", pch = NA)
        graphics::legend("left", legend = c("Before", "After"), lty = 1, col = c("black", "red"))
        rc <- grDevices::dev.off()
        writeLines('```{r, fig.align=\'center\'}', con)
        writeLines(paste0('knitr::include_graphics("', png_cdfransac, '")'), con)
        writeLines('```', con)
        writeLines('', con)
      
        # sensors with correlation below threshold, by interval
        png_ransacsens <- tempfile(pattern = "img", fileext = ".png")
        grDevices::png(png_ransacsens, width = 14, height = 10, units = "cm", res = 600)
        mx <- ceiling(max(nsPre, nsPost))
        plot(nsPre, pch = "+", ylim = c(0, mx), 
             xlab = "Intervals", ylab = paste0("# sensors (out of ", ns, ")"),
             main = paste("Number of sensors below RANSAC threshold\n",
                        "(by intervals of", noisyPre$rwin, "s)"))
        graphics::points(nsPost, pch = "o", col = "red")
        graphics::legend("topleft", legend = c("Before", "After"), pch = c("+", "o"),
                         col = c("black", "red"))
        rc <- grDevices::dev.off()
        writeLines('```{r, fig.align=\'center\'}', con)
        writeLines(paste0('knitr::include_graphics("', png_ransacsens, '")'), con)
        writeLines('```', con)
        writeLines('', con)
      
        writeLines('', con)
        writeLines('\\clearpage', con)
        writeLines('', con)
        
      } else {
        writeLines('RANSAC not performed both before and after referencing', con)
      }
    }
    
    writeLines('', con)
    writeLines('\\clearpage', con)
    writeLines('', con)
    
    # HF noise
    writeLines('## Noisiness criterion: HF noise z-score', con)

    png_zHF_before <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_zHF_before, width = 16, height = 14, units = "cm", res = 600)
    topoplot(noisyPre$zHFnoise, sl,
             plot = c("sensorlocs", "legend"), main = "Before referencing",
             col = grDevices::colorRampPalette(c("blue", "white", "red"))(50))
    graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
    plot_sens(sl, noisyPre)
    rc <- grDevices::dev.off()
    png_zHF_after <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_zHF_after, width = 16, height = 14, units = "cm", res = 600)
    topoplot(noisyPost$zHFnoise, sl,
             plot = c("sensorlocs", "legend"), main = "After referencing",
             col = grDevices::colorRampPalette(c("blue", "white", "red"))(50))
    graphics::title(sub = "n:NA, z:const, d:dropout, +:dev, x:HFnoise, c:corr, ?:RANSAC")
    plot_sens(sl, noisyPost)
    rc <- grDevices::dev.off()
    writeLines('```{r, fig.align=\'center\'}', con)
    writeLines(paste0('knitr::include_graphics("', png_zHF_before, '")'), con)
    writeLines(paste0('knitr::include_graphics("', png_zHF_after, '")'), con)
    writeLines('```', con)
    writeLines('', con)
    
    # HF noise statistics table
    
    thrHF_before <- noisyPre$cor_noise > noisyPre$HFnoisethr
    thrHF_after <- noisyPost$cor_noise > noisyPost$HFnoisethr
    mHF_before <- apply(thrHF_before, 2, mean)
    mHF_after <- apply(thrHF_after, 2, mean)
    nsPre <- apply(thrHF_before, 1, sum)
    nsPost <- apply(thrHF_after, 1, sum)
    
    mx <- max(mHF_before, mHF_after)
    
    col1 <- c(
      "High frequency noise threshold:",
      "Bad sensors because of HF noise above this threshold:",
      "Number of segments for which HF noise is calculated:",
      paste("Mean number of sensors with HF noise above",
            "threshold in 1-second segments:"),
      paste("Mean proportion of sensors with HF noise above",
            "threshold in 1-second segments:"),
      "Median HF noise:",
      "Standard Deviation HF noise:",
      paste("Number of windows with more than 25% of",
            "sensors above threshold:"),
      paste("Number of windows with more than 50% of",
            "sensors above thresholod:"),
      "Median HF noise over segments:",
      "Standard deviation HF noise over segments:"
    )
    col2 <- c(
      noisyPre$HFnoisethr,
      wrs(colnames(x)[noisyPre$BadFromHFnoise]),
      nrow(noisyPre$cor_noise),
      mnsPre <- mean(thrHF_before),
      mnsPre / ns,
      noisyPre$medHFnoise,
      noisyPre$sdHFnoise,
      sum(nsPre > round(0.25 * ns)),
      sum(nsPre > round(0.50 * ns)),
      stats::median(noisyPre$cor_noise),
      stats::sd(noisyPre$cor_noise)
    )
    col3 <- c(
      noisyPost$HFnoisethr,
      wrs(colnames(x)[noisyPost$BadFromHFnoise]),
      nrow(noisyPost$cor_noise),
      mnsPost <- mean(thrHF_after),
      mnsPost / ns,
      noisyPre$medHFnoise,
      noisyPre$sdHFnoise,
      sum(nsPost > round(0.25 * ns)),
      sum(nsPost > round(0.50 * ns)),
      stats::median(noisyPost$cor_noise),
      stats::sd(noisyPost$cor_noise)
    )
    df <- data.frame(col1, col2, col3)
    colnames(df) <- c(" ", "Before", "After")
    table <- huxtable::hux(df)
    huxtable::valign(table) <- "bottom"
    huxtable::width(table) <- 0.9
    huxtable::wrap(table) <- TRUE
    huxtable::caption(table) <- paste('High frequency noise window statistics', x_name)
    table <- huxtable::set_bottom_border(table, huxtable::brdr(1, "solid", "black"))

    writeLines(huxtable::to_latex(table), con)
    writeLines('', con)

    writeLines('\\', con)
    writeLines('\\clearpage', con)
    writeLines('', con)

    # Cumulative distribution function of noisiness scores
    zHFpre <- (noisyPre$cor_noise - mean(noisyPre$cor_noise)) / noisyPre$sdHFnoise
    zHFpost <- (noisyPost$cor_noise - mean(noisyPost$cor_noise)) / noisyPost$sdHFnoise
    png_cdfHFnoise <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_cdfHFnoise, width = 14, height = 10, units = "cm", res = 600)
    plot(stats::ecdf(zHFpre), xlim = c(-5, 5), pch = NA,
         xlab = "HF noise", ylab = "Cumulative probability",
         main = "High-frequency noise distribution")
    graphics::lines(stats::ecdf(zHFpost), col = "red", pch = NA)
    graphics::legend("left", legend = c("Before", "After"), lty = 1, col = c("black", "red"))
    rc <- grDevices::dev.off()
    writeLines('```{r, fig.align=\'center\'}', con)
    writeLines(paste0('knitr::include_graphics("', png_cdfHFnoise, '")'), con)
    writeLines('```', con)
    writeLines('', con)

    # sensors with correlation below threshold, by interval
    png_HFsens <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_HFsens, width = 14, height = 10, units = "cm", res = 600)
    mx <- ceiling(max(nsPre, nsPost))
    plot(nsPre, pch = "+", ylim = c(0, mx),
         xlab = "Intervals", ylab = paste0("# sensors (out of ", ns, ")"),
         main = paste("Number of sensors above HF noise threshold\n",
                      "(by intervals of", noisyPre$corsecs, "s)"))
    graphics::points(nsPost, pch = "o", col = "red")
    graphics::legend("topleft", legend = c("Before", "After"), pch = c("+", "o"),
                     col = c("black", "red"))
    rc <- grDevices::dev.off()
    writeLines('```{r, fig.align=\'center\'}', con)
    writeLines(paste0('knitr::include_graphics("', png_HFsens, '")'), con)
    writeLines('```', con)
    writeLines('', con)
    
    # Comparison of original and referenced data
    
    writeLines('## Comparison of data before and after referencing', con)
    
    tmp <- removetrend(x, fs)$y
    refPre <- rowMeans(tmp)
    tmp <- removetrend(rr$y, fs)$y
    refPost <- rowMeans(tmp)
    corprepost <- stats::cor(refPost, refPre)
    png_corpp <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_corpp, width = 14, height = 10, units = "cm", res = 600)
    plot(refPost, refPre, xlab = "After referencing", ylab = "Before referencing",
         main = paste0("Comparison before and after referencing\n",
                      "(corr = ", round(corprepost, 7), ")"))
    rc <- grDevices::dev.off()
    writeLines('```{r, fig.align=\'center\'}', con)
    writeLines(paste0('knitr::include_graphics("', png_corpp, '")'), con)
    writeLines('```', con)
    writeLines('', con)
    
    png_corpp_time <- tempfile(pattern = "img", fileext = ".png")
    grDevices::png(png_corpp_time, width = 14, height = 10, units = "cm", res = 600)
    plot(seq(0, (npts - 1) / fs, length.out = npts), refPre - refPost,
         type = "l", xlab = "Time (s)", ylab = "Pre - Post",
         main = "Comparison before and after referencing")
    rc <- grDevices::dev.off()
    writeLines('```{r, fig.align=\'center\'}', con)
    writeLines(paste0('knitr::include_graphics("', png_corpp_time, '")'), con)
    writeLines('```', con)
    writeLines('', con)
    
  } else {
    warning("Robust rereferencing unsuccessful. Error message:")
    warning(result)
    if (rep) {
      writeLines('\\', con)
      writeLines('Robust referencing unsuccessful. Error message:', con)
      writeLines(paste(result), con)
    }
  }

  toc <- Sys.time()
  tdiff <- toc - tic
  if (rep) {
    writeLines('', con)
    writeLines(paste('Elapsed time:', format(tdiff, digits = 4, scientific = FALSE)), con)
    writeLines('\\clearpage', con)
    writeLines('', con)
  }
  
  #############################################################################
  # render report
  
  if (rep) {
    close(con)
    suppressWarnings(rmarkdown::render(input = Rmd, output_file = report, quiet = TRUE))
  }
  
  #############################################################################
  # return value
  
  list(
    y = ctd(rr$y, fs = rr$fs, ref = "average"),
    noisyPre = noisyPre,
    noisyPost = noisyPost
  )
  
  
}
