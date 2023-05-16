# readEDF.R - Function to read an EDF/BDF file into an R data frame
# Copyright (C) 2012, 2019  Geert van Boxtel,
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
# 20120207  GvB     Original version
# 20190701  GvB     complete overhaul for eegr v0.1.0
# 20201104  GvB     link to gsignal instead of signal package
# 20230516  GvB     bugfixes in calculating record size and up/downsampling
#
#-----------------------------------------------------------------------------------------------

#' readEDF
#' 
#' Read data from a European Data Format (EDF, EDF+, BDF, BDF+) file.
#' 
#' The function readEDF reads files containing electrophysiological data stored in the European Data
#' Format (EDF) or its variants.
#' 
#' EDF was published in 1992 and stores multichannel data, allowing different sample rates for each signal.
#' Internally it includes a header and one or more data records. The header contains some general information 
#' (patient identification, start time...) and technical specs of each signal (calibration, sampling rate, filtering, ...),
#' coded as ASCII characters. The data records contain samples as little-endian 16-bit integers. 
#' 
#' EDF+ was published in 2003 and is largely compatible to EDF, but EDF+ files also allow coding discontinuous recordings
#' as well as annotations, stimuli and events in UTF-8 format.
#' 
#' Biosemi developed an EEG acquisition system with 24-bit A/D converters, which could not be stored in the 16-bit
#' EDF format. They introduced a 24-bit variant of EDF called BDF (Biosemi Data Format), and a BDF+ variant is sometimes also
#' found.
#' 
#' The function \code{readEDF} can read all of these variants and produce different output that can be used for further processing
#' by function in the \code{eegr} package.
#' 
#' @references Kemp B., Värri, A., Rosa, A.C., Nielsen, K.D. and Gade, J. (1992). A simple format for exchange of
#' digitized polygraphic recordings. \emph{Electroencephalogr Clin Neurophysiol}, \emph{82(5)}, 391-393.
#' @references Kemp, B. and Olivan, J. (2003). European data format ’plus’ (EDF+), an EDF alike standard format
#' for the exchange of physiological data. \emph{Clin Neurophysiol}, \emph{114(9)}, 1755-1761.
#' @references European Data Format. \url{https://www.edfplus.info/}
#' @references Which file format does Biosemi use? \url{https://www.biosemi.com/faq/file_format.htm}
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @param file Character string. The name of the file which the data are to be read from. \code{file} 
#' can also be a complete URL. (For the supported URL schemes,see the ‘URLs’ section of the help for url.)
#' @param records Numeric or \code{'all'}. In EDF, signals are stored in 'records', usually containing 0.5 or 1
#' second of data. This function either reads all records (default: \code{all}), or a subset of the records, thus 
#' limiting the time interval of the signals. \code{records} can be a range (e.g., \code{1:50}) or a random sequence 
#' of records (e.g., \code{c(1,5,12,3)}). In the latter case, the records are returned in the specified order (only
#' useful for discontinuous data). If \code{length(records) == 0} or \code{is.na(records)}, then no records are read
#' and only the file header is returned.
#' @param signals Numeric, character string, or \code{'all'}. Specifies the signals to read from the input file. 
#' \code{signals} can be a range (e.g., \code{1:21}) a random sequence of signals (e.g., \code{c(1,5,12,3)}),
#' or a sequence of signal labels (e.g., \code{c("F3", "Fz", "F4")}. In the latter two cases, the signals are returned 
#' in the specified order. If \code{length(signals) == 0} or \code{is.na(signals)}, then no records are read
#' and only the file header is returned.
#' @param multifreq Character string. Specifies what to do when signals with different sampling frequencies are
#' present in the file. Options are:
#' \describe{
#'   \item{\code{"group"}}{Signals with identical sampling frequencies are grouped together in a single \code{ctd} object
#'   (Default)}
#'   \item{\code{"separate"}}{Each signal will have its own \code{ctd} object}
#'   \item{\code{"downsample"}}{Signals with higher sampling frequencies are downsampled to the lowest frequency in the file}
#'   \item{\code{"upsample"}}{Signals with lower sampling frequencies are upsampled to the highest frequency in the file}
#' }
#' Up- and downsampling are performed by the \code{resample} function in the \code{\link{gsignal}} package. Note that this can 
#' take a lot of time, especially in case of upsampling to the frequency of the annotation signal (usually 1-2 kHz).
#' @param physical Logical. If true, convert the digital data values to physical values. \eqn{physval = gain * digival + offset},
#' where \eqn{gain = (physmax - physmin) / (digimax - digimin)}, and \eqn{offset = physmax - gain * digimax}. See the description
#' of the header below for the definition of the variables used in these equations.
#' 
#' @return \code{list (head, ctd, ..., anno)}. The function returns a \code{list} containing three possible elements: the
#' file header, the data, and annotations, if present. The header is always returned; the other two elements depend on the
#' input arguments.
#' \describe{
#'   \item{\strong{header}}{A list named \code{'header'} consisting of the following elements:
#'     \describe{
#'       \item{General header:\cr}{
#'         \describe{
#'           \item{\code{id0}}{First identification byte. For BDF: -1 (0xFF), for EDF: ASCII "0" (0x30, old) or " " (0x32, new)}
#'           \item{\code{id1 }}{Identification string, "" for EDF, "BIOSEMI" for BDF (7 bytes ASCII)}
#'           \item{\code{lsi}}{Local subject identification (80 bytes ASCII)}
#'           \item{\code{lri}}{Local recording identification (80 bytes ASCII)}
#'           \item{\code{date}}{Starting date of recording (8 bytes ASCII - "dd.mm.yy")}
#'           \item{\code{time}}{Starting time of recording (8 bytes ASCII - "hh.mm.ss")}
#'           \item{\code{nbytes}}{Number of bytes in header record (converted from 8 bytes ASCII to numeric)}
#'           \item{\code{version}}{Reserved. Used for version of format. Empty for EDF, "EDF+C" for EDF continuous data,
#'             "EDF+D" for EDF discontinuous data, "BIOSEMI" or "24BIT" for BDF (44 bytes ASCII)}
#'           \item{\code{nrec}}{Number of data records (converted from 8 bytes ASCII to numeric), -1 if unknown}
#'           \item{\code{duration}}{Duration of one data record in seconds (converted from 8 bytes ASCII to numeric)}
#'           \item{\code{ns}}{Number of signals in each data record (converted from 4 bytes ASCII to numeric)}
#'         }
#'       }
#'       \item{Signal-specific header (one for each signal in the file):\cr}{
#'         \describe{
#'           \item{\code{label}}{Signal labels, e.g. "Fpz", "Fz" (ns * 16 bytes ASCII)}
#'           \item{\code{transducer}}{Transducer type, e.g. "AgAgCl electrode" (ns * 80 bytes ASCII)}
#'           \item{\code{dimension}}{Physical dimension of signal, e.g., "uV", "Ohm" (ns * 8 bytes ASCII)}
#'           \item{\code{physmin}}{Physical minimum in units of dimension (ns * 8 bytes ASCII converted to numeric)}
#'           \item{\code{physmax}}{Physical maximum in units of dimension (ns * 8 bytes ASCII converted to numeric)}
#'           \item{\code{digimin}}{Digital minimum ((ns * 8 bytes ASCII converted to numeric)}
#'           \item{\code{digimax}}{Digital maximum (ns * 8 bytes ASCII converted to numeric)}
#'           \item{\code{gain}}{Gain factor used for mapping digital to physical values (ns * numeric)}
#'           \item{\code{offset}}{Offset value used for mapping digital to physical values (ns * numeric)}
#'           \item{\code{prefilt}}{Prefiltering, e.g., "TC:3s;LP:70Hz" (ns * 80 bytes ASCII)}
#'           \item{\code{nsamp}}{Number of samples in one data record (ns * 8 bytes ASCII converted to numeric)}
#'           \item{\code{resvd}}{Reserved (ns * 32 bytes ASCII)}
#'         }
#'       }
#'     }
#'   }
#'   \item{\strong{data}}{One or more Continuous Time Domain (\code{\link{ctd}}) objects. If the input file contains signals
#'     which all have the same sampling frequency, or if either the \code{"upsample"} or \code{"downsample"} options were 
#'     specified in the \code{multifreq} parameter, then there is a single \code{ctd} object, and it is named '\code{ctd}'.\cr\cr
#'     In case of multiple sampling frequencies, then if the \code{multifreq = "group"}, as many \code{ctd} objects are
#'     returned as there are sampling frequencies, and they are named \code{ctdF}, where F equals the sampling frequency
#'     of the data in that object. For instance, for a data file with sampling frequencies of 256 and 1024 Hz, two \code{ctd}
#'     objects will be returned, named \code{ctd256} and \code{ctd1024}.\cr\cr
#'     In case \code{multifreq = "separate"} is specified, there will be as many \code{ctd} objects as there are signals, 
#'     and they will be named \code{'ctd1...ctdN'}, where \code{N} equals the number of signals in the file.
#'   }
#'   \item{\strong{annotations}}{A (list of) data frame(s) named \code{'anno'} containing the annotations parsed from
#'   the EDF+/BDF+ annotation channels, or from the Biosemi 'status' channel. Each annotation consists of the following fields:
#'     \describe{
#'       \item{\code{sample}}{sample number of annotation}
#'       \item{\code{time}}{corresponding time, relative to beginning of file}
#'       \item{\code{value}}{value of annotation, as per EDF+ specification, see \url{https://www.edfplus.info/specs/edfplus.html}.
#'         For Biosemi this ia a 16 bit number, as only bits 0-15 of 24 are used as triggers, see
#'         \url{http://biosemi.com/faq/trigger_signals.htm}}.
#'     }
#'     For EDF+ and BDF+ files, the Time-stamped Annotation Lists (TALs) will also be returned. In this case, sample numbers and
#'     times are relative to the beginning of the file, irrespective of the records read.
#'   }
#' }
#' 
#' @examples 
#' \dontrun{
#'   edf <- readEDF('test.edf', signals = 1:10, records = "all")
#'   summary(edf$ctd)
#'   header <- readEDF(file.bdf, records = NA)
#'   header
#' }
#' @importFrom mmap mmap munmap int16 int24 char
#' @importFrom gsignal resample
#' @export


# file <- "/media/geert/DATA/Data/Stefan/Robin/edf/p05.edf"
# signals=28:30
# signals<-c("F3","Fz","F4")
# #signals=c(1:21, 23,24)
# #signals=22
# records=NULL

readEDF <- function (file, records = "all", signals = "all", 
                     multifreq = c("group", "separate", "downsample", "upsample"),
                     physical = TRUE) {

  # helper function to trim leading and trailing spaces from header strings
  trim.spaces <- function(a) gsub("^ *","",gsub(" *$","",a))
  
  # check parameters
  if (missing(file) || is.null(file) || !is.character(file)) stop(paste("Invalid 'file' parameter:", file))
  if (!file.exists(file)) stop(paste(file, "does not exist"))

  # Open input file and read header variables
  fp <- file(file, "rb")
  if (!isOpen(fp)) stop(paste("Unable to open input file:", file))

  result <- tryCatch({
    id0 <- readBin(fp, integer(), n = 1L, size = 1L, endian = "little")
    id1 <- readChar(fp, 7)
    lsi <- trim.spaces(readChar(fp, 80))
    lri <- trim.spaces(readChar(fp, 80))
    date <- readChar(fp, 8)
    time <- readChar(fp, 8)
    nbytes <- as.numeric(readChar(fp, 8))
    version <- trim.spaces(readChar(fp, 44))
    nrec <- as.numeric(readChar(fp, 8))
    if (nrec <= 0) stop(paste("\nInvalid number of records, nrec =", nrec))
    duration <- as.numeric(readChar(fp, 8))
    ns <- as.numeric(readChar(fp, 4))
    if (ns <= 0) stop(paste("\nInvalid number of signals, ns =", ns))
    label <- array('', ns)
    for (is in seq_len(ns)) {
      label[is] <- trim.spaces(readChar(fp, 16))
    }
    transducer <- array('', ns)
    for (is in seq_len(ns)) {
      transducer[is] <- trim.spaces(readChar(fp, 80))
    }
    dimension <- array('', ns)
    for (is in seq_len(ns)) {
      dimension[is] <- trim.spaces(readChar(fp, 8))
    }
    physmin <- array(0, ns)
    for (is in seq_len(ns)) {
      physmin[is] <- as.numeric(readChar(fp, 8))
    }
    physmax <- array(0, ns)
    for (is in seq_len(ns)) {
      physmax[is] <- as.numeric(readChar(fp, 8))
    }
    digimin <- array(0, ns)
    for (is in seq_len(ns)) {
      digimin[is] <- as.numeric(readChar(fp, 8))
    }
    digimax <- array(0, ns)
    for (is in seq_len(ns)) {
      digimax[is] <- as.numeric(readChar(fp, 8))
    }
    gain <- (physmax - physmin) / (digimax - digimin)
    offset <- physmax - gain * digimax
    prefilter <- array('', ns)
    for (is in seq_len(ns)) {
      prefilter[is] <- trim.spaces(readChar(fp, 80))
    }
    nsamp <- array(0, ns)
    for (is in seq_len(ns)) {
      nsamp[is] <- as.numeric(readChar(fp, 8))
    }
    resvd <- array(0, ns)
    for (is in seq_len(ns)) {
      resvd[is] <- trim.spaces(readChar(fp, 32))
    }
  },  #end of tryCatch expression
    warning = function(w) warning(paste("Warning reading header:", w)),
    error = function(e) stop(paste("Error reading header:", e))
  ) #end of tryCatch block

  # check header size
  if (nbytes != ns*256 + 256) warning("inconsistent header size, attempting to continue...")
  
  head <- list(id0 = id0, id1 = id1, lsi = lsi, lri = lri, date = date, time = time,
               nbytes = nbytes, version = version, nrec = nrec, duration = duration,
               ns = ns, label = label, transducer = transducer, dimension = dimension,
               physmin = physmin, physmax = physmax, digimin = digimin,
               digimax = digimax, gain = gain, offset = offset, prefilter = prefilter,
               nsamp = nsamp, resvd = resvd)

  
  # Determine how many records and signals to read. Return header if either
  # no records or no signals are requested
  if (missing(records) || is.null(records)  || trim.spaces(records) == 'all') records <- seq_len(nrec)
  if (missing(signals) || is.null(signals) || trim.spaces(signals) == 'all') signals <- seq_len(ns)
  if (length(records) == 0 || (length(records) == 1 && (is.na(records) || records == 0)) ||
      length(signals) == 0 || (length(signals) == 1 && (is.na(signals) || signals == 0))) {
    close(fp)
    return(head)
  }
  
  # Translate named signals as numbers
  if (typeof(signals) == "character") signals <- match(signals, label)

  # Avoid array out of bounds conditions
  records <- unique(abs(records[!is.na(records) & records > 0 & records <= nrec]))
  signals <- unique(abs(signals[!is.na(signals) & signals > 0 & signals <= ns]))
  
  # Read all samples as raw() and write them back to a temporary file which will
  # then be mapped into memory by the mmap function. This is a bit clumsy,
  # but using mmap is the fastest way to read 24-bit integers. Unfortunately,
  # I have found no way to map the entire original file (including the header) into memory.
  # The mmap function does have an offset variable, but it can only  be equal to the system pagesize.
  # So, the solution is to write back the samples to disk withoug the header, and to map this temporary
  # file to memory starting with the first byte. Tests have shown that this is faster
  # than fooling around with bytes in a big loop.

  # Determine whether to read 16-bit (EDF) or 24-bit (BDF) integers
  size <- ifelse(id0==-1 && substr(id1,1,7)=='BIOSEMI' && 
                   (substr(version,1,5)=='24BIT' || substr(version,1,3) == 'BDF'), 3, 2)

  # open temporary file for writing single bytes
  tmpf <- tempfile()
  tmp <-  file(tmpf, "wb")
  if (!isOpen(tmp)) {
    message (paste("Unable to open temporary file for writing:", tmpf))
    close (fp)
    return(head)
  }

  # read everything at once
  bytes2read <- size * sum(nsamp) * nrec
  sample <- readBin (fp, raw(), n=bytes2read, size=1L, endian="little")
  writeBin (sample, tmp, raw())
  close (fp)
  close (tmp)

  # Map the temporary file into memory, either as 24-bit int or as 16-bit int
  if (size == 3) {
    smp <- mmap::mmap(tmpf, mode=mmap::int24())
  } else {
    smp <- mmap::mmap(tmpf, mode=mmap::int16())
  }

  recsize <- sum(nsamp)                         # total record size including the signals we do not want to read
  if (length(smp) != recsize * nrec) {          # check against length of smp
    message (paste("Header/data mismatch",
                   "\nHeader indicates to read", recsize * nrec, "bytes,",
                   "\nbut", length(smp), "bytes were read"))
    return(head)
  }
  
  # Get the signals we want from the records we want and move them into temporay arrays
  # these arrays can be of different length
  dat <- ptr <- list()
  for (s in signals) {
    dat[[s]] <- array(0, length(records) * nsamp[s])
    ptr[[s]] <- 1
  }
  for (r in records) {
    for (s in signals) {
      smp.f <- 1 + (r - 1) * recsize + ifelse(s == 1, 0, sum(nsamp[1:(s-1)]))
      smp.t <- smp.f + nsamp[s] - 1
      dat.f <- ptr[[s]]
      dat.t <- ptr[[s]] + nsamp[s] - 1
      dat[[s]][dat.f:dat.t] <- smp[smp.f:smp.t]
      ptr[[s]] <- ptr[[s]] + nsamp[s]
    }
  }

  # convert digital to physical values, if requested
  if (physical) {
    for (s in signals) {
      if (label[s] != 'EDF Annotations' && label[s] != 'BDF Annotations' && label[s] != 'Status') {
        dat[[s]] <- dat[[s]] * gain[s] + offset[s]
      }
    }
  }

  # Determine how the data part of the output should look like.
  multifreq <- match.arg(multifreq)
  
  # 1. If there is a single sampling frequency, then make a single ctd object
  freq <- nsamp / duration
  frex <- unique(freq[signals])
  nfrex <- length(frex)
  if (nfrex == 1) {
    ctd <- ctd(matrix(unlist(dat[signals]), nrow = frex * duration * length(records), ncol = length(signals)), frex)
    colnames(ctd)[1:ns(ctd)] <- label[signals]
  }
  
  # 2. If there are multiple frequencies and multifreq == 'group', then group all signals with
  # the same sampling frequencies into one ctd object. Name the ctd objects ctdF where F is the
  # sampling frequency.
  if (nfrex > 1 && multifreq == 'group') {
    ctd <- mtx <- list()
    m <- match(freq, frex)[signals]     #match the selected signals to the frex
    tab <- table(m)                     #number of signal per freq
    ptr <- array(1, nfrex)              # keep track of column
    for (f in seq_along(frex)) {
      mtx[[f]] <- matrix(0, nrow = frex[f] * duration * length(records), ncol = tab[f])
      cn <- NULL
      for (s in signals) {
        if (freq[s] == frex[f]) {
          mtx[[f]][, ptr[f]] <- dat[[s]]
          ptr[f] <- ptr[f] + 1
          cn <- c(cn, label[s])
        }
      }
      colnames(mtx[[f]]) <- cn
      ctd[[f]] <- ctd(mtx[[f]], frex[f])
    }
    names(ctd) <- paste0('ctd', frex)
  }
  
  # 3.   If there are multiple frequencies and multifreq == 'separate', then there will be a ctd object for
  # each signal, and they will be named ctd1...ctdN, where N is the number of signals
  if (nfrex > 1 && multifreq == 'separate') {
    ctd <- list()
    ptr <- 1
    for (s in signals) {
      ctd[[ptr]] <- ctd(dat[[s]], freq[s])
      colnames(ctd[[ptr]])[1] <- label[s]
      ptr <- ptr + 1
    }
    names(ctd) <- paste0('ctd', 1:length(ctd))
  }
  
  # 4. If there are multiple frequencies and "upsample" is specified, then determine the highest
  # sampling rate and upsample every signal to that frequency
  # 5. If there are multiple frequencies and "downsample" is specified, then determine the lowest
  # sampling rate and downsample every signal to that frequency
  if (nfrex > 1 && (multifreq == "upsample" || multifreq == "downsample")) {
    newf <- ifelse(multifreq == "upsample", max(frex), min(frex))
    mtx <- matrix(0, nrow = newf * duration * length(records), ncol = length(signals))
    col <- 1
    for (s in signals) {
      if (freq[s] != newf) mtx[, col] <- gsignal::resample(as.vector(dat[[s]]), newf, freq[s])
      else mtx[, col] <- as.vector(dat[[s]])
      col <- col + 1
    }
    ctd <- ctd(mtx, newf)
    colnames(ctd)[1:ns(ctd)] <- label[signals]
  }
  
  # Annotations are ALWAYS parsed if present in the file. For EDF files, there is no annotation channel. 
  # For Biosemi generated BDF files, there is always a 'Status' channel which is the last signal in the file. 
  #Try to find this channel in the selected signals.
  anno <- NULL
  if (id0 == -1 && substr(id1, 1, 7) == 'BIOSEMI' && substr(version,1,5)=='24BIT') {
    annosig <- which(label[signals] == 'Status')
    nanno <- length(annosig)
    if (nanno > 1) annosig <- annosig[1]
    if (nanno != 0) {
      df <- diff(dat[[annosig]])
      trig <- which(df>0 & df <= 2^16)          # Biosemi uses 16-bit triggers
      if (length(trig) > 0) {
        as <- trig + 1                          # high level is moment of trigger
        at <- (trig / (nsamp[annosig] / duration))
        av <- df[trig]
        anno <- data.frame(sample=as, time=at, value=av);
      }
    }
  } else {
    # For B/EDF+C and B/EDF+D files, try to find the annotation channels among the selected signals.
    # These signals will contain Time-stamped Annotation Lists (TALs), containing single-byte raw characters.
    # Re-open the temporary file as char(). 
    getTALs <- function (chan, sf) {
      TALs <- NULL
      start0 <- which(chan == 0x2b | chan == 0x2d)     #TAL starts with '+' or '-'
      # bugfix 20191210
      if (length(start0) == 0) return(NULL)
      start1 <- start0 - 1                             #preceded by 0x00
      if (start0[1] == 1) start <- c(1, start0[2:length(start0)][which(chan[start1] == 0)])
      else start <- start0[which(chan[start1] == 0)]
      end20 <- which(chan == 20)                       #TAL ends with 20
      end0 <- end20 + 1                                #followed by 0
      end <- end20[which(chan[end0] == 0)]             #last sample handled implicitly
      n <- length(start)
      if (n == length(end)) {
        for (tal in seq_along(start)) TALs <- rbind(TALs, rawToChar(chan[start[tal]:end[tal]]))
      }
      parseTALs(TALs, sf)
    }
    parseTALs <- function(tals, sf) {
      if(nrow(tals) <= 0) return(NULL)
      s <- t <- d <- v <- array(0, nrow(tals))
      for (r in seq_len(nrow(tals))) {
        split1 <- unlist(strsplit(tals[r], split = "\024"))    #octal 24 = hex 14 = dec 20
        split2 <- unlist(strsplit(split1, split = "\025"))     #octal 25 = hex 15 = dec 21
        suppressWarnings(t[r] <- as.numeric(split2[1]))        #allow NAs by coercion
        s[r] <- round(t[r] * sf) + 1
        suppressWarnings(d[r] <- as.numeric(split2[2]))        #allow NAs by coercion
        v[r] <- trim.spaces(split1[2])
      }
      data.frame(sample = s, time = t, value = v, duration = d, tal = tals,
                 stringsAsFactors = FALSE)
    }
    m <- match(label, c('EDF Annotations', 'BDF Annotations'))
    annosig <- which(m > 0)
    sig <- anno <- list()
    if (!is.null(annosig)) {
      smp <- mmap::mmap(tmpf, mode=mmap::char())
      for (s in seq_along(annosig)) {
        sig[[s]] <- raw(length = size * length(records) * nsamp[annosig[s]])
        ptr <- 1
        for (r in seq_len(nrec)) {
          smp.f <- 1 + (r - 1) * recsize * size + sum(nsamp[1:(annosig[s]-1)] * size)
          smp.t <- smp.f + nsamp[annosig[s]] * size - 1
          dat.f <- ptr
          dat.t <- ptr + nsamp[annosig[s]] * size - 1
          sig[[s]][dat.f:dat.t] <- smp[smp.f:smp.t]
          ptr <- ptr + nsamp[annosig[s]] * size
        }
        anno[[s]] <- getTALs(sig[[s]], nsamp[annosig[s]] / duration)
      }
    }
    if(length(anno) == 1) anno <- anno[[1]]    # no point in returning a list if there is only one anno channel
  }
  
  # Clean up
  mmap::munmap(smp)
  unlink(tmpf)

  # Ready
  return(list(header = head, ctd = ctd, anno = anno))

}


# s <- sprintf("%04X", ctd$data[1,1])
# h <- sapply(seq(1, nchar(s), by=2), function(x) substr(s, x, x+1))
# rawToChar(as.raw(strtoi(h[1])))
# 
# t <- tempfile()
# fp <- file(t, "wb")
# writeBin(ctd$data[1:5,1], fp, numeric())
# close(fp)
# fp <- file(t, "rb")
# asc <- readBin(fp, character(), n=10)
# close(fp)