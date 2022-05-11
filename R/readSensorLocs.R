# readSensorLocs.R
# Read sensor locations from an external file or data frame.
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
# Version history:
# 20160622    GvB           initial version
# 20190912    GvB           version for eegr v0.1.0
# 20220511    GvB           use inherits() instead of class(type/sens)
#------------------------------------------------------------------------------

#' Read Sensor Locations
#'
#' Read sensor locations from an external file or select locations from a data
#' frame.
#'
#' The function \code{readSensorLocs} reads sensor coordinates from a variety of
#' input sources. Some well-known file types are supported, such as BESA
#' spherical coordinates, and EEGLAB and Polhemus Cartesian coordinates. These
#' coordinate systems are converted by this function to the coordinate system
#' used by \code{eegr} (see \code{\link{sensorlocs}}).
#' 
#' Alternatively, a \code{data.frame} may be specified as the source. This may
#' be useful if standard sensor locations are used, of which a particular
#' experiment uses a selection (see the \code{select} parameter).
#' 
#' In addition, sensor locations may be read from an external file with a custom
#' format. In this case (\code{type = 'custom'}), the function expects a format
#' string that specifies the coordinates (see the \code{format} argument).
#' 
#' @param source Character string. The source file or data frame containing the
#'   sensor location coordinates. \code{source} can also be a complete URL. (For
#'   the supported URL schemes, see the ‘URLs’ section of the help for
#'   \code{url}). Alternatively, source may be an object of class
#'   \code{\link{sensorlocs}} or a compatible data frame from which sensor
#'   locations may be selected.\cr Default: \code{EEGlocations}, the
#'   \code{sensorlocs} object that comes with \code{eegr}.
#' @param type Character string, only used if \code{source} denotes a filename.
#'   \code{type} describes the type of file that the sensor locations are read
#'   from. Can be any of:
#' \describe{
#'   \item{besa}{BESA spherical coordinates (3,4, or columns)}
#'   \item{elp}{alias for \code{besa}}
#'   \item{eeglab}{Cartesian coordinates as used in the EEGLAB package}
#'   \item{xyz}{Cartesian Polhemus coordinates}
#'   \item{custom}{custom input type, see the \code{format} parameter}
#' }
#' @param format Character string. Format of data to read when using a
#'   \code{custom} source type. The format string should consist of
#'   keyword-expression combinations of the form \code{keyword = (expr)}. The
#'   keyword-expression combinations recognized are:
#' \describe{
#'   \item{\code{vars = (var1, var2, ...)}}{The variables to read from the file.
#'   Variables can be either label, theta, phi, x, y, z, or dummy, separated by
#'   commas. The variables should at least contain either theta and phi, or x,
#'   y, and z. A label is not required, but recognized.}
#'   \item{\code{adjust = (expr, expr, ...)}}{used for adjusting \code{vars}.
#'   \code{expr} is either an expression applied to the value read (e.g.,
#'   '*-1'), or NULL. There should be just as many expressions as there are
#'   \code{vars}.}
#'   \item{\code{skiplines = (x)}}{number of header lines to skip}
#' }
#' For example: \code{format='vars = (theta, phi, label),}\cr
#' \code{adjust = (*90/72, NULL, NULL), skip = (2)'}
#' @param select Numeric or character vector, only used if \code{source} denotes
#'   a data frame. Specify '\code{all}' (default) to select all sensors from the
#'   \code{source} data frame. Alternatively, a character vector indicating the
#'   labels of the sensors to be selected (a variable named \code{label} should
#'   then be present in the data frame), or a numeric vector specifying the rows
#'   from the data frame to be selected.
#' @param plot logical, default: FALSE. If \code{TRUE}, make a 2D plot of all
#'   sensor locations (top view, nose pointing upward). Alternatively
#'   \code{plot} may be an array of character strings denoting the sensor labels
#'   to plot.
#' 
#' @examples 
#' \dontrun{
#'   readSensorLocs ('/path/to/filename', type = 'besa', plot = TRUE)
#'   readSensorLocs ('/path/to/filename', type = 'xyz', plot = c('Fz', 'Cz', 'Pz'))
#'   readSensorLocs (source = '/somestrangeformat.txt', type = 'custom',
#'                format='vars = (theta, phi, dummy), adjust=(*90/72, NULL, NULL), skip = (2)')
#' }
#' 
#' @return A data frame of class \code{\link{sensorlocs}} containing x, y, z,
#'   x2d, y2d, theta, phi in the eegr internal format. Other variables in the
#'   input data frame or matrix are copied unchanged to the outpuut data frame.
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname readSensorLocs
#' @export

readSensorLocs <- function (source = 'EEGlocations',
                         type = c('besa', 'elp', 'eeglab', 'xyz', 'custom'),
                         format = '',
                         select = 'all',
                         plot = FALSE) {
  
  # define some functions and variables
  deg2rad <- (2 * pi) / 360                                       # degrees to radians
  remsp <- function(str) gsub('[[:space:]]', '', str)             # remove spaces from string
  index <- function (subs, str, from = 1) {
    as.numeric(regexpr(subs, substring(str, from), fixed = TRUE)) # find first occurence of 'substr' in 'str' starting from 'from'
  }
  
  # source should be a character string specifying an external file or to a data frame
  if (is.null(source) || !is.character(source)) {
    stop ("'source' parameter should be a character string")
  }
  # If no source was specified, then assume 'EEGlocations' dataset
  if (missing (source) || source == '' || length(source) == 0 || length(remsp(source)) == 0) {
    source <- 'EEGlocations'
  }
  # test existence of data frame or external file
  df <- exists(source)
  ef <- file.exists(source)
  if (!df && !ef) stop(paste(source, "does not exist"))

  # if source is a data.frame, then we should have a select argument (character or numeric)
  if (df) {
    if (missing(select) || is.null(select) ||
        length(select) == 0 || length(remsp(select)) == 0 ||
        (!is.character(select) && !is.numeric(select))) {
      stop(paste("Invalid 'select' variable: ", select))
    }
  } else if (ef) {
  # if source is external file, then we should either have a matching type, and if type=='custom'
  # then we should also have a 'format' (character). If no type is specified, then we try to guess
  # the file type from the file name extension
    if (!missing(type) && !is.null(type)) {
      type <- try(match.arg(type), silent = TRUE)
      if (inherits(type, "try-error")) {
        stop (paste("invalid 'type' value:", type))
      }
      if (type == 'custom' && (missing(format) || is.null(format))) {
        stop ("custom type specified, but no format parameter present")
      } else {
        if (!is.character(format)) {
          stop ("'format' parameter should be a character string");
        }
      }
    } else {
      # no type specified; try to guess it from the file extension
      ext <- unlist(strsplit(basename(source), '.', fixed = TRUE));
      ext <- tolower(ext[length(ext)])
      if (ext == 'elp' || ext == 'eps') {
        type <= 'besa'
      } 
      if (ext == 'loc' || ext == 'locs') {
        type <- 'eeglab'
      }
      if (ext == 'xyz') {
        type <- 'xyz'
      }
      # if still not found, then issue an error
      if (is.null(type)) {
        stop ('type could not be established - use a custom format')
      }
    }
  }
  
  # Read BESA format.
  # BESA uses an angular coordinate system that is virtually the same as the
  # internal coordinate system. There, however, one important difference. In
  # BESA the X-axis points through T4, whereas in our internal system, the
  # X-axis points to the right periauricular point. Similarly, the Y-axis in
  # BESA runs through Fpz, and in our system through the nasion. So theta has
  # to be adjusted by a factor 72/90 degrees. For instance, Fpz has a theta of
  # 90 degrees in BESA, which will become 72 degrees in the present system.
  # All angles are rounded to units of 1 degree.
  
  if (ef && (type == 'besa' || type == 'elp')) {
    
    # try to open file
    try(fp <- file(source, "r"), silent = TRUE)
    if (!isOpen(fp)) stop ('unable to open source file', source)
    lines <- readLines(fp)
    close (fp)
    
    # make a data frame, skip blank lines and lines with less than 3 columns
    # also skip lines without three numeric variables
    sens <- NULL
    for (l in seq_along(lines)) {
      if (nchar(lines[l]) > 0) {
        t <- utils::read.table(text = lines[l], stringsAsFactors = FALSE)
        if (ncol(t) > 2 && length(grep("[[:digit:]]", t)) >= 3) {
          sens <- rbind(sens, t)
        }
      }
    }
    # no lines added to data frame: unknown error
    if (is.null(sens) || nrow(sens) <= 0) stop ("unknown error while reading BESA electrode file.")
    
    # BESA format can be 3, 4 or 5 columns
    nc <- ncol(sens)
    if (nc < 3 || nc > 5) {
      stop ("invalid number of columns in BESA electrode file.")
    }
    if (ncol(sens) == 5) {
      colnames(sens) <- c('type', 'label', 'theta', 'phi', 'radius')
    }
    if (ncol(sens) == 4) {
      colnames(sens) <- c('label', 'theta', 'phi', 'radius')
    }
    if (ncol(sens) == 3) {
      colnames(sens) <- c('theta', 'phi', 'radius')
    }
    sens$theta <- sens$theta * 72/90
    sens <- sensorlocs(sens)
  }
  
  # Read EEGlab format.
  # EEGlab used a Cartesian coordinate system consisting of x, y, z,
  # and a label. In this system, x points toward the nose (our y-axis),
  # y toward the left ear (our -x), and z toward the vertex (same).
  # The X and Y axes run through Fpz and T5, not though Nz and the preauricular point,
  # so the coordinates need to be adapted just as in the besa format.
  # (factor 72/90 degrees).
  
  if (ef && type == 'eeglab') {
    
    # try to read file
    try(fp <- file(source, "r"), silent = TRUE)
    if (!isOpen(fp)) stop ('unable to open source file', source)
    lines <- readLines(fp)
    close (fp)
    
    # make a data frame, skip blank lines and lines with less than 5 columns
    # also skip lines without four numeric variables
    sens <- NULL
    for (l in seq_along(lines)) {
      if (nchar(lines[l]) > 0) {
        t <- utils::read.table(text=lines[l], stringsAsFactors = FALSE)
        if (ncol(t) == 5 && length(grep("[[:digit:]]", t)) >= 4) {
          sens <- rbind(sens, t)
        }
      }
    }
    # no lines added to data frame: unknown error
    if (is.null(sens) || nrow(sens) <= 0) stop ("unknown error while reading EEGLAB electrode file.")
    nc <- ncol(sens)
    if (nc != 5) {
      stop ("invalid number of columns in EEGLAB electrode file.")
    }
    
    x <- -sens[, 3]
    y <- sens[, 2]
    z <- sens[, 4]
    label <- sens[, 5]
    # compute theta and phi (needed for 2D projection)
    r <- sqrt(x^2 + y^2 + z^2)
    theta <- acos(z / r) * 72/90 / deg2rad
    theta[which(x < 0)] <- -theta[which(x < 0)]
    phi <- atan(y / x) / deg2rad
    phi[which(x == 0 && y == 0)] <- 0
    phi[which(x == 0 && y > 0)] <- 90
    phi[which(x == 0 && y < 0)] <- -90
    # pass theta and phi to sensorlocs and compute corrected x y z coordinates
    sens <- sensorlocs(cbind(label,theta, phi))
  }
  
  # 
  # Read xyz format. Not sure where this orginates from, probably EGI xyz.
  # It is a Cartesian format in which x points toward inion (our -y),
  # y toward the right ear (our x), and z toward the vertex (same).
  # Unsure at this point where the origin of the coordinate system is;
  # probably in the middle of the head as in our system
  
  if (ef && type == 'xyz') {
    
    # try to read file
    try(fp <- file(source, "r"), silent = TRUE)
    if (!isOpen(fp)) stop ('unable to open source file', source)
    lines <- readLines(fp)
    close (fp)
    
    # make a data frame, skip blank lines and lines with less than 5 columns
    # also skip lines without at least 4 numeric variables
    sens <- NULL
    for (l in seq_along(lines)) {
      if (nchar(lines[l]) > 0) {
        t <- utils::read.table(text = lines[l], stringsAsFactors = FALSE)
        if (ncol(t) == 5 && length(grep("[[:digit:]]",t)) >= 4) {
          sens <- rbind(sens, t)
        }
      }
    }
    # no lines added to data frame: unknown error
    if (is.null(sens) || nrow(sens) <= 0) stop ("unknown error while reading XYZ electrode file.")
    nc <- ncol(sens)
    if (nc != 5) {
      stop ("invalid number of columns in XYZ electrode file.")
    }
    x <- sens[, 3]
    y <- -sens[ ,2]
    z <- sens[ ,4]
    label <- sens[, 5]
    # compute theta and phi (needed for 2D projection)
    r <- sqrt(x^2 + y^2 + z^2)
    theta <- acos(z / r) * 72/90 / deg2rad
    theta[which(x < 0)] <- -theta[which(x < 0)]
    phi <- atan(y / x) / deg2rad
    phi[which(x == 0 && y == 0)] <- 0
    phi[which(x == 0 && y > 0)] <- 90
    phi[which(x == 0 && y < 0)] <- -90
    # pass theta and phi to sensorlocs and compute corrected x y z coordinates
    sens <- sensorlocs(cbind(label,theta, phi))
  }
  
  # Custom type. Parse the format string
  if (ef && type == 'custom') {
    
    #remove spaces
    format <- tolower(remsp(format))
    
    #extract everything between parentheses
    m <- gregexpr("\\([^)]*\\)", format)
    left <- unlist(regmatches(format, m, invert = TRUE))
    right <- unlist(regmatches(format, m, invert = FALSE))
    # remove leading and trailing parentheses
    for (i in seq_along(right)) {
      right[i] <- sub("^\\(", '', right[i])
      right[i] <- sub("\\)$", '', right[i])
    }
    
    # find 'vars' (required)
    vars <- NA
    for (i in seq_along(left)) {
      if (regexpr('vars=', left[i]) > 0) vars <- i
    }
    if (is.na(vars)) stop ('custom format requested, but no variables detected.')
    
    # either theta and phi should be present, or x, y, z
    theta <- phi <- x <- y <- z <- FALSE
    split.vars <- unlist(strsplit (right[vars], "[;|,][[:blank:]]*"))
    for (i in seq_along(split.vars)) assign((split.vars[i]), TRUE)
    if (!((theta&&phi) || (x&&y&&z))) stop ('custom format requested, but invalid variables string detected.')

    # find adjust (optional)
    adjust <- NA
    for (i in seq_along(left)) {
      if (regexpr('adjust=', left[i]) > 0) adjust <- i
    }
    if (!is.na(adjust)) {
      split.adj <- unlist(strsplit(right[adjust], "[;|,][[:blank:]]*"))
      ls <- length (split.adj)
      adj <- rep('*1', ls)
      if (ls > 1) {
        for (i in 1:ls) {
          if (tolower(split.adj[i]) != 'null') adj[i] <- split.adj[i];
        }
      }
    }
    
    # find skiplines (optional)
    skip <- NA
    for (i in seq_along(left)) {
      if (regexpr('skip=', left[i]) > 0) skip <- i
    }
    if (!is.na(skip)) skipl <- as.numeric(right[skip])
    else skipl <- 0
    
    # Try to read the file and make the elecs data frame
    sens <- try(utils::read.table(source, stringsAsFactors=FALSE, skip=skipl), silent = TRUE)
    if (inherits(sens, "try-error")) {
      stop (paste("unable to read sensor locations file", source));
    }
    
    # Does the number of columns match the number of variables?
    # If not, just ignore the variables not listed
    for (i in seq_along(sens)) {
      names(sens)[i] <- split.vars[i]
    }
    # Do we have theta and phi?
    if (theta && phi) {
      # determine location of theta
      for (i in seq_along(sens)) {
        if (colnames(sens[i]) == 'theta') break
      }
      # apply adjust
      eval(parse(text=paste('sens$theta <- sens$theta',adj[i])))
      # same for phi
      for (i in seq_along(sens)) {
        if (colnames(sens[i]) == 'phi') break
      }
      eval(parse(text=paste('sens$phi <- sens$phi',adj[i])))
    }
    # Or x, y, and z?
    if (x && y && z) {
      for (i in 1:ncol(sens)) {
        if (colnames(sens[i]) == 'x') break
      }
      eval(parse(text=paste('sens$x <- sens$x',adj[i])));
      for (i in 1:ncol(sens)) {
        if (colnames(sens[i]) == 'y') break
      }
      eval(parse(text=paste('sens$y <- sens$y',adj[i])))
      for (i in 1:ncol(sens)) {
        if (colnames(sens[i]) == 'z') break
      }
      eval(parse(text=paste('sens$z <- sens$z',adj[i])))
    }
    sens <- sensorlocs (sens)
  }
  
  # if source is a data frame, then select a subset (or all sensors if requested)
  if (df) {
    source <- eval(parse(text = source))
    if (length(select) == 1 && select == 'all') sens <- source
    else {
      sens <- NULL
      for (i in seq_along(select)) {
        if (is.numeric(select[i])) sens <- rbind(sens, source[as.numeric(select[i]),])    #number in select
        else sens <- rbind(sens, source[which(source$label==select[i]),])                 #string in select
      }
      sens <- sensorlocs(sens)
    }
  }

  # Make plot if requested
  if (is.logical(plot)) {
    if (plot) plot (sens)
  } else if (is.character(plot)) {
    ss <- subset(sens, label %in% plot)
    plot(ss)
  } else stop (paste("invalid 'plot' value:", plot))
  
  # ready
  sens
}
  
#' @rdname readSensorLocs
#' @param labels  Character vector indicating the labels of the sensors to be selected
#' @examples
#' 
#' locs <- getLocationsfromLabels(c('Fz', 'Cz', 'Pz'), plot = TRUE)
#' 
#' @export

getLocationsfromLabels <- function (labels, plot = FALSE) {
  readSensorLocs (source = 'EEGlocations', select = labels, plot = plot)
}

