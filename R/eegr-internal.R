# eegr-internal.R - internal or barely commented
# functions not exported from the namespace
#
# Copyright (C) 2022  Geert van Boxtel,
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Version history
# 20220122  GvB       Initial setup, ckks()
#------------------------------------------------------------------------------

#' Internal functions not exported to the namespace
#'
#' @keywords internal
#' @noRd

# check if sensors are valid for
# x = data matrix, e.g. ctd
# s = sensors
# return value: TRUE or FALSE
chks <- function(x, s) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  if (!is.matrix(x)) {
    stop("x must be a matrix")
  }
  ns <- ncol(x)
  npts <- nrow(x)
  cn <- colnames(x)
  
  tf <- TRUE
  if (is.numeric(s)) {
    if (any(s <= 0 || s > ns)) {
      tf <- FALSE
    }
  } else if (is.character(s)) {
    if (length(s) <= 0) {
      tf <- FALSE
    } else if (length(s) == 1) {
      s <- unlist(strsplit(s, "[ ,;]"))
    }
    if (!(s == "avg" || s == "average" ||
          all(s %in% colnames(x)))) {
      tf <- FALSE
    }
  } else if (!is.na(s)) {
    tf <- FALSE
  }
  tf
}

# split sensor names specified as a single string
# into a character vector
splitsens <- function(s) {
  s <- unlist(strsplit(s, "[ ,;]"))
  s <- s[s != ""]
  s
}

# function that rounds half up (Matlab) instead of even (R)
# see https://theobligatescientist.blogspot.com/2010/02/r-i-still-love-you-but-i-hate-your.html
mround <- function (x) 
  ifelse(round(abs(x-trunc(x)), 1) == 0.5, trunc(x + 0.5), round(x))
