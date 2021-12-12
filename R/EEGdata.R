# EEGdata.R
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
# 20201229  GvB       Initial setup (v0.3-0)
#---------------------------------------------------------------------------------------------------------------------

#' EEGdata
#'
#' Sample EEG, EOG, EMG, ECG and respiration data.
#' 
#' @author Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @examples
#' data(EEGdata)
#' plot(EEGdata)
#' 
#' @format A \code{ctd} object with 30 seconds of data sampled at 200 Hz
#'   containing the following signals, all measured in microvolts, except
#'   respiration:
#' \describe{
#'   \item{AF7 ... O2}{28 EEG channels, referenced to M1 (left mastoid)}
#'   \item{M2}{right mastoid,referenced to M1}
#'   \item{EOGh, EOGl, EOGr}{bipolar horizontal and vertical left and right EOG channels}
#'   \item{Resp}{Respiration measured by a chest belt (arbitrary units)}
#'   \item{ECG}{Electrocardiogram}
#'   \item{AgL, AgR}{Agonist and Antagonist EMG activity (finger flexion and extension)}
#' }
#'
"EEGdata"
