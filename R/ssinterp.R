# ssinterp.R
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
# 20220119  GvB           Initial setup v0.3-2
# 20220208  GvB           bugfix in interpMx (is.na(dG), is.na(dH))
#---------------------------------------------------------------------------------------------------------------------

#' ssinterp
#' 
#' Spherical spline interpolation of sensor locations
#' 
#' @param src source sensor positions, specified as a data frame (or coerceable
#'   to a data frame), with 3 columns labeled x, y, and z.
#' @param dest destination sensor positions, specified as a data frame (or
#'   coerceable to a data frame) with 3 columns labeled x, y, and z.
#' @param lambda regularization parameter for smoothing the estimates, specified
#'   as a numeric scalar value. Default: 1e-5.
#' @param order order of the polynomial interpolation to use, specified as a
#'   positive integer scalar value. Default: 4
#' @param type character string specifying which interpolation method to use.
#'   One of \code{"spline"} - spherical spline interpolation) or \code{"slap"} -
#'   surface Laplacian (CSD). Default: \code{"spline"}
#' @param tol tolerance for the Legendre polynomial approximation, specified as a
#'   numeric scalar value. Default: 1e-7.
#' 
#' @return linear mapping matrix between old and new coordinates; a matrix with
#'   as many rows as the number of rows in \code{dest} and as many columns as
#'   the number of rows in \code{src}
#'
#' @examples
#' src <- matrix(runif(300), 100, 3)
#' src[, 3] <- abs(src[, 3])
#' dest <- matrix(runif(90), 30, 3)
#' dest[, 3] <- abs(dest[, 3])
#' W <- ssinterp(src, dest)
#' 
#' @author Matlab code 2009 by Jason D.R. Farquhar (\email{jdrf@@zepler.org}),
#'   adapted 2014 by Kay Robbins. R code by Geert van Boxtel,
#'   \email{G.J.M.vanBoxtel@@gmail.com}
#' 
#' @references Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F. (1989).
#'   Spherical splines for scalp potential and current density mapping.
#'   Electroencephalography and clinical neurophysiology, 72(2), 184-187.
#' 
#' @export

ssinterp <- function(src, dest, lambda = 1e-5, order = 4,
                     type = c("spline", "slap"), tol = 1e-7) {
  
  #parameter checking
  src <- t(as.matrix(src))    # procedure assumed x y z in rows
  dest <- t(as.matrix(dest))
  if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop('lambda must be a numeric scalar value')
  }
  if (!is.numeric(order) || length(order) != 1L ||
      order <= 0 || !order == trunc(order)) {
    stop('order must be a positive numeric scalar integer value')
  }
  type = match.arg(type)
  if (!is.numeric(tol) || length(tol) != 1L) {
    stop('lambda must be a numeric scalar value')
  }
  
  # Map the positions onto the unit sphere 
  src  <- src / pracma::repmat(sqrt(apply(src^2, 2, sum)), nrow(src), 1)
  dest <- dest / pracma::repmat(sqrt(apply(dest^2, 2, sum)), nrow(dest), 1)
  
  # Calculate the cosine of the angle between the new and old sensors. If
  # the vectors are on top of each other, the result is 1, if they are
  # pointing the other way, the result is -1
  cosSS <- t(src) %*% src  # angles between source positions
  cosDS <- t(dest) %*% src # angles between destination positions
  
  # Compute the interpolation matrix to tolerance tol
  Gss <- interpMx(cosSS, order, tol)$G
  GHds <- interpMx(cosDS, order, tol)
  Gds <- GHds$G
  
  # Include the regularisation
  if (lambda > 0) { 
    Gss <- Gss + lambda * diag(1, nrow(Gss))
  }
  
  # Compute the mapping to the polynomial coefficients space % [nSrc+1 x nSrc+1]
  # N.B. this can be numerically unstable so use the PINV to solve..
  muGss <- 1  # Used to improve condition number when inverting. Probably unnecessary
  C <- rbind(cbind(Gss, muGss * matrix(1, nrow(Gss), 1)),
             c(rep(muGss, ncol(Gss)), 0))
  # kludge
  #C[which(is.na(C))] <- 0
  #C[which(C == Inf)] <- 0
  
  iC <- pracma::pinv(C)

  # Compute the mapping from source measurements and positions to destination positions
  if (type == "spline") {
    W <- cbind(Gds, muGss * matrix(1, nrow(Gds), 1)) %*% iC[, 1:(ncol(iC) - 1)]
  } else if (pracma::strcmpi(type, 'slap')) {
    Hds <- GHds$H
    W <- Hds %*% iC[1:(nrow(iC) - 1), 1:(ncol(iC) - 1)]
  }
  W
}

interpMx <- function(cosEE, order, tol = 1e-10) {
  # compute the interpolation matrix for this set of point pairs
  nr <- nrow(cosEE)
  nc <- ncol(cosEE)
  numel <- nr * nc
  G <- H <- matrix(0, nr, nc)
  for (i in seq_len(numel)) {
    x <- cosEE[i]
    n <- 1; Pns1 <- 1; Pn <- x;        # seeds for the legendre ploy recurrence
    tmp  <- ((2 * n + 1) * Pn) / ((n * n + n)^order)
    G[i] <- tmp                        # 1st element in the sum
    H[i] <- (n * n + n) * tmp          # 1st element in the sum
    oGi <- Inf; dG <- abs(G[i]); oHi <- Inf; dH <- abs(H[i])
    for (n in 2:500) {                 # do the sum
      Pns2 <- Pns1; Pns1 <- Pn; 
      Pn <- ((2 * n - 1) * x * Pns1 - (n - 1) * Pns2) / n; # legendre poly recurrence
      oGi <- G[i];  oHi <- H[i];
      tmp <- ((2 * n + 1) * Pn) / ((n * n + n)^order)
      G[i] <- G[i] + tmp               # update function estimate, spline interp     
      H[i] <- H[i] + (n * n + n) * tmp # update function estimate, SLAP
      dG <- (abs(oGi - G[i]) + dG) / 2 
      dH <- (abs(oHi - H[i]) + dH) / 2 # moving ave gradient est for convergence
      if ((!is.na(dG) && dG < tol) && (!is.na(dH) && dH < tol)) {
        break
      }           # stop when tol reached
    }
  }
  G <- G / (4 * pi)
  H <- H / (4 * pi)
  
  list(G = G, H = H)
  
}
