#' Integrate the multivariate normal intensity with local evaluation
#' 
#' Calculate the integral of the intensity function with an isotropic multivariate normal
#' kernel over a set of windows. The local evaluation technique is implemented. 
#' 
#' 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (AC location).
#' @param baseIntensities Vector of baseline intensities for all windows.
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param numLocalWindows Number of windows that are close to the activity center
#' @param localWindows Vector of indices of the windows that are close to the activity center.
#'
#' @return A vector of integrated intensities over all local windows.
#' 
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#'
#' @author Cyril Milleret and Wei Zhang
#' 
#' @import nimble
#' @importFrom stats pnorm
#' 
#' @examples 
#' 
#'  lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#'  upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#'  s <- c(0.1, 0.9)
#'  sd <- 0.1
#'  baseIntensities <- c(1:4)
#'  numLocalWindows <- 2
#'  localWindows <- c(1, 3)
#'  integrateIntensityLocal_normal(lowerCoords, upperCoords, s,
#'                                 baseIntensities, sd,
#'                                 numLocalWindows, localWindows)
#' 
#' @export
integrateIntensityLocal_normal <- nimbleFunction(
  run = function(
    lowerCoords      = double(2),
    upperCoords      = double(2),
    s                = double(1),
    baseIntensities  = double(1),
    sd               = double(0),
    numLocalWindows  = double(0),
    localWindows     = double(1)
  ) {
    # numDims <- 2 ## We only consider 2D models for now
    res <- rep(0.0, numLocalWindows)
    constant <- 6.283185 * sd^2  # (2.0 * pi)^(numDims / 2.0) * (sd^numDims)
    for(i in 1:numLocalWindows) {
      res[i] <- constant * baseIntensities[localWindows[i]] *
        prod(pnorm((upperCoords[localWindows[i],] - s) / sd) - pnorm((lowerCoords[localWindows[i],] - s) / sd))
    }
    returnType(double(1))
    return(res)
  }
)