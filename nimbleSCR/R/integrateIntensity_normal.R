#' Integrate the multivariate normal intensity
#' 
#' Calculate the integral of the intensity function with an isotropic multivariate normal
#' kernel over a set of windows.
#' 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (AC location).
#' @param baseIntensities Vector of baseline intensities for all windows.
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param numWindows Total number of windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored. 

#' @return A vector of integrated intensities over all windows.
#' 
#' @author Wei Zhang
#' 
#' @import nimble
#' @importFrom stats pnorm
#' 
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' @examples 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' s <- c(1, 1)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' numWindows <- 4
#' integrateIntensity_normal(lowerCoords, upperCoords, s, baseIntensities, sd, numWindows)
#' 
#' @rdname integrateIntensity_normal
#' @export
integrateIntensity_normal <- nimbleFunction(
  run = function(
    lowerCoords      = double(2),
    upperCoords      = double(2),
    s                = double(1),
    baseIntensities  = double(1),
    sd               = double(0),
    numWindows       = double(0)
  ) {
    # numDims <- 2 ## We only consider 2D models for now
    # constant <- (2.0 * pi)^(numDims / 2.0) * (sd^numDims) ## A constant used repeatedly
    constant <- (6.283185)*(sd^2) 
    res <- constant * baseIntensities 
    for(i in 1:numWindows) {
      ## See Eqn 27 in Zhang et al. (2020, DOI:10.1101/2020.10.06.325035)
      res[i] <- res[i] * prod(pnorm((upperCoords[i,] - s) / sd) - pnorm((lowerCoords[i,] - s) / sd))
    }
    returnType(double(1))
    return(res)
  }
)
