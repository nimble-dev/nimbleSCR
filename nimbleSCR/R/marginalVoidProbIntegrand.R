#' Integrand of the marginal void probability integral
#' 
#' Integrand of the marginal void probability integral. The domain of this function is the habitat domain.   
#'
#' @param x Matrix of x- and y-coordinates of a set of spatial points. One row corresponds to one point.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of detection windows. One row for each window.
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param numPoints Number of points that should be considered. This value (positive integer) is used to truncate \code{x} 
#' so that extra rows beyond \code{numPoints} are ignored. 
#' @param numWindows Number of windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored.
#' 
#' @return A vector of values of the integrand evaluated at each point of \code{x}.
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
#' x <- matrix(runif(10, 0, 2), nrow = 5)
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' numPoints <- 5
#' numWindows <- 4
#' marginalVoidProbIntegrand(x, lowerCoords, upperCoords, sd, baseIntensities, numPoints, numWindows)
#' 
#' @rdname marginalVoidProbIntegrand
#' @export
marginalVoidProbIntegrand <- nimbleFunction(
  run = function(
    x               = double(2), 
    lowerCoords     = double(2), 
    upperCoords     = double(2), 
    sd              = double(0), 
    baseIntensities = double(1),
    numPoints       = integer(0),
    numWindows      = integer(0)
  ) {
    ## numDims <- 2 ## Consider 2D models for now.
    constant <- 6.283185 * sd^2 ## Need this constant multiple times below
    expNumDetections <- rep(constant, numPoints)
    for(i in 1:numPoints){
      ## Compute the expected number of detections for the individual with AC location being x[i,]
      tmp <- 0
      for(j in 1:numWindows){
        ## See Eqn 28 of Zhang et al.(2020, DOI:10.1101/2020.10.06.325035)
        tmp <-  tmp + baseIntensities[j] * 
          prod(pnorm((upperCoords[j,] - x[i,]) / sd) - pnorm((lowerCoords[j,] - x[i,]) / sd))
      }
      expNumDetections[i] <- expNumDetections[i] * tmp
    }
    returnType(double(1))
    return(exp(-expNumDetections)) 
  }
)
