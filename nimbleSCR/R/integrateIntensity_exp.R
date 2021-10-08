#' Integrate the multivariate exponential intensity
#' 
#' Calculate the integral of the intensity function with an isotropic multivariate exponential
#' kernel over a set of windows.
#' 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the AC location.
#' @param baseIntensities Vector of baseline intensities for all windows.
#' @param lambda Rate parameter of the isotropic multivariate exponential distribution.
#' @param numWindows Total number of windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored. 

#' @return A vector of integrated intensities over all windows.
#' 
#' @author Wei Zhang
#' 
#' @import nimble
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
#' lambda <- 1.0
#' baseIntensities <- c(1:4)
#' numWindows <- 4
#' integrateIntensity_exp(lowerCoords, upperCoords, s, baseIntensities, lambda, numWindows)
#' 
#' @export
integrateIntensity_exp <- nimbleFunction(
  run = function(
    lowerCoords      = double(2),
    upperCoords      = double(2),
    s                = double(1),
    baseIntensities  = double(1),
    lambda           = double(0),
    numWindows       = double(0)
  ) {
    res <- baseIntensities 
    for(i in 1:numWindows) { 
      for(j in 1:2){
        ## Lower and upper bounds for integration
        uppIntBound <- upperCoords[i,j] - s[j]
        lowIntBound <- lowerCoords[i,j] - s[j] 
        ## Efficiency is not considered here
        if(uppIntBound > 0 & lowIntBound > 0) 
          res[i] <- res[i] * (exp(-lambda*lowIntBound) - exp(-lambda*uppIntBound)) / lambda
        else if(uppIntBound < 0 & lowIntBound < 0) 
          res[i] <- res[i] * (exp(lambda*uppIntBound) - exp(lambda*lowIntBound)) / lambda
        else
          res[i] <- res[i] * (2 - exp(lambda*lowIntBound) - exp(-lambda*uppIntBound)) / lambda
      }
    }
    returnType(double(1))
    return(res)
  }
)

