#' Integrate the multivariate exponential intensity with local evaluation
#' 
#' Calculate the integral of the intensity function with an isotropic multivariate exponential
#' kernel over a set of windows. The local evaluation technique is implemented. 
#' 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate exponential distribution mean (AC location).
#' @param baseIntensities Vector of baseline intensities for all windows.
#' @param lambda Rate parameter of the isotropic multivariate exponential distribution.
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
#' @examples 
#' 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' s <- c(0.1, 0.9)
#' lambda <- 0.1
#' baseIntensities <- c(1:4)
#' numLocalWindows <- 2
#' localWindows <- c(1, 3)
#' integrateIntensityLocal_exp(lowerCoords, upperCoords,
#'                             s, baseIntensities, lambda, 
#'                             numLocalWindows, localWindows)
#' 
#' @keywords internal
#' @export
integrateIntensityLocal_exp <- nimbleFunction(
  run = function(
    lowerCoords      = double(2),
    upperCoords      = double(2),
    s                = double(1),
    baseIntensities  = double(1),
    lambda           = double(0),
    numLocalWindows  = double(0),
    localWindows     = double(1)
  ) {
    
    res <- rep(0.0, numLocalWindows)
    for(i in 1:numLocalWindows) { 
      res[i] <- baseIntensities[localWindows[i]] / lambda^2
      for(j in 1:2){
        ## Lower and upper bounds for integration
        uppIntBound <- upperCoords[localWindows[i],j] - s[j]
        lowIntBound <- lowerCoords[localWindows[i],j] - s[j] 
        ## Need to ensure upperCoords > lowerCoords for input data
        if(uppIntBound*lowIntBound > 0)
          res[i] <- res[i] * abs(exp(-lambda*abs(lowIntBound)) - exp(-lambda*abs(uppIntBound)))
        else
          res[i] <- res[i] * (2 - exp(lambda*lowIntBound) - exp(-lambda*uppIntBound)) / lambda
      }
    }
    returnType(double(1))
    return(res)
    
    
  }
)