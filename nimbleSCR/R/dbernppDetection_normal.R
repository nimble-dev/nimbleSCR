#' Bernoulli point process detection model
#' 
#' Density and random generation functions of the Bernoulli point process for detection. 
#' An isotropic multivariate normal distribution is used as the decay kernel.
#' 
#' @name dbernppDetection_normal
#' 
#' @param x Vector of x- and y-coordinates of a single spatial point (detection location).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (i.e. the AC location).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param windowIndex Index of the detection window where \code{x} falls.
#' @param numPoints Number (0 or 1) of points that should be considered. 
#' @param numWindows Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored.
#' @param indicator Binary variable (0 or 1). \code{indicator = 0} means the individual is not available for detection 
#' and thus the probability of no detection is 1. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dbernppDetection_normal} gives the (log) probability density of the observation vector \code{x}. 
#' \code{rbernppDetection_normal} gives coordinates of a randomly generated spatial point.  
#' 
#' @author Wei Zhang and Cyril Milleret
#' 
#' @import nimble
#' @importFrom stats dnorm 
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
#' windowIndex <- 4
#' numPoints <- 1
#' numWindows <- 4
#' indicator <- 1
#' dbernppDetection_normal(c(1.1, 1.2), lowerCoords, upperCoords,
#'                         s, sd, baseIntensities, 
#'                         windowIndex, numPoints, numWindows,
#'                         indicator, log = TRUE)
#'
NULL
#' @rdname dbernppDetection_normal
#' @export
dbernppDetection_normal <- nimbleFunction(
  run = function(
    x               = double(1),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    windowIndex     = double(0),
    numPoints       = integer(0),
    numWindows      = integer(0),
    indicator       = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## If the individual does not exists
    if(indicator == 0) {
      if(numPoints == 0){
        if(log) return(0.0) else return(1.0) 
      } else {
        if(log) return(-Inf) else return(0.0) 
      }
    }
    ## Make sure individuals with indicator==1 have one detection (necessary for dead recovery)
    if(indicator == 1) {
      if(numPoints == 0) {
        if(log) return(-Inf) else return(0.0)
      }
    }
    ## Integrate the detection intensity over all detection windows
    windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE], 
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                   s = s,
                                                   baseIntensities = baseIntensities[1:numWindows], 
                                                   sd = sd,
                                                   numWindows =numWindows)
                                                 
    sumIntensity <- sum(windowIntensities) ## Expected total number of detections
    ## Make sure sumIntensity is positive 
    if(sumIntensity <= 0.0) {
      if(log) return(-Inf)
      else return(0.0)
    }
    # numDims <- 2 
    # logPointIntensity <- (numDims / 2.0) * log(2.0 * pi) + log(baseIntensities[windowIndex]) + sum(dnorm((x - s) / sd, log = 1))
    logPointIntensity <- 1.837877 + log(baseIntensities[windowIndex]) + sum(dnorm((x - s) / sd, log = 1))
    ## Log probability density
    logProb <- logPointIntensity - log(sumIntensity)
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dbernppDetection_normal
#' @export
rbernppDetection_normal <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    windowIndex     = double(0),
    numPoints       = integer(0),
    numWindows      = integer(0),
    indicator       = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernppDetection_normal only allows n = 1; using n = 1")
    }
    if(indicator==0){return(c(0,0))}else{
      
    ## Integrate the detection intensity over all detection windows
    windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                   s = s,
                                                   baseIntensities = baseIntensities[1:numWindows],
                                                   sd = sd,
                                                   numWindows = numWindows)
    ## Call the statified rejection sampler
    outCoordinates <- stratRejectionSampler_normal(numPoints = 1,
                                                   lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                   s = s,
                                                   windowIntensities = windowIntensities[1:numWindows],
                                                   sd = sd)
    return(outCoordinates[1,])
    }
    returnType(double(1))
  }
)
