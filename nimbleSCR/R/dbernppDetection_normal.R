#' Bernoulli point process detection model
#' 
#' Density and random generation functions of the Bernoulli point process for detection.  
#' 
#' The \code{dbernppDetection_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' Bernoulli observations (\emph{x}) of a single individual in continuous space over a set of detection windows defined by their upper and lower
#' coordinates (\emph{lowerCoords,upperCoords}). The distribution assumes that an individual’s detection probability 
#' follows an isotropic multivariate normal centered on the individual's activity center (\emph{s}) with standard deviation (\emph{sd}).
#' 
#' @name dbernppDetection_normal
#' 
#' @param x Vector with three elements representing the x- and y-coordinates  and the corresponding id the detection window of a single spatial point (detection location).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (i.e. the AC location).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
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
#' coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
#'                                     1.5, 3.5,
#'                                     2.5, 3.5,
#'                                     3.5, 3.5,
#'                                     0.5, 2.5,
#'                                     1.5, 2.5,
#'                                     2.5, 2.5,
#'                                     3.5, 2.5,
#'                                     0.5, 1.5,
#'                                     1.5, 1.5,
#'                                     2.5, 1.5,
#'                                     3.5, 1.5,
#'                                     0.5, 0.5,
#'                                     1.5, 0.5,
#'                                     2.5, 0.5,
#'                                     3.5, 0.5), ncol = 2,byrow = TRUE)
#' colnames(coordsHabitatGridCenter) <- c("x","y")
#' # Create observation windows
#'  lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#'  upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#'  colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
#' # Rescale coordinates
#' ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData =  lowerCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords$coordsDataScaled[,2] <- ScaledUpperCoords$coordsDataScaled[,2] + 1.5
#' ScaledLowerCoords$coordsDataScaled[,2] <- ScaledLowerCoords$coordsDataScaled[,2] - 1.5
#' 
#' 
#' s <- c(1, 1)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' windowIndex <- 4
#' numPoints <- 1
#' numWindows <- 4
#' indicator <- 1
#' x <- c(0.5, 2)
#' windowIndex <- getWindowIndex(curCoords = x,
#'                               lowerCoords = ScaledLowerCoords$coordsDataScaled,
#'                               upperCoords =ScaledUpperCoords$coordsDataScaled)
#' x <- c(x, windowIndex)
#' 
#' dbernppDetection_normal(x, lowerCoords, upperCoords,
#'                         s, sd, baseIntensities
#'                         , numWindows,
#'                         indicator, log = TRUE)
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
    numWindows      = integer(0),
    indicator       = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## If the individual does not exists
    if(indicator == 0) {
      if(x[3] == 0){
        if(log) return(0.0) else return(1.0) 
      } else {
        if(log) return(-Inf) else return(0.0) 
      }
    }
    ## Make sure individuals with indicator==1 have one detection (necessary for dead recovery)
    if(indicator == 1) {
      if(x[3] == 0) {
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
    logPointIntensity <- 1.837877 + log(baseIntensities[x[3]]) + sum(dnorm((x[1:2] - s) / sd, log = 1))
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
    numWindows      = integer(0),
    indicator       = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernppDetection_normal only allows n = 1; using n = 1")
    }
    
    if(indicator==0){return(c(0,0,0))}else{
      
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
      
      windowIndex <- getWindowIndex(curCoords = outCoordinates[1, 1:2],
                                    lowerCoords = lowerCoords[1:numWindows,,drop = FALSE] ,
                                    upperCoords = upperCoords[1:numWindows,,drop = FALSE] )
      
      return(c(outCoordinates[1,],windowIndex))
    }
    returnType(double(1))
  }
)
