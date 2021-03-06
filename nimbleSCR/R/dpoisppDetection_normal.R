#' Poisson point process detection model
#' 
#' Density and random generation functions of the Poisson point process for detection. 
#' An isotropic multivariate normal distribution is used as the decay kernel.
#' 
#' @name dpoisppDetection_normal
#' 
#' @param x Matrix of x- and y-coordinates of a set of spatial points (detection locations). One row corresponds to one point.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (the AC location).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param windowIndices Vector of indices of detection windows where each point of \code{x} falls. 
#' @param numPoints Number of points that should be considered. For \code{dpoisppDetection_normal}, this should be a non-negative integer and is used to truncate \code{x} 
#' so that extra rows beyond \code{numPoints} are ignored. For \code{rpoisppDetection_normal}, this can be any integer. 
#' Usually it is set to be a negative value, which indicates that the number of points is generated from a Poisson model. 
#' Otherwise, the number of points generated will be equal to the specified value.
#' @param numWindows Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored.  
#' @param indicator Binary variable (0 or 1) used for data augmentation. \code{indicator = 0} means the individual does not exist 
#' and thus the probability of no detection is 1.  
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dpoisppDetection_normal} gives the (log) probability density of the observation matrix \code{x}.
#' \code{rpoisppDetection_normal} gives coordinates of a set of randomly generated spatial points.  
#' 
#' @author Wei Zhang 
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
#' x <- matrix(c(0.5, 0.2, 1.1, 0.4, 1.4, 0.3, 0.1, 1.3, 1, 1.5), nrow = 5, byrow = TRUE)
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' s <- c(1, 1)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' windowIndices <- c(1, 2, 2, 3, 4)
#' numPoints <- 5
#' numWindows <- 4
#' indicator <- 1
#' dpoisppDetection_normal(x, lowerCoords, upperCoords, s, sd, baseIntensities,
#'                         windowIndices, numPoints, numWindows, indicator, log = TRUE)

NULL
#' @rdname dpoisppDetection_normal
#' @export
dpoisppDetection_normal <- nimbleFunction(
  run = function(
    x               = double(2),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    windowIndices   = double(1),  
    numPoints       = integer(0),
    numWindows      = integer(0),
    indicator       = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## If the individual does not exists 
    if(indicator == 0) {
      if(numPoints == 0){
        if(log) return(0.0)
        else return(1.0)
      } else{
        if(log) return(-Inf)
        else return(0.0) 
      }
    }
    ## Integrate the detection intensity function over all detection windows
    windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE], 
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                   s = s,
                                                   baseIntensities = baseIntensities[1:numWindows], 
                                                   sd = sd,
                                                   numWindows = numWindows)
    sumIntensity <- sum(windowIntensities) ## Expected total number of detections
    ## Calculate the probability
    if(numPoints == 0) {
      logProb <- -sumIntensity
    } else{
      ## 1.837877 <- log(2.0 * pi) 
      logPointIntensity <- rep(1.837877, numPoints)
      for(i in 1:numPoints) {
        ## Log intensity at the ith point: see Eqn 24 of Zhang et al.(2020, DOI:10.1101/2020.10.06.325035)
        pointBaseIntensity <- baseIntensities[windowIndices[i]]
        logPointIntensity[i] <- logPointIntensity[i] + log(pointBaseIntensity) + sum(dnorm((x[i,] - s) / sd, log = 1))
      }
      ## Log probability density
      logProb <- sum(logPointIntensity) - sumIntensity
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dpoisppDetection_normal
#' @export
rpoisppDetection_normal <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    windowIndices   = double(1),  
    numPoints       = integer(0),
    numWindows      = integer(0),
    indicator       = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rpoisppDetection_normal only allows n = 1; using n = 1")
    }
    numDims <- 2 
    ## If the individual does not exist in the population, we return an empty matrix directly
    if(indicator == 0) return(matrix(nrow = numPoints, ncol = numDims))
    
    ## Initialize an output matrix
    ## This empty matrix will be returned if numWindows==0 or numPoints==0
    outCoordinates <- matrix(nrow = numPoints, ncol = numDims)

    ## Generate random points below:
    ## Current implementation uses a stratified rejection sampler. 
    ## This can become inefficient if the observation window becomes very large compared to the 
    ## standard deviation of the multivariate normal distribution. 
    if((numWindows > 0) & (numPoints != 0)) {
      ## Integrate the detection intensity function over all detection windows
      windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                     upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                     s = s,
                                                     baseIntensities = baseIntensities[1:numWindows],
                                                     sd = sd,
                                                     numWindows = numWindows)
                                                  
      sumIntensity <- sum(windowIntensities)
      if(sumIntensity > 0.0) {
       
        numPoints1 <- rpois(1, sumIntensity)
        if(numPoints1>numPoints){stop("There are more simulated individual detections than what can be stored within x.\n You may need to increase the size of the 'x' object and make sure that 'numPoints'is equal to dim(x)[2]" )}
        if(numPoints1 > 0) {
          outCoordinates[1:numPoints1, ] <- stratRejectionSampler_normal(numPoints = numPoints1,
                                                         lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                         upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                         s = s,
                                                         windowIntensities = windowIntensities[1:numWindows],
                                                         sd = sd)
        }
      }
    }
    return(outCoordinates)
    returnType(double(2))
  }
)



