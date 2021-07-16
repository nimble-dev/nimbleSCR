#' Bernoulli point process for activity center movement
#' 
#' Density and random generation functions of the Bernoulli point process for activity center movement. Intended for use in open population models. 
#' An isotropic multivariate normal distribution is used as the decay kernel. 
#' 
#' @name dbernppACmovement_normal
#' 
#' @param x Vector of x- and y-coordinates of a single spatial point (typically AC location at time t+1). 
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows. One row for each window. 
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (AC location at time t).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline habitat intensities for all habitat windows.
#' @param habitatGrid Matrix of habitat window indices. Baseline habitat intensities should match the order in \code{lowerCoords} and \code{upperCoords}. 
#' When the grid has only one row/column, artificial indices have to be provided to inflate \code{habitatGrid} in order to be able to use the distribution in \code{nimble} model code.     
#' @param numGridRows,numGridCols Numbers of rows and columns of the habitat grid. 
#' @param numWindows Number of habitat windows. This value (positive integer) can be used to truncate \code{lowerCoords} and \code{upperCoords} so that extra rows beyond \code{numWindows} are ignored.  
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dbernppACmovement_normal} gives the (log) probability density of the observation vector \code{x}. 
#' \code{rbernppACmovement_normal} gives coordinates of a randomly generated spatial point.  
#' 
#' @author Wei Zhang and Cyril Milleret
#'
#' @import nimble
#' @importFrom stats dnorm  
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' 
#' @examples 
#' # Use the distribution in R
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' s <- c(1, 1) # Currrent activity center location
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numRows <- nrow(habitatGrid)
#' numCols <- ncol(habitatGrid)
#' numWindows <- 4
#' # The log probability density of moving from (1,1) to (1.2, 0.8) 
#' dbernppACmovement_normal(c(1.2, 0.8), lowerCoords, upperCoords, s, sd, baseIntensities, 
#'                          habitatGrid, numRows, numCols, numWindows, log = TRUE)

NULL
#' @rdname dbernppACmovement_normal
#' @export  
dbernppACmovement_normal <- nimbleFunction(
  run = function(
    x               = double(1),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0),
    numWindows      = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## Check if the point x falls within the habitat
    if(min(x) < 0 | x[2] >= numGridRows | x[1] >= numGridCols) {
      if(log) return(-Inf) 
      else return(0.0)
    }
    ## Index of the window where x falls
    windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1]
    ## windowInd == 0 means this window is not defined as habitat
    if(windowInd == 0) {
        if(log) return(-Inf)
        else return(0.0)
    }
    ## Integrate the intensity function over all habitat windows
    windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE], 
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                   s = s,
                                                   baseIntensities = baseIntensities[1:numWindows], 
                                                   sd = sd,
                                                   numWindows = numWindows)
    sumIntensity <- sum(windowIntensities)
    ## Log intensity at x
    # numDims <- 2
    # logPointIntensity <- (numDims / 2.0) * log(2.0 * pi) + log(baseIntensities[windowInd]) + sum(dnorm((x - s) / sd, log = 1))
    logPointIntensity <-  1.837877 + log(baseIntensities[windowInd]) + sum(dnorm((x - s) / sd, log = 1))
    
    ## Log probability density under the Bernoulli point process
    logProb <- logPointIntensity - log(sumIntensity)
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dbernppACmovement_normal
#' @export
rbernppACmovement_normal <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0),
    numWindows      = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernppACmovement only allows n = 1; using n = 1")
    }
    ## Integrate the intensity function over all habitat windows
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
    returnType(double(1))
  }
)



