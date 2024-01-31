#' Poisson point process for the distribution of activity centers
#' 
#' Density and random generation functions of the Poisson point process for the distribution of activity centers. 
#
#' The \code{dpoisppAC} distribution is a NIMBLE custom distribution which 
#' can be used to model and simulate activity center locations (\emph{x}) of multiple individual in 
#' continuous space over a set of habitat windows defined by their upper and lower coordinates (\emph{lowerCoords,upperCoords}). 
#' The distribution assumes that activity centers follow a Poisson point process with intensity = exp(logIntensities).  
#' All coordinates (\emph{s} and \emph{trapCoords}) should be scaled to the habitat (\code{\link{scaleCoordsToHabitatGrid}}).
#' 
#' @name dpoisppAC
#' 
#' @param x Matrix of x- and y-coordinates of a set of spatial points (AC locations) scaled to the habitat (\code{\link{scaleCoordsToHabitatGrid}}). Each row corresponds to a point.
#' @param n Integer specifying the number of realisations to generate. Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}). One row for each window. Each window should be of size 1x1.
#' @param logIntensities Vector of log habitat intensities for all habitat windows. 
#' @param sumIntensity Sum of the habitat intensities over all windows. Provided as an argument for computational speed, instead of calculating it in the function. 
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in
#' \code{lowerCoords}, \code{upperCoords}, and \code{logIntensities}. 
#' When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.
#' @param numGridRows,numGridCols Numbers of rows and columns of the habitat grid. 
#' @param numPoints Number of points in the Poisson point process. This value (non-negative integer) is used to truncate \code{x} 
#' so that extra rows beyond \code{numPoints} are ignored. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dpoisppAC} gives the (log) probability density of the observation matrix \code{x}. 
#' \code{rpoisppAC} gives coordinates of a set of randomly generated spatial points.  
#' 
#' @author Wei Zhang
#' @import nimble
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#'
#' @examples 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' logIntensities <- log(c(1:4))
#' logSumIntensity <- sum(exp(logIntensities))
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' #Simulate data
#' x <- rpoisppAC(1, lowerCoords, upperCoords, logIntensities, logSumIntensity, habitatGrid,
#'                numGridRows, numGridCols, -1)
#' numPoints <- nrow(x)
#' dpoisppAC(x, lowerCoords, upperCoords, logIntensities, logSumIntensity,
#'           habitatGrid, numGridRows, numGridCols, numPoints, log = TRUE)
#'           
NULL
#' @rdname dpoisppAC
#' @export
dpoisppAC <- nimbleFunction(
  run = function(
    x              = double(2),
    lowerCoords    = double(2),
    upperCoords    = double(2),
    logIntensities = double(1),
    sumIntensity   = double(0),
    habitatGrid    = double(2),
    numGridRows    = integer(0),
    numGridCols    = integer(0),
    numPoints      = integer(0),
    log            = integer(0, default = 0)
  ) {
    if(numPoints == 0) {
      ## No points in the Poisson process
      logProb <- -sumIntensity
    } else { 
      ## Check if all points fall within the habitat 
      if(min(x[1:numPoints,]) < 0 | max(x[1:numPoints, 2]) >= numGridRows | max(x[1:numPoints, 1]) >= numGridCols) {
        if(log) return(-Inf) 
        else return(0.0)
      }
      logProb <- -sumIntensity
      for(i in 1:numPoints) {
        windowInd <- habitatGrid[trunc(x[i,2])+1, trunc(x[i,1])+1] ## Window index of the ith point
        ## windowInd == 0 means this window is not defined as habitat
        if(windowInd == 0) {
            if(log) return(-Inf)
            else return(0.0)
        }
        logProb <- logProb + logIntensities[windowInd]
      }
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dpoisppAC
#' @export
rpoisppAC <- nimbleFunction(
  run = function(
    n              = integer(0),
    lowerCoords    = double(2),
    upperCoords    = double(2),
    logIntensities = double(1),
    sumIntensity   = double(0),
    habitatGrid    = double(2),
    numGridRows    = integer(0),
    numGridCols    = integer(0),
    numPoints      = integer(0)
  ) {
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rpoisppAC only allows n = 1; using n = 1")
    }
    numDims <- 2 
    ## If numPoints is negative, generate it from a Poisson model
    if(numPoints < 0) numPoints <- rpois(1, sumIntensity)
    outCoordinates <- matrix(nrow = numPoints, ncol = numDims)
    if(numPoints == 0) return(outCoordinates) 
    ## When numPoints > 0:
    for(i in 1:numPoints) {
      windowInd <- rcat(1, exp(logIntensities))
      ## Within each window, it follows a uniform distribution
      outCoordinates[i,] <- lowerCoords[windowInd,] + 
        runif(numDims, 0.0, 1.0) * (upperCoords[windowInd,] - lowerCoords[windowInd,])
    }
    return(outCoordinates)
    returnType(double(2))
  }
)



