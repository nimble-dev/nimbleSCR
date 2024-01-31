#' Bernoulli point process for the distribution of activity centers
#' 
#' Density and random generation functions of the Bernoulli point process for the distribution of activity centers. 
#'
#' The \code{dbernppAC} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' the activity center location (\emph{x}) of a single individual in continuous space over a set of habitat windows defined by their upper and lower
#' coordinates (\emph{lowerCoords,upperCoords}). The distribution assumes that the activity center  
#' follows a Bernoulli point process with intensity = \emph{exp(logIntensities)}.
#' 
#' 
#' @name dbernppAC 
#' 
#' @param x Vector of x- and y-coordinates of a single spatial point (i.e. AC location) scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}). One row for each window. 
#' Each window should be of size 1x1.
#' @param logIntensities Vector of log habitat intensities for all habitat windows. 
#' @param logSumIntensity Log of the sum of habitat intensities over all windows.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in
#' \code{lowerCoords}, \code{upperCoords}, and \code{logIntensities}. 
#' When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.
#' @param numGridRows,numGridCols Numbers of rows and columns of the habitat grid.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dbernppAC} gives the (log) probability density of the observation vector \code{x}. 
#' \code{rbernppAC} gives coordinates of a randomly generated spatial point.  
#' 
#' @author Wei Zhang
#' 
#' @import nimble
#' @importFrom stats runif
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' @examples
#' # Use the distribution in R
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(c(1:4))
#' logSumIntensity <- log(sum(c(1:4)))  
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' dbernppAC(c(0.5, 1.5), lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'           habitatGrid, numGridRows, numGridCols, log = TRUE)

NULL
#' @rdname dbernppAC
#' @export
dbernppAC <- nimbleFunction(
  run = function(
    x               = double(1),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## Check if the point falls within the habitat: 
    ## Getting numGridRows and numGridCols using the following code takes some time
    ## and may cause inefficiency if the function is called repeatedly in a loop.
    ## numGridRows <- dim(habitatGrid)[1]
    ## numGridCols <- dim(habitatGrid)[2]
    ## In addition, when the true habitat grid has one row/column we need to inflate it for use in NIMBLE model code. 
    ## In this case, we have problems using the code above.
    ## Note that we need to rescale the habitat gird to ensure x and y coordinates start from 0
    ## and each window is of size 1x1. So the following code works correctly. 
    if(min(x) < 0 | x[2] >= numGridRows | x[1] >= numGridCols) {
      if(log) return(-Inf) 
      else return(0.0)
    }
    ## Find which window point x falls within
    windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1]
    ## windowInd == 0 means this window is not defined as habitat
    if(windowInd == 0) {
        if(log) return(-Inf)
        else return(0.0)
    }
    ## Log probability density 
    logProb <- logIntensities[windowInd] - logSumIntensity
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dbernppAC
#' @export
rbernppAC <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0)
  ) {
    if(n <= 0) stop("The number of requested samples must be above zero")
    else if(n > 1) print("rbernppAC only allows n = 1; using n = 1")
    ## Simulate window index
    windowInd <- rcat(1, exp(logIntensities))
    numDims <- 2 
    ## A uniform distribution is used within the target window
    outCoordinates <- lowerCoords[windowInd,] + 
      runif(numDims, 0.0, 1.0) * (upperCoords[windowInd,] - lowerCoords[windowInd,])
    return(outCoordinates)
    returnType(double(1)) 
  }
)


