#' @title NIMBLE Function to calculate the density of individuals alive in each habitat cell.
#'
#' @description
#' \code{calculateDensity} is a NIMBLE Function to calculate number of individual activity centers (s) in each habitat cell.
#' 
#' @param s \code{Vector} of x- and y-coordinates of individual AC location. 
#' @param habitatGrid \code{Matrix} Matrix of habitat window indices. 
#' @param indicator \code{Vector} specifying whether (1) or not (0) individuals are considered alive.
#' @param numWindows \code{Scalar} Number of habitat windows. 
#' @param nIndividuals \code{Scalar} Total number of individuals .
#'
#' @author Cyril Milleret
#'
#' @examples
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(rep(1,4))
#' logSumIntensity <- log(sum(c(1:4))) 
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' 
#' s <- matrix(NA,nrow=10,ncol=2)
#' for(i in 1:10){
#'   s[i,] <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'                      habitatGrid, numGridRows, numGridCols)
#' }
#' 
#' calculateDensity(s = s,
#'                  habitatGrid = habitatGrid,
#'                  indicator = rep(1, 10),
#'                  numWindows = prod(dim(habitatGrid)),
#'                  nIndividuals = 10
#' )
#' @rdname calculateDensity
#' @export
calculateDensity <- nimbleFunction(run = function(  s = double(2)
                                                  , habitatGrid = double(2)
                                                  , indicator = double(1)
                                                  , numWindows = double(0)
                                                  , nIndividuals = double(0)
){
  # Return type declaration
  returnType(double(1))
  
  dens <- numeric(length = numWindows, value = 0)
  for(i in 1:nIndividuals){
    if(indicator[i]==1){
      windowInd <- habitatGrid[trunc(s[i,2]) + 1, trunc(s[i,1]) + 1]
      dens[windowInd] <- 1 + dens[windowInd]
    }
  }
  
  return(dens)
  
})


