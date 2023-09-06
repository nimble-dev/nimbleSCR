#' @title NIMBLE function to calculate the density of individuals alive in each habitat cell.
#'
#' @description
#' \code{calculateDensity} is a NIMBLE function to calculate number of individual activity centers (s) in each habitat cell.
#' 
#' @param s \code{Matrix} of x- and y-coordinates of individual AC locations. 
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the 
#' order of habitat windows in spatial probabilities (e.g. \code{prob1To2Hab} as used in the function \code{dcatState1Alive2Dead}) or in \code{lowerCoords} and \code{upperCoords} as used in the \code{dbernppAC} function. #' @param indicator \code{Vector} of binary arguments specifying whether the individuals are considered alive (indicator = 1) or not (indicator = 0).
#' @param indicator \code{Vector} of binary arguments specifying whether the individuals are considered alive (indicator = 1) or not (indicator = 0).
#' @param numWindows \code{Scalar} Number of habitat windows. 
#' @param nIndividuals \code{Scalar} Number of individuals.
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
  numGridRows <- dim(habitatGrid)[1]
  for(i in 1:nIndividuals){
    if(indicator[i]==1){
      windowInd <- habitatGrid[numGridRows-trunc(s[i,2]), trunc(s[i,1]) + 1]
      dens[windowInd] <- 1 + dens[windowInd]
    }
  }
  return(dens)
  
})


