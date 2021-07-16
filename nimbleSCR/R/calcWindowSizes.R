#' Window size calculation
#' 
#' Calculates the sizes of a set of windows based on their lower and upper coordinates of each dimension.
#' Can be applied to detection and habitat windows.
#'  
#' @param lowerCoords Matrix of lower x- and y-coordinates of all windows. One row for each window.
#' @param upperCoords Matrix of upper x- and y-coordinates of all windows. One row for each window.
#' 
#' @return A vector of window sizes.
#' 
#' @author Wei Zhang
#' 
#' @examples
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 3, 1, 1, 4, 3, 4), nrow = 4, byrow = TRUE)
#' calcWindowSizes(lowerCoords, upperCoords)
#' 
#'
#' @rdname calcWindowSizes
#' @export
calcWindowSizes <- nimbleFunction(
  run = function(
    lowerCoords = double(2),
    upperCoords = double(2)
  ) {
    numDims <- 2 ## We consider 2D models for now
    numWindows <- dim(lowerCoords)[1]
    sizes <- rep(1.0, numWindows)
    for(i in 1:numDims) {
      ## Writing the code below into two lines avoids compilation problems
      tmp <-  upperCoords[,i] - lowerCoords[,i]
      sizes <- sizes * tmp
    }
    returnType(double(1))
    return(sizes)
  }
)
