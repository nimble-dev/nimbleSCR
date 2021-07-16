#' Get window index
#'
#' From a set of windows, find the index of the window into which a given point falls. 
#' Can be applied to detection and habitat windows.
#'
#' @param curCoords Vector of coordinates of a single spatial point
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' 
#' @return Index of the window where the given point falls.
#'
#' @author Pierre Dupont
#'
#' @examples
#' 
#' sourceCoords <- c(1.5,2.2)
#' lowerCoords <- cbind(c(0,1,3,0),c(0,1,2,2))
#' upperCoords <- cbind(c(1,3,5,3),c(1,2,4,4))
#' getWindowIndex(sourceCoords, lowerCoords, upperCoords)
#' 
#' @rdname getWindowIndex
#' @export
getWindowIndex <- nimbleFunction(
  run = function(
    curCoords   = double(1),
    lowerCoords = double(2),
    upperCoords = double(2)
  ) {
    numDims <- dim(lowerCoords)[2]
    numWindows <- dim(lowerCoords)[1]
    chosenWindow <- -1
    for(winIter in 1:numWindows) {
      test <- 1
      for(dimIter in 1:numDims) {
        tmp <- (lowerCoords[winIter, dimIter] <= curCoords[dimIter]) * (upperCoords[winIter, dimIter] > curCoords[dimIter])
        test <- test * tmp
      }
      if(test == 1) {
        chosenWindow <- winIter
        return(chosenWindow)
      }
    }
    return(chosenWindow) ## -1 means the point does not fall within any of the given windows
    returnType(double(0))
  })
