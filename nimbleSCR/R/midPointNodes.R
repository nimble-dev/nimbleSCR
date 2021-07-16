#' Generate midpoint integration nodes 
#' 
#' Generate midpoint nodes and weights for integrating a function numerically over a set of windows. 
#' For each window, generate a set of equally spaced nodes and weights.
#' 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of windows. One row for each window.
#' @param numSubintervals Number of subintervals each dimension of a window is divided into.
#' 
#' @return A list of midpoint nodes and weights.
#' @author Wei Zhang 
#' 
#' @examples 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' midPointNodes(lowerCoords, upperCoords, 5)
#' 
#' @rdname midPointNodes
#' @export
midPointNodes <- function(lowerCoords, upperCoords, numSubintervals = 10) {
  numDims <- ncol(lowerCoords)
  numWindows <- nrow(lowerCoords)
  numNodes <- numSubintervals^numDims ## Number of nodes for each window
  quadNodes <- array(0, dim = c(numNodes, numDims, numWindows))
  quadWeights <- rep(0, numWindows) ## All nodes in the same window have the same weight
  allCoords <- list()
  for(i in 1:numWindows){
    loCoords <- lowerCoords[i,]
    upCoords <- upperCoords[i,]
    cellWidths <- (upCoords - loCoords) / numSubintervals
    quadWeights[i] <- prod(cellWidths) ## Weight is the area of each cell
    for(j in 1:numDims){
      tmp <- 0.5 * cellWidths[j]
      allCoords[[j]] <- seq(loCoords[j] + tmp, upCoords[j] - tmp, length = numSubintervals)
    }
    quadNodes[,,i] <- as.matrix(expand.grid(allCoords))
  }
  ## Return an array of midpoints and a vector of weights
  return(list(quadNodes = quadNodes, quadWeights = quadWeights))
}
