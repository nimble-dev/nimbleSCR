#' Marginal void probability
#' 
#' Calculate the marginal void probability using the midpoint integration method. 
#' 
#' @param quadNodes Three-dimensional array of nodes for midpoint integration. The dimension sizes are equal to 
#' the number of nodes per habitat window (1st), 2 (2nd), and the number of habitat windows (3rd).
#' @param quadWeights Vector of weights for midpoint integration. 
#' @param numNodes Vector of numbers of nodes for all habitat windows.
#' @param lowerCoords,upperCoords Matrix of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param habIntensities Vector of habitat intensities for all habitat windows.
#' @param sumHabIntensity Total habitat selection intensity over all windows.
#' @param numObsWindows Number of detection windows.
#' @param numHabWindows Number of habitat windows.
#'
#' @return The marginal void probability. 
#' @author Wei Zhang
#' 
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' @examples 
#' lowerHabCoords <- matrix(c(0, 0, 0, 1), nrow = 2, byrow = TRUE)
#' upperHabCoords <- matrix(c(2, 1, 2, 2), nrow = 2, byrow = TRUE)
#' lowerObsCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperObsCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' nodesRes <- getMidPointNodes(lowerHabCoords, upperHabCoords, 10)
#' quadNodes <- nodesRes$quadNodes
#' quadWeights <- nodesRes$quadWeights
#' numNodes <- rep(100, 2)
#' sd <- 0.1
#' baseDetIntensities <- c(1:4)
#' habIntensities <- c(1:2)
#' sumHabIntensity <- sum(habIntensities * c(2, 2))
#' numObsWindows <- 4
#' numHabWindows <- 2
#' marginalVoidProbNumIntegration(quadNodes, quadWeights, numNodes,
#'                                lowerObsCoords, upperObsCoords, sd, 
#'                                baseDetIntensities, habIntensities,
#'                                sumHabIntensity, numObsWindows, numHabWindows)
#'
#' @rdname marginalVoidProbNumIntegration
#' @export
marginalVoidProbNumIntegration <- nimbleFunction(
  run = function(
    quadNodes       = double(3),
    quadWeights     = double(1),
    numNodes        = double(1),
    lowerCoords     = double(2), 
    upperCoords     = double(2), 
    sd              = double(0), 
    baseIntensities = double(1),
    habIntensities  = double(1),
    sumHabIntensity = double(0),
    numObsWindows   = integer(0),
    numHabWindows   = integer(0)
  ) {
    tmp <- rep(0.0, numHabWindows)
    for(i in 1:numHabWindows){
      ## Evaluate the integrand at all nodes in the ith habitat window 
      integrandVals <- marginalVoidProbIntegrand(x = quadNodes[,,i],
                                                 lowerCoords = lowerCoords[1:numObsWindows,,drop = FALSE],
                                                 upperCoords = upperCoords[1:numObsWindows,,drop = FALSE],
                                                 sd = sd,
                                                 baseIntensities = baseIntensities[1:numObsWindows],
                                                 numPoints = numNodes[i],
                                                 numWindows = numObsWindows)
      tmp[i] <- habIntensities[i] * quadWeights[i] * sum(integrandVals)
    }
    outProb <- sum(tmp) / sumHabIntensity
    return(outProb)
    returnType(double(0))
  }
)
