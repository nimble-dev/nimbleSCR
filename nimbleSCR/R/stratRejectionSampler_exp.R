#' Stratified rejection sampler for multivariate exponential point process
#' 
#' Simulate data using a stratified rejection sampler from a point process with an isotropic multivariate exponential decay kernel.  
#' 
#' @param numPoints Number of spatial points to generate. 
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of a set of detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of of the isotropic multivariate exponential distribution mean.
#' @param windowIntensities Vector of integrated intensities over all detection windows.
#' @param lambda Rate parameter of the isotropic multivariate exponential distribution.
#' 
#' @return A matrix of x- and y-coordinates of the generated points. One row corresponds to one point.
#' 
#' @import nimble
#' @importFrom stats qexp runif
#' @author Wei Zhang
#' 
#' @examples 
#' numPoints <- 10
#' lowerObsCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperObsCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' s <- c(1, 1)
#' windowIntensities <- c(1:4)
#' lambda <- 0.1
#' stratRejectionSampler_exp(numPoints, lowerObsCoords, upperObsCoords, s, windowIntensities, lambda)
#' 
#' @rdname stratRejectionSampler_exp
#' @export
stratRejectionSampler_exp <- nimbleFunction(
  run = function(
    numPoints         = integer(0),
    lowerCoords       = double(2),
    upperCoords       = double(2),
    s                 = double(1),
    windowIntensities = double(1),
    lambda            = double(0)
  ) {
    numWindows <- dim(lowerCoords)[1]
    numDims <- 2 
    if(numPoints <= 0) return(matrix(nrow = 0, ncol = numDims))
    sumIntensity <- sum(windowIntensities)
    outCoordinates <- matrix(nrow = numPoints, ncol = numDims)
    ## Test the edge case where the sum of the intensities is zero
    ## This seems unnecessary, but not a big deal
    if(sumIntensity <= 0.0) {
      windowIntensities <- calcWindowSizes(lowerCoords, upperCoords)
      sumIntensity <- sum(windowIntensities)
      if(sumIntensity <= 0.0) return(matrix(nrow = 0, ncol = numDims))
    }
    ## Find the limits at which the max truncation probability is obtained
    maxWidth <- qexp(0.999, rate = lambda)
    ## Ensure that cells further than the designated maximum width are removed as candidates for sampling points from
    withinMaxWidth <- numeric(length = numWindows)
    for(i in 1:numWindows) {
      withinMaxWidth[i] <- prod(s[1:numDims] - maxWidth < upperCoords[i, 1:numDims]) * 
        prod(s[1:numDims] + maxWidth > lowerCoords[i, 1:numDims])
    }
    correctedIntensities <- windowIntensities[1:numWindows] * withinMaxWidth[1:numWindows]
    for(i in 1:numPoints) {
      ## Randomly sample an index from a categorical distribution weighted according to 
      ## the integral of the intensity function
      curInd <- rcat(1, correctedIntensities[1:numWindows])
      ## Calculate the nearest point to the source coordinate in the chosen detection window
      nearestPoint <- pmin(pmax(s[1:numDims], lowerCoords[curInd, 1:numDims]), upperCoords[curInd, 1:numDims])
      ## Calculate the intensity at the nearest point (omitting a baseline detection intensity value)
      ## (this is the maximum possible intensity in the selected detection window)
      maxIntensity <- exp(-lambda * sum(abs(nearestPoint[1:numDims] - s[1:numDims])))
      ## Create a set of test coordinates (truncating the uniform distribution when it is too far from 
      ## the source coordinates to avoid inefficient rejection sampling)
      testCoords <- runif(numDims,
                          min = pmax(lowerCoords[curInd, 1:numDims], s[1:numDims] - maxWidth),
                          max = pmin(upperCoords[curInd, 1:numDims], s[1:numDims] + maxWidth))
      ## Calculate the intensity at the test coordinates
      testIntensity <- exp(-lambda * sum(abs(testCoords[1:numDims] - s[1:numDims])))
      randVal <- runif(1, min = 0.0, max = 1.0)
      ## Perform rejection sampling...
      while(randVal > (testIntensity / maxIntensity)) {
        ## Create a set of test coordinates (truncating the uniform distribution when it is too
        ## far from the source coordinates to avoid inefficient rejection sampling)
        testCoords <- runif(numDims,
                            min = pmax(lowerCoords[curInd, 1:numDims], s[1:numDims] - maxWidth),
                            max = pmin(upperCoords[curInd, 1:numDims], s[1:numDims] + maxWidth))
        ## Calculate the intensity at the test coordinates
        testIntensity <- exp(-lambda * sum(abs(testCoords[1:numDims] - s[1:numDims])))
        randVal <- runif(1, min = 0.0, max = 1.0)
      }
      outCoordinates[i, 1:numDims] <- testCoords[1:numDims]
    }
    return(outCoordinates)
    returnType(double(2))
  }
)
