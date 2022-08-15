#' Poisson point process detection model
#' 
#' Density and random generation functions of the Poisson point process for detection. 
#
#' The \code{dpoisppDetection_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' Poisson observations (\emph{x}) of a single individual in continuous space over a set of detection windows defined by their upper and lower
#' coordinates (\emph{lowerCoords,upperCoords}). The distribution assumes that an individual’s detection intensity 
#' follows an isotropic multivariate normal centered on the individual's activity center (\emph{s}) with standard deviation (\emph{sd}).
#' 
#' 
#' @name dpoisppDetection_normal
#' 
#' @param x Array containing the total number of detections (x[1,1]), the x- and y-coordinates (x[2:(x[1,1]+1),1:2]), 
#' and the corresponding detection window indices (x[2:(x[1,1]+1),3]) for a set of spatial points (detection locations). 
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (the AC location).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param numMaxPoints Maximum number of points. This value (non-negative integer) is only used when simulating detections to constrain the maximum number of detections. 
#' @param numWindows Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored.  
#' @param indicator Binary variable (0 or 1) used for data augmentation. \code{indicator = 0} means the individual does not exist 
#' and thus the probability of no detection is 1.  
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return 
#' \code{dpoisppDetection_normal} gives the (log) probability density of the observation matrix \code{x}.
#' \code{rpoisppDetection_normal} gives coordinates of a set of randomly generated spatial points.  
#' 
#' @author Wei Zhang 
#'
#' @import nimble
#' @importFrom stats dnorm
#'
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' @examples 
#' coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
#'                                     1.5, 3.5,
#'                                     2.5, 3.5,
#'                                     3.5, 3.5,
#'                                     0.5, 2.5,
#'                                     1.5, 2.5,
#'                                     2.5, 2.5,
#'                                     3.5, 2.5,
#'                                     0.5, 1.5,
#'                                     1.5, 1.5,
#'                                     2.5, 1.5,
#'                                     3.5, 1.5,
#'                                     0.5, 0.5,
#'                                     1.5, 0.5,
#'                                     2.5, 0.5,
#'                                     3.5, 0.5), ncol = 2,byrow = TRUE)
#' colnames(coordsHabitatGridCenter) <- c("x","y")
#' # Create observation windows
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
#' 
#' # Rescale coordinates
#' ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData =  lowerCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords$coordsDataScaled[,2] <- ScaledUpperCoords$coordsDataScaled[,2] + 1.5
#' ScaledLowerCoords$coordsDataScaled[,2] <- ScaledLowerCoords$coordsDataScaled[,2] - 1.5
#' 
#' 
#' 
#' # Detection locations
#' x <- matrix(c(1.5, 2, 1.1, 1.5,0.6, 2.1, 0.5, 2, 1, 1.5), nrow = 5, byrow = TRUE)
#' 
#' # get the window indeces on the third dimension of x
#' 
#' windowIndexes <- 0
#' for(i in 1:nrow(x)){
#'   windowIndexes[i] <- getWindowIndex(curCoords = x[i,],
#'                                      lowerCoords = ScaledLowerCoords$coordsDataScaled,
#'                                      upperCoords = ScaledUpperCoords$coordsDataScaled)
#' }
#' x <- cbind(x, windowIndexes)
#' # get the total number of detections on x[1,1]
#' x <- rbind(c(length(windowIndexes),0,0) ,x )
#' 
#' 
#' s <- c(1, 1)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' windowIndices <- c(1, 2, 2, 3, 4)
#' numPoints <- 5
#' numWindows <- 4
#' indicator <- 1
#' dpoisppDetection_normal(x, lowerCoords, upperCoords, s, sd, baseIntensities,
#'                         numMaxPoints = dim(x)[1] , numWindows, indicator, log = TRUE)

NULL
#' @rdname dpoisppDetection_normal
#' @export
dpoisppDetection_normal <- nimbleFunction(
  run = function(
    x               = double(2),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    numMaxPoints       = integer(0),
    numWindows    = double(0),
    indicator       = integer(0),
    log             = integer(0, default = 0)
  ) {
    ## If the individual does not exists 
    if(indicator == 0) {
      #if(numPoints == 0){
      if(x[1,1] == 0){
        if(log) return(0.0)
        else return(1.0)
      } else{
        if(log) return(-Inf)
        else return(0.0) 
      }
    }
    ## Integrate the detection intensity function over all detection windows
    windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE], 
                                                   upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                   s = s,
                                                   baseIntensities = baseIntensities[1:numWindows], 
                                                   sd = sd,
                                                   numWindows = numWindows)
    sumIntensity <- sum(windowIntensities) ## Expected total number of detections
    ## Calculate the probability
    if(x[1,1] == 0) {
      logProb <- -sumIntensity
    } else{
      ## 1.837877 <- log(2.0 * pi) 
      logPointIntensity <- rep(1.837877, x[1,1])
      for(i in 1:x[1,1]) {
        ## Log intensity at the ith point: see Eqn 24 of Zhang et al.(2020, DOI:10.1101/2020.10.06.325035)
        pointBaseIntensity <- baseIntensities[x[i+1,3]]
        logPointIntensity[i] <- logPointIntensity[i] + log(pointBaseIntensity) + sum(dnorm((x[i+1, 1:2] - s) / sd, log = 1))
      }
      ## Log probability density
      logProb <- sum(logPointIntensity) - sumIntensity
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dpoisppDetection_normal
#' @export
rpoisppDetection_normal <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    s               = double(1),
    sd              = double(0),
    baseIntensities = double(1),
    numMaxPoints    = double(0),
    numWindows      = integer(0),
    indicator       = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rpoisppDetection_normal only allows n = 1; using n = 1")
    }
    numDims <- 3 
    ## If the individual does not exist in the population, we return an empty matrix directly
    if(indicator == 0) return(matrix(nrow = numMaxPoints, ncol = numDims))
    
    ## Initialize an output matrix
    ## This empty matrix will be returned if numWindows==0 or numPoints==0
    outCoordinates <- matrix(0, nrow = numMaxPoints, ncol = numDims)
    
    ## Generate random points below:
    ## Current implementation uses a stratified rejection sampler. 
    ## This can become inefficient if the observation window becomes very large compared to the 
    ## standard deviation of the multivariate normal distribution. 
    if((numWindows > 0) & (numMaxPoints != 0)) {
      ## Integrate the detection intensity function over all detection windows
      windowIntensities <- integrateIntensity_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                     upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                     s = s,
                                                     baseIntensities = baseIntensities[1:numWindows],
                                                     sd = sd,
                                                     numWindows = numWindows)
      
      sumIntensity <- sum(windowIntensities)
      if(sumIntensity > 0.0) {
        
        outCoordinates[1,1] <-  rpois(1, sumIntensity)
        numPoints1 <- outCoordinates[1,1]+1
        if(outCoordinates[1,1] >numMaxPoints){stop("There are more simulated individual detections than what can be stored within x.\n You may need to increase the size of the 'x' object and make sure that 'numMaxPoints'is equal to dim(x)[2]" )}
        if(outCoordinates[1,1] > 0) {
          outCoordinates[2:numPoints1,1:2] <- stratRejectionSampler_normal(numPoints = outCoordinates[1,1],
                                                                           lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                                           upperCoords = upperCoords[1:numWindows,,drop = FALSE],
                                                                           s = s,
                                                                           windowIntensities = windowIntensities[1:numWindows],
                                                                           sd = sd)
          ## Get detection window index
          for(i in 2:numPoints1){
            outCoordinates[i,3] <- getWindowIndex(curCoords = outCoordinates[i, 1:2],
                                                  lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                  upperCoords = upperCoords[1:numWindows,,drop = FALSE])
          }
        }
      }
    }
    return(outCoordinates)
    returnType(double(2))
  }
)



