#' Poisson point process detection model with local evaluation
#' 
#' Density and random generation functions of the Poisson point process for detection. The local evaluation technique is implemented. 
#' An isotropic multivariate normal distribution is used as the decay kernel. 
#' 
#' @name dpoisppLocalDetection_normal
#' 
#' @param x Matrix of x- and y-coordinates and the corresponding id the detection window of a set of spatial points (detection locations). One row corresponds to one point.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' @param s Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (the AC location).
#' @param sd Standard deviation of the isotropic multivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param habitatGridLocal Matrix of rescaled habitat grid cells indices, as returned by the \code{getLocalObjects} function (object named \code{habitatGrid}).      
#' @param resizeFactor Scalar (aggregation factor) for rescaling habitat windows as used in \code{getLocalObjects}.   
#' @param localObsWindowIndices Matrix of indices of local observation windows around each rescaled habitat grid cell, as returned by the getLocalObjects function (object named \code{localIndices}).        
#' @param numLocalObsWindows Vector of numbers of local observation windows around all habitat grid cells, as returned by the getLocalObjects function (object named \code{numLocalIndices}). 
#' The ith number gives the number of local (original) observation windows for the ith (rescaled) habitat window.
#' @param numMaxPoints Maximum number of points. This value (non-negative integer) is only used when simulating detections to constrain the maximum number of detections. 
#' @param numWindows Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored. 
#' @param indicator Binary variable (0 or 1) used for data augmentation. \code{indicator = 0} means the individual does not exist 
#' and thus the probability of no detection is 1. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return The (log) probability density of the observation matrix \code{x}.
#' 
#' @author Wei Zhang and Cyril Milleret
#' 
#' @import nimble
#' @importFrom stats dnorm
#' 
#' @references
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' C. Milleret, P. Dupont, C. Bonenfant, H. Brøseth, Ø. Flagstad, C. Sutherland and R. Bischof. 2019. 
#' A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363
#'
#'#' @examples 
#' # Create habitat grid
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
#'                                     3.5, 0.5), ncol = 2, byrow = TRUE)
#' colnames(coordsHabitatGridCenter) <- c("x","y")
#' # Create observation windows
#' lowerCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(2, 2, 3, 2, 2, 3, 3, 3), nrow = 4, byrow = TRUE)  
#' colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
#' # Plot check
#' plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
#' points(lowerCoords[,"y"]~lowerCoords[,"x"],col="red",pch=16) 
#' points(upperCoords[,"y"]~upperCoords[,"x"],col="red",pch=16) 
#' 
#' s <- c(1, 1)
#' sd <- 0.1
#' baseIntensities <- c(1:4)
#' windowIndex <- 4
#' numPoints <- 1
#' numWindows <- 4
#' indicator <- 1
#' 
#' # Rescale coordinates
#' ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData =  lowerCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
#' ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
#' ScaledUpperCoords[,2] <- ScaledUpperCoords[,2] + 1.5
#' ScaledLowerCoords[,2] <- ScaledLowerCoords[,2] - 1.5
#' habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
#' # Create local objects 
#' ObsWindowsLocal <- getLocalObjects(habitatMask = habitatMask,
#'                                    coords = ScaledLowerCoords,
#'                                    dmax=3,
#'                                    resizeFactor = 1,
#'                                    plot.check = TRUE
#' )
#' 
#' # Detection locations
#' x <- matrix(c(1.5, 2, 1.1, 1.5, 1.4, 0.7, 2, 1.3, 1, 1.5), nrow = 5, byrow = TRUE)
#' 
#' # get the window indeces on the third dimension of x
#' windowIndexes <- 0
#' for(i in 1:nrow(x)){
#'   windowIndexes[i] <- getWindowIndex(curCoords = x[i,],
#'                                      lowerCoords = ScaledLowerCoords,
#'                                      upperCoords =ScaledUpperCoords)
#' }
#' x <- cbind(x, windowIndexes)
#' # get the total number of detections on x[1,1]
#' x <- rbind(c(length(windowIndexes),0,0) ,x )
#' dpoisppLocalDetection_normal(x, ScaledLowerCoords, ScaledUpperCoords, s, sd, baseIntensities,  
#'                              ObsWindowsLocal$habitatGrid, ObsWindowsLocal$resizeFactor,
#'                              ObsWindowsLocal$localIndices,ObsWindowsLocal$numLocalIndices,
#'                              numMaxPoints = dim(x)[1], numWindows, indicator, log = TRUE)
NULL
#' @rdname dpoisppLocalDetection_normal
#' @export
dpoisppLocalDetection_normal <- nimbleFunction(
  run = function(
    x                     = double(2),
    lowerCoords           = double(2),
    upperCoords           = double(2),
    s                     = double(1),
    sd                    = double(0),
    baseIntensities       = double(1),
    habitatGridLocal      = double(2),
    resizeFactor          = double(0),
    localObsWindowIndices = double(2),
    numLocalObsWindows    = double(1),
    numMaxPoints          = integer(0),
    numWindows            = integer(0),
    indicator             = integer(0),
    log                   = integer(0, default = 0)
  ) {
    
    ## If the individual does not exists
    if(indicator == 0) {
      if(x[1,1] == 0) {
        if(log) return(0.0)
        else return(1.0)
      }
      else{
        if(log) return(-Inf)
        else return(0.0)
      }
    }
    ## Local evaluation:
    ## Get the habitat window index where the AC is located
    habWindowInd <- habitatGridLocal[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    ## Get the indices of detection windows that are close to the AC
    localWindows <- localObsWindowIndices[habWindowInd, 1:numLocalObsWindows[habWindowInd]]
    
    numPoints1 <- x[1,1]+1
    ## Check if any point is out of the local area of the AC
    if(x[1,1] > 0) {
      for(i in 2:numPoints1) {
        #if(sum(localWindows == windowIndices[i]) == 0){
        if(sum(localWindows == x[i,3]) == 0){
          if(log) return(-Inf)
          else return(0.0)
        }
      }
    }
    ## Integrate the detection intensity function over all detection windows
    windowIntensities <- integrateIntensityLocal_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                        upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                        s = s,
                                                        baseIntensities = baseIntensities[1:numWindows], 
                                                        sd = sd,
                                                        numLocalWindows = numLocalObsWindows[habWindowInd],
                                                        localWindows = localWindows)
    sumIntensity <- sum(windowIntensities)
    ## Calculate the probability
    if(x[1,1] == 0) {
      logProb <- -sumIntensity
    } else{
      #numDims <- 2 ## We consider 2D models for now
      #constant <- (numDims / 2.0) * log(2.0 * pi) # 1.837877
      logPointIntensity <- rep(1.837877, x[1,1])
      for(i in 1:x[1,1]) {
        ## Log intensity at the ith point: see Eqn 24 of Zhang et al.(2020, DOI:10.1101/2020.10.06.325035)
        pointBaseIntensity <- baseIntensities[x[i+1,3]]
        logPointIntensity[i] <- logPointIntensity[i] + log(pointBaseIntensity) + sum(dnorm((x[i+1,1:2] - s) / sd, log = 1))
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
#' @rdname dpoisppLocalDetection_normal
#' @export
rpoisppLocalDetection_normal <- nimbleFunction(
  run = function(
    n                     = integer(0),
    lowerCoords           = double(2),
    upperCoords           = double(2),
    s                     = double(1),
    sd                    = double(0),
    baseIntensities       = double(1),
    habitatGridLocal      = double(2),
    resizeFactor          = double(0),
    localObsWindowIndices = double(2),
    numLocalObsWindows    = double(1),
    numMaxPoints          = integer(0),
    numWindows            = integer(0),
    indicator             = integer(0)
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
      ## Local evaluation:
      ## Get the habitat window index where the AC is located
      habWindowInd <- habitatGridLocal[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
      ## Get the indices of detection windows that are close to the AC
      numWindowsLoc <- numLocalObsWindows[habWindowInd]
      localWindows <- localObsWindowIndices[habWindowInd, 1:numWindowsLoc]      
      
      ## Integrate the detection intensity function over all detection windows
      windowIntensities <- integrateIntensityLocal_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                          upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                          s = s,
                                                          baseIntensities = baseIntensities[1:numWindows], 
                                                          sd = sd,
                                                          numLocalWindows = numWindowsLoc,
                                                          localWindows = localWindows)
      sumIntensity <- sum(windowIntensities)
      
      if(sumIntensity > 0.0) {
        
        outCoordinates[1,1] <-  rpois(1, sumIntensity)
        numPoints1 <- outCoordinates[1,1]+1
        #numPoints1 <- rpois(1, sumIntensity)
        if(outCoordinates[1,1] >numMaxPoints){stop("There are more simulated individual detections than what can be stored within x.\n You may need to increase the size of the 'x' object and make sure that 'numPoints'is equal to dim(x)[2]" )}
        if(outCoordinates[1,1]  > 0) {
          
          ## DO THE SUBSETTING OF THE LOWER AND UPPER COORDS HERE. 
          lowerCoords1 <- nimMatrix(nrow = numWindowsLoc, ncol = 2)
          upperCoords1 <- nimMatrix(nrow = numWindowsLoc, ncol = 2)
          
          for(i in 1:numWindowsLoc){
            lowerCoords1[i,1:2] <- lowerCoords[localWindows[i],,drop = FALSE]
            upperCoords1[i,1:2] <- upperCoords[localWindows[i],,drop = FALSE]
          }
          
          
          
          outCoordinates[2:numPoints1 , 1:2 ] <- stratRejectionSampler_normal(numPoints = outCoordinates[1,1],
                                                                              lowerCoords = lowerCoords1[,,drop = FALSE],
                                                                              upperCoords = upperCoords1[,,drop = FALSE],
                                                                              s = s,
                                                                              windowIntensities = windowIntensities,
                                                                              sd = sd)
          # GET WINDOW INDEX
          for(i in 2:numPoints1 ){
            outCoordinates[i, 3 ] <- getWindowIndex(curCoords = outCoordinates[i, 1:2],
                                                    lowerCoords = lowerCoords[1:numWindows,,drop = FALSE] ,
                                                    upperCoords = upperCoords[1:numWindows,,drop = FALSE] )
          }
        }
      }
    }
    return(outCoordinates)
    returnType(double(2))
  }
)



