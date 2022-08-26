#' Local evaluation for a Bernoulli point process detection model
#' 
#' Density and random generation functions of the Bernoulli point process for detection based on a bivariate normal distribution. 
#' 
#' The \code{dbernppDetection_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' Bernoulli observations (\emph{x}) of a single individual in continuous space over a set of detection windows defined by their upper and lower
#' coordinates (\emph{lowerCoords,upperCoords}). The distribution assumes that an individual’s detection probability 
#' follows an isotropic multivariate normal centered on the individual's activity center (\emph{s}) with standard deviation (\emph{sd}).
#' The local evaluation approach is implemented.
#' 
#' 
#' @name dbernppLocalDetection_normal
#' 
#' @param x Vector with three elements representing the x- and y-coordinates (x[1:2]), and the corresponding id the detection window (x[3]) of a single spatial point (detection location) scaled to the habitat (see \code{scaleCoordsToHabitatGrid}).
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param s Vector of x- and y-coordinates of the isotropic bivariate normal distribution mean (i.e. the AC location).
#' @param sd Standard deviation of the bivariate normal distribution.
#' @param baseIntensities Vector of baseline detection intensities for all detection windows.
#' @param habitatGridLocal Matrix of rescaled habitat grid cells indices, as returned by the \code{getLocalObjects} function (object named \code{habitatGrid}).        
#' @param resizeFactor Aggregation factor used in the  \code{getLocalObjects} function to reduce the number of habitat grid cells.   
#' @param localObsWindowIndices Matrix of indices of local observation windows around each local habitat grid cell (habitatGridLocal), from localIndices returned by the \code{getLocalObjects} function.
#' @param numLocalObsWindows Vector of numbers of local observation windows around all habitat grid cells, as returned by the getLocalObjects function (object named \code{numLocalIndices}). 
#' The ith number gives the number of local (original) observation windows for the ith (rescaled) habitat window.
#' @param numWindows Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored. 
#' @param indicator Binary argument specifying whether the individual is available for detection (indicator = 1) or not (indicator = 0).
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return The (log) probability density of the observation vector \code{x}.
#' 
#' @author Wei Zhang and Cyril Milleret
#' 
#' 
#' @references 
#'
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' C. Milleret, P. Dupont, C. Bonenfant, H. Brøseth, Ø. Flagstad, C. Sutherland and R. Bischof. 2019. 
#' A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363
#'
#' @examples 
#' 
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
#'                                     3.5, 0.5), ncol = 2,byrow = TRUE)
#' colnames(coordsHabitatGridCenter) <- c("x","y")
#' # Create observation windows
#' lowerCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(2, 2, 3, 2, 2, 3, 3, 3), nrow = 4, byrow = TRUE)  
#' colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
#' # Plot check
#' plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
#' points(lowerCoords[,"y"]~lowerCoords[,"x"],col="red",pch=16) 
#' points(upperCoords[,"y"]~upperCoords[,"x"],col="red",pch=16) 
#' #'
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
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords$coordsDataScaled[,2] <- ScaledUpperCoords$coordsDataScaled[,2] + 1.5
#' ScaledLowerCoords$coordsDataScaled[,2] <- ScaledLowerCoords$coordsDataScaled[,2] - 1.5
#' habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
#' # Create local objects 
#' ObsWindowsLocal <- getLocalObjects(habitatMask = habitatMask,
#'                                    coords = ScaledLowerCoords$coordsDataScaled,
#'                                    dmax=3,
#'                                    resizeFactor = 1,
#'                                    plot.check = TRUE
#' )
#' x <- c(1.1, 1.2)
#' windowIndex <- getWindowIndex(curCoords = x,
#'                               lowerCoords = ScaledLowerCoords$coordsDataScaled,
#'                               upperCoords =ScaledUpperCoords$coordsDataScaled)
#' x <- c(x, windowIndex)
#' dbernppLocalDetection_normal(x, ScaledLowerCoords$coordsDataScaled,
#'                              ScaledUpperCoords$coordsDataScaled,
#'                              s, sd, baseIntensities,  
#'                              ObsWindowsLocal$habitatGrid, ObsWindowsLocal$resizeFactor,
#'                              ObsWindowsLocal$localIndices,ObsWindowsLocal$numLocalIndices,
#'                              numWindows, indicator, log = TRUE)
NULL
#' @rdname dbernppLocalDetection_normal
#' @export
dbernppLocalDetection_normal <- nimbleFunction(
  run = function(
    x                     = double(1),
    lowerCoords           = double(2),
    upperCoords           = double(2),
    s                     = double(1),
    sd                    = double(0),
    baseIntensities       = double(1),
    habitatGridLocal      = double(2),
    resizeFactor          = double(0),
    localObsWindowIndices = double(2),
    numLocalObsWindows    = double(1),
    numWindows            = integer(0),
    indicator             = integer(0),
    log                   = integer(0, default = 0)
  ) {
    ## If the individual does not exists
    if(indicator == 0) {
      if(x[3] == 0){
        if(log) return(0.0) else return(1.0) 
      } else {
        if(log) return(-Inf) else return(0.0) 
      }
    }
    ## Make sure individuals with indicator==1 have one detection (necessary for dead recovery)
    if(indicator == 1) {
      if(x[3] == 0) {
        if(log) return(-Inf) else return(0.0)
      }
    }
    ## Get the habitat window index where the AC is located
    habWindowInd <- habitatGridLocal[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    ## Get the indices of detection windows that are close to the AC
    localWindows <- localObsWindowIndices[habWindowInd, 1:numLocalObsWindows[habWindowInd]]
    
    ## Check if the point is out of the local area of the AC
    ## windowIndex is the index of the detection window where the point x falls
    if(sum(localWindows == x[3]) == 0){
      if(log) return(-Inf)
      else return(0.0)
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
    ## Make sure sumIntensity is positive 
    if(sumIntensity <= 0.0) {
      if(log) return(-Inf)
      else return(0.0)
    }
    # numDims <- 2 
    ## Log intensity at x
    #logPointIntensity <- (numDims / 2.0) * log(2.0 * pi) + log(baseIntensities[windowIndex]) + sum(dnorm((x - s) / sd, log = 1))
    logPointIntensity <- 1.837877 + log(baseIntensities[x[3]]) + sum(dnorm((x[1:2] - s) / sd, log = 1))
    logProb <- logPointIntensity - log(sumIntensity)
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dbernppLocalDetection_normal
#' @export
rbernppLocalDetection_normal <- nimbleFunction(
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
    numWindows            = integer(0),
    indicator             = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernppDetection_normal only allows n = 1; using n = 1")
    }
    
    if(indicator==0){return(c(0,0,0))}else{
      
      
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
      
      
      ## DO THE SUBSETTING OF THE LOWER AND UPPER COORDS HERE. 
      lowerCoords1 <- nimMatrix(nrow = numWindowsLoc, ncol = 2)
      upperCoords1 <- nimMatrix(nrow = numWindowsLoc, ncol = 2)
      
      for(i in 1:numWindowsLoc){
        lowerCoords1[i,1:2] <- lowerCoords[localWindows[i],,drop = FALSE]
        upperCoords1[i,1:2] <- upperCoords[localWindows[i],,drop = FALSE]
      }
      
      outCoordinates <- stratRejectionSampler_normal(numPoints = 1,
                                                     lowerCoords = lowerCoords1[,,drop = FALSE],
                                                     upperCoords = upperCoords1[,,drop = FALSE],
                                                     s = s,
                                                     windowIntensities = windowIntensities,
                                                     sd = sd)
      windowIndex <- getWindowIndex(curCoords = outCoordinates[1, 1:2],
                                    lowerCoords = lowerCoords[1:numWindows,,drop = FALSE] ,
                                    upperCoords = upperCoords[1:numWindows,,drop = FALSE] )
      
      return(c(outCoordinates[1,],windowIndex))    }
    returnType(double(1))
  }
)


