#' Local evaluation of a Bernoulli point process for activity center movement (normal kernel)
#' 
#' Density and random generation functions of the Bernoulli point process for activity center movement between occasions based on a bivariate normal distribution and local evaluation.
#' 
#' The \code{dbernppLocalACmovement_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' movement of activity centers between consecutive occasions in open population models.
#' The distribution assumes that the new individual activity center location (\emph{x})
#' follows an isotropic multivariate normal centered on the previous activity center (\emph{s}) with standard deviation (\emph{sd}).
#' The local evaluation technique is implemented.
#' 
#' 
#' @name dbernppLocalACmovement_normal
#' 
#' @param x Vector of x- and y-coordinates of a single spatial point (typically AC location at time t+1) scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}). 
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows. One row for each window.
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param s Vector of x- and y-coordinates of the isotropic bivariate normal distribution mean (i.e. the AC location).
#' @param sd Standard deviation of the isotropic bivariate normal distribution.
#' @param baseIntensities Vector of baseline habitat intensities for all habitat windows.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in  \code{lowerCoords} and \code{upperCoords}. 
#'  When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.     
#' @param habitatGridLocal Matrix of rescaled habitat grid cells indices, from localIndices returned by the  \code{getLocalObjects} function (        
#' @param resizeFactor Aggregation factor used in the  \code{getLocalObjects} function to reduce the number of habitat grid cells.   
#' @param localHabWindowIndices Matrix of indices of local habitat windows around each local habitat grid cell (\code{habitatGridLocal}), from localIndices returned by the \code{getLocalObjects} function.        
#' @param numLocalHabWindows Vector of numbers of local habitat windows around all habitat grid cells, as returned by the getLocalObjects function (object named \code{numLocalIndices}). 
#' The ith number gives the number of local (original) habitat windows for the ith (rescaled) habitat window.
#' @param numGridRows,numGridCols Numbers of rows and columns of the \code{habitatGrid}.
#' @param numWindows Number of habitat windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
#' so that extra rows beyond \code{numWindows} are ignored. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' 
#' @return The (log) probability density of the observation vector \code{x}.
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
#' @examples 
#' 
#' # Creat habitat grid
#' habitatGrid <- matrix(c(1:(4^2)), nrow = 4, ncol=4, byrow = TRUE)
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
#' # Create habitat windows
#' lowerCoords <- coordsHabitatGridCenter-0.5
#' upperCoords <- coordsHabitatGridCenter+0.5
#' colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
#' # Plot check
#' plot(lowerCoords[,"y"]~lowerCoords[,"x"],pch=16, xlim=c(0,4), ylim=c(0,4),col="red") 
#' points(upperCoords[,"y"]~upperCoords[,"x"],col="red",pch=16) 
#' points(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
#' 
#' # Rescale coordinates 
#' ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData =  lowerCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
#'                                               coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledUpperCoords$coordsDataScaled[,2] <- ScaledUpperCoords$coordsDataScaled[,2] + 1
#' ScaledLowerCoords$coordsDataScaled[,2] <- ScaledLowerCoords$coordsDataScaled[,2] - 1
#' habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
#' # Create local objects 
#' HabWindowsLocal <- getLocalObjects(habitatMask = habitatMask,
#'                                    coords = coordsHabitatGridCenter,
#'                                    dmax=4,
#'                                    resizeFactor = 1,
#'                                    plot.check = TRUE
#' )
#' 
#' s <- c(1, 1) # Currrent activity center location
#' sd <- 0.1
#' numWindows <- nrow(coordsHabitatGridCenter)
#' baseIntensities <- rep(1,numWindows)
#' numRows <- nrow(habitatGrid)
#' numCols <- ncol(habitatGrid)
#' 
#' # The log probability density of moving from (1,1) to (1.2, 0.8) 
#' dbernppLocalACmovement_normal(x = c(1.2, 0.8), lowerCoords, upperCoords, s,
#'                               sd, baseIntensities, habitatGrid, 
#'                               HabWindowsLocal$habitatGrid, HabWindowsLocal$resizeFactor,
#'                               HabWindowsLocal$localIndices, HabWindowsLocal$numLocalIndices,
#'                               numRows, numCols, numWindows, log = TRUE)
#' 
#'
NULL
#' @rdname dbernppLocalACmovement_normal
#' @export
dbernppLocalACmovement_normal <- nimbleFunction(
  run = function(
    x                      = double(1),
    lowerCoords            = double(2),
    upperCoords            = double(2),
    s                      = double(1),
    sd                     = double(0),
    baseIntensities        = double(1),
    habitatGrid            = double(2),
    habitatGridLocal       = double(2),
    resizeFactor           = double(0, default = 1),
    localHabWindowIndices  = double(2),
    numLocalHabWindows     = double(1),
    numGridRows            = integer(0),
    numGridCols            = integer(0),
    numWindows             = integer(0),
    log                    = integer(0, default = 0)
  ) {
    ## Check if the point x falls within the habitat
    if(min(x) < 0 | x[2] >= numGridRows | x[1] >= numGridCols) {
      if(log) return(-Inf) 
      else return(0.0)
    }
    ## Index of the window where x falls
    windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1]
    ## windowInd == 0 means this window is not defined as habitat
    if(windowInd == 0) {
        if(log) return(-Inf)
        else return(0.0)
    }
    ## Find in which habitat window (from the rescaled habitat grid) the s (source AC location) falls in
    sourceAC <- habitatGridLocal[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    ## Get local windows ids within a close distance from the source AC  
    localWindows <- localHabWindowIndices[sourceAC, 1:numLocalHabWindows[sourceAC]]
    ## Check if the point x is out of the local area of the AC
    if(sum(localWindows == windowInd) == 0){
      if(log) return(-Inf)
      else return(0.0)
    }
    ## Integrate the intensity function over all habitat windows
    windowIntensities <- integrateIntensityLocal_normal(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                        upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                        s = s,
                                                        baseIntensities = baseIntensities[1:numWindows], 
                                                        sd = sd,
                                                        numLocalWindows = numLocalHabWindows[sourceAC], 
                                                        localWindows = localWindows)
    sumIntensity <- sum(windowIntensities)
    ## Log intensity at x
    # numDims <- 2
    # logPointIntensity <- (numDims / 2.0) * log(2.0 * pi) + log(baseIntensities[windowInd]) + sum(dnorm((x - s) / sd, log = 1))
    logPointIntensity <- 1.837877 + log(baseIntensities[windowInd]) + sum(dnorm((x - s) / sd, log = 1))
    
    ## Log probability density under the Bernoulli point process
    logProb <- logPointIntensity - log(sumIntensity)
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
  }
)

NULL
#' @rdname dbernppLocalACmovement_normal
#' @export 
#' 
rbernppLocalACmovement_normal <- nimbleFunction(
  run = function(
    n                      = integer(0),
    lowerCoords            = double(2),
    upperCoords            = double(2),
    s                      = double(1),
    sd                     = double(0),
    baseIntensities        = double(1),
    habitatGrid            = double(2),
    habitatGridLocal       = double(2),
    resizeFactor           = double(0, default = 1),
    localHabWindowIndices  = double(2),
    numLocalHabWindows     = double(1),
    numGridRows            = integer(0),
    numGridCols            = integer(0),
    numWindows             = integer(0)
  ) {
    ## Ensure that only one sample is requested
    if(n <= 0) {
      stop("The number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernppACmovement only allows n = 1; using n = 1")
    }
    ## Find in which habitat window (from the rescaled habitat grid) the s (source AC location) falls in
    sourceAC <- habitatGridLocal[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    ## Get local windows ids within a close distance from the source AC  
    numWindowsLoc <- numLocalHabWindows[sourceAC] 
    localWindows <- localHabWindowIndices[sourceAC, 1:numWindowsLoc]
    
    ## Integrate the intensity function over all habitat windows
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
    
    ## Call the statified rejection sampler
    outCoordinates <- stratRejectionSampler_normal(numPoints = 1,
                                                   lowerCoords = lowerCoords1[,,drop = FALSE],
                                                   upperCoords = upperCoords1[,,drop = FALSE],
                                                   s = s,
                                                   windowIntensities = windowIntensities[1:numWindowsLoc],
                                                   sd = sd)
    return(outCoordinates[1,])
    returnType(double(1))
    
    
  }
  
)





