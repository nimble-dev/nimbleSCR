#' Local evaluation of a Bernoulli point process for activity center movement (exponential kernel)
#' 
#' Density and random generation functions of the Bernoulli point process for activity center movement. 
#'
#' The \code{dbernppLocalACmovement_exp} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' movement of activity centers between consecutive occasions in open population models.
#' The distribution assumes that the new individual activity center location (\emph{x})
#' follows an isotropic exponential normal centered on the previous activity center (\emph{s}) with rate (\emph{lambda}).
#' The local evaluation approach is implemented.
#' 
#' @name dbernppLocalACmovement_exp
#' 
#' @param x Vector of x- and y-coordinates of a single spatial point (typically AC location at time t+1).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param lowerCoords,upperCoords Matrices of lower and upper x- and y-coordinates of all habitat windows. One row for each window.
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param s Vector of x- and y-coordinates of the isotropic multivariate exponential distribution mean (AC location at time t).
#' @param lambda Rate parameter of the isotropic multivariate exponential distribution.
#' @param baseIntensities Vector of baseline habitat intensities for all habitat windows.
#' @param habitatGrid Matrix of habitat window indices. When the grid has only one row/column, artificial indices have to be provided to inflate \code{habitatGrid} 
#' in order to be able to use the distribution in \code{nimble} model code.     
#' @param habitatGridLocal Matrix of rescaled habitat grid cells indices, as returned by the \code{getLocalObjects} function (object named \code{habitatGrid}).        
#' @param resizeFactor Scalar (aggregation factor) for rescaling habitat windows as used in \code{getLocalObjects}.   
#' @param localHabWindowIndices Matrix of indices of local habitat windows around each rescaled habitat grid cell, as returned by the getLocalObjects function (object named \code{localIndices}).        
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
#' @references 
#' 
#' W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
#' A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 
#' 
#' C. Milleret, P. Dupont, C. Bonenfant, H. Broseth, O. Flagstad, C. Sutherland and R. Bischof. 2019. 
#' A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363
#' 
#' @examples 
#' 
#' # Creat habitat grid
#' habitatGrid <- matrix(c(1:(4^4)), nrow = 4, ncol=4, byrow = TRUE)
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
#' lambda <- 0.1
#' numWindows <- nrow(coordsHabitatGridCenter)
#' baseIntensities <- rep(1,numWindows)
#' numRows <- nrow(habitatGrid)
#' numCols <- ncol(habitatGrid)
#' 
#' # The log probability density of moving from (1,1) to (1.2, 0.8) 
#' dbernppLocalACmovement_exp(x = c(1.2, 0.8), lowerCoords, upperCoords, s,
#'                               lambda, baseIntensities, habitatGrid, 
#'                               HabWindowsLocal$habitatGrid, HabWindowsLocal$resizeFactor,
#'                               HabWindowsLocal$localIndices, HabWindowsLocal$numLocalIndices,
#'                               numRows, numCols, numWindows, log = TRUE)
#' 
#'
#'
NULL
#' @rdname dbernppLocalACmovement_exp
#' @export  
dbernppLocalACmovement_exp <- nimbleFunction(
  run = function(
    x                      = double(1),
    lowerCoords            = double(2),
    upperCoords            = double(2),
    s                      = double(1),
    lambda                 = double(0),
    baseIntensities        = double(1),
    habitatGrid            = double(2),
    habitatGridLocal       = double(2),
    resizeFactor           = double(0),
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
    windowInd <- habitatGrid[trunc(x[2]/resizeFactor)+1, trunc(x[1]/resizeFactor)+1]
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
    windowIntensities <- integrateIntensityLocal_exp(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                     upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                     s = s,
                                                     baseIntensities = baseIntensities[1:numWindows], 
                                                     lambda = lambda,
                                                     numLocalWindows = numLocalHabWindows[sourceAC], 
                                                     localWindows = localWindows)
    sumIntensity <- sum(windowIntensities)
    
    ## Log intensity at x
    logPointIntensity <-  log(baseIntensities[windowInd]) - lambda * sum(abs(x - s)) 
    
    ## Log probability density under the Bernoulli point process
    logProb <- logPointIntensity - log(sumIntensity)
    
    if(log) return(logProb)
    else return(exp(logProb))
    returnType(double(0))
    
  }
)

NULL
#' @rdname dbernppLocalACmovement_exp
#' @export
rbernppLocalACmovement_exp <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords            = double(2),
    upperCoords            = double(2),
    s                      = double(1),
    lambda                 = double(0),
    baseIntensities        = double(1),
    habitatGrid            = double(2),
    habitatGridLocal       = double(2),
    resizeFactor           = double(0),
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
    windowIntensities <- integrateIntensityLocal_exp(lowerCoords = lowerCoords[1:numWindows,,drop = FALSE],
                                                     upperCoords = upperCoords[1:numWindows,,drop = FALSE], 
                                                     s = s,
                                                     baseIntensities = baseIntensities[1:numWindows], 
                                                     lambda = lambda,
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
    outCoordinates <- stratRejectionSampler_exp(numPoints = 1,
                                                lowerCoords = lowerCoords1[,,drop = FALSE],
                                                upperCoords = upperCoords1[,,drop = FALSE],
                                                s = s,
                                                windowIntensities = windowIntensities[1:numWindowsLoc],
                                                lambda = lambda)
    return(outCoordinates[1,])
    returnType(double(1))
    
    
  }
)



