#' Local evaluation of a Poisson SCR observation process 
#'
#' The \code{dpoisLocal_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' Poisson observations (\emph{x}) of a single individual over a set of detectors defined by their 
#' coordinates (\emph{trapCoords}). The distribution assumes that an individual’s detection probability at any detector
#' follows a half-normal function of the distance between the  individual's activity center (\emph{s}) and the detector location.
#'
#'
#' The \code{dpoisLocal_normal} distribution incorporates three features to increase computation efficiency (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item A local evaluation of the detection probability calculation (see Milleret et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item A sparse matrix representation (\emph{x}, \emph{detIndices} and \emph{detNums}) of the observation data to reduce the size of objects to be processed.
#' \item An indicator (\emph{indicator}) to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' The \code{dpoisLocal_normal} distribution requires x- and y- detector coordinates (\emph{trapCoords}) to be scaled to the habitat grid (\emph{habitatGrid}) using the (\code{\link{scaleCoordsToHabitatGrid}} function.)
#'
#' When the aim is to simulate detection data: 
#' \enumerate{
#' \item \emph{x} should be provided using the \emph{yCombined} object as returned by \code{\link{getSparseY}}, 
#' \item arguments \emph{detIndices} and \emph{detNums} should not be provided, 
#' \item argument \emph{lengthYCombined} should be provided using the \emph{lengthYCombined} object as returned by  \code{\link{getSparseY}}.
#' }
#' 
#' 
#' 
#' @name dpoisLocal_normal
#'
#' @param x Vector of individual detection frequencies. This argument can be provided in two formats: (i) with the \emph{y} object as returned by the \code{\link{getSparseY}} function; (ii) with the \emph{yCombined} object as returned by \code{\link{getSparseY}}. 
#' Note that when the random generation functionality is used (\code{rpoisLocal_normal}), only the \emph{yCombined} format can be used. 
#' The \emph{yCombined} object combines \emph{detNums}, \emph{x}, and \emph{detIndices} (in that order).  When such consolidated representation of the detection data \emph{x} is used, \emph{detIndices} and \emph{detNums} arguments shouldn’t be specified.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param detIndices Vector of indices of traps where the detections in x were recorded, as returned by the \emph{detIndices} object from the \code{\link{getSparseY}} function. This argument should not be specified when \emph{x} is provided as the \emph{yCombined} object (returned by \code{\link{getSparseY}}) and when detection data are simulated.
#' @param detNums Number of detections recorded in \emph{x}, as returned by the \emph{detNums} object from the \code{\link{getSparseY}} function. This argument should not be specified when the \emph{yCombined} object (returned by \code{\link{getSparseY}}) is provided as \emph{x}, and when detection data are simulated.
#' @param lambda Baseline detection rate used in the half-normal detection function.
#' @param lambdaTraps Vector of baseline detection rate for each trap used in the half-normal detection function. When \emph{lambdaTraps} is used, \emph{lambda} should not be provided. 
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Individual activity center x- and y-coordinates.
#' @param trapCoords Matrix of x- and y-coordinates of all traps.
#' @param localTrapsIndices Matrix of indices of local traps around each habitat grid cell, as returned by the \code{\link{getLocalObjects}} function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the \code{\link{getLocalObjects}} function.
#' @param resizeFactor Aggregation factor used in the \code{\link{getLocalObjects}} function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid Matrix of habitat grid cells indices, as returned by the \code{\link{getLocalObjects}} function.
#' @param indicator Logical argument specifying whether the individual is available for detection.
#' @param lengthYCombined The length of the  x argument when the (\emph{yCombined}) format of the detection data is provided (as returned by the \emph{lengthYCombined} object from \code{\link{getSparseY}}). 
#' 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = lambda * exp(-d^2 / \sigma^2)}.
#'
#' @author Cyril Milleret, Soumen Dey
#'
#' @import nimble
#' @importFrom stats dpois
#' @importFrom stats rpois
#'
#' @examples
#' # I. DATA SET UP 
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
#'                                     3.5, 0.5), ncol=2,byrow = TRUE)
#' colnames(coordsHabitatGridCenter) <- c("x","y")
#' # CREATE OBSERVATION WINDOWS
#' trapCoords <- matrix(c(1.5, 1.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5), nrow = 4, byrow = TRUE)
#' colnames(trapCoords) <- c("x","y")
#' # PLOT CHECK
#' plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
#' points(trapCoords[,"y"]~trapCoords[,"x"],col="red",pch=16) 
#' 
#' # PARAMETERS
#' lambda <- 0.2
#' sigma <- 2
#' indicator <- 1 
#' # WE CONSIDER 2 INDIVIDUALS
#' y <- matrix(c(0, 1, 1, 0,
#'               0, 1, 0, 1),ncol=4,nrow=2)
#' s <- matrix(c(0.5, 1,
#'               1.6, 2.3),ncol=2,nrow=2)
#' 
#' # RESCALE COORDINATES 
#' ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
#'                                              coordsHabitatGridCenter = coordsHabitatGridCenter)
#' ScaledtrapCoords<- ScaledtrapCoords$coordsDataScaled
#' habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
#' 
#' 
#' # CREATE LOCAL OBJECTS 
#' TrapLocal <- getLocalObjects(habitatMask = habitatMask,
#'                                    coords = ScaledtrapCoords,
#'                                    dmax=2.5,
#'                                    resizeFactor = 1,
#'                                    plot.check = TRUE
#' )
#' 
#' # GET SPARSE MATRIX 
#' SparseY <- getSparseY(y)
#' 
#' # II. USING THE DENSITY FUNCTION 
#'  # WE TAKE THE FIRST INDIVIDUAL
#' i=1
#'   # OPTION 1: USING THE RANDOM GENERATION FUNCTIONNALITY 
#' dpoisLocal_normal(x=SparseY$y[i,,1],
#'                    detNums=SparseY$detNums[i],
#'                    detIndices=SparseY$detIndices[i,,1],
#'                    lambda = lambda,
#'                    sigma= sigma, 
#'                    s=s[i,1:2],
#'                    trapCoords=ScaledtrapCoords,
#'                    localTrapsIndices=TrapLocal$localIndices,
#'                    localTrapsNum=TrapLocal$numLocalIndices,
#'                    resizeFactor=TrapLocal$resizeFactor,
#'                    habitatGrid=TrapLocal$habitatGrid,
#'                    indicator=indicator)
#'                                                                 
#'   # OPTION 2: USING RANDOM GENERATION FUNCTIONNALITY 
#'   # WE DO NOT PROVIDE THE detNums AND detIndices ARGUMENTS
#' dpoisLocal_normal(x=SparseY$yCombined[i,,1],
#'                    lambda = lambda,
#'                    sigma= sigma, 
#'                    s=s[i,1:2],
#'                    trapCoords=ScaledtrapCoords,
#'                    localTrapsIndices=TrapLocal$localIndices,
#'                    localTrapsNum=TrapLocal$numLocalIndices,
#'                    resizeFactor=TrapLocal$resizeFactor,
#'                    habitatGrid=TrapLocal$habitatGrid,
#'                    indicator=indicator,
#'                    lengthYCombined = SparseY$lengthYCombined)
#' 
#' # III. USING THE RANDOM GENERATION FUNCTION 
#' rpoisLocal_normal(n=1,
#'                    lambda = lambda,
#'                    sigma= sigma, 
#'                    s=s[i,1:2],
#'                    trapCoords=ScaledtrapCoords,
#'                    localTrapsIndices=TrapLocal$localIndices,
#'                    localTrapsNum=TrapLocal$numLocalIndices,
#'                    resizeFactor=TrapLocal$resizeFactor,
#'                    habitatGrid=TrapLocal$habitatGrid,
#'                    indicator=indicator,
#'                    lengthYCombined = SparseY$lengthYCombined)
#' 
#' @export
NULL


#' @rdname dpoisLocal_normal
#' @export
dpoisLocal_normal <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  lambda = double(0, default = -999),
                  lambdaTraps = double(1),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0),
                  lengthYCombined = double(0, default = 0),
                  log = integer(0, default = 0)
  ) {
    ## Specify return type
    returnType(double(0))
    
    ## Deal with cases where detection info is combined in one vector 
    if(detNums==-999){
      detNums <- x[1]
      nMaxDetectors <- (lengthYCombined-1)/2
      detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
      x1 <- x[2:(nMaxDetectors+1)]
    }else{
      x1 <- x
      detIndices1 <- detIndices
    }
    
    ## Shortcut if the current individual is not available for detection
    if(indicator == 0){
      if(detNums == 0){
        if(log == 0) return(1.0)
        else return(0.0)
      } else {
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
    
    ## Retrieve the index of the habitat cell where the current AC is
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## Retrieve the indices of the local traps surrounding the selected habita grid cell
    theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
    
    ## CHECK IF DETECTIONS ARE WITHIN THE LIST OF  LOCAL TRAPS
    if(detNums > 0){
      for(r in 1:detNums){
        if(sum(detIndices1[r] == theseLocalTraps) == 0){
          if(log == 0) return(0.0)
          else return(-Inf)
        }
      }
    }
    
    ## Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logProb <- 0.0 
    detIndices1 <- c(detIndices1,0)
    count <- 1 
    
    
    if(lambda==-999){# when lambda is provide through lambdaTraps
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- lambdaTraps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <- logProb + dpois(x1[count], p, log = TRUE)
          count <- count + 1
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- lambdaTraps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <- logProb + dpois(0, p, log = TRUE)
          
        }
      }
    }else{# when lambda is provide through lambda
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- lambda * exp(alpha * d2)
          logProb <- logProb + dpois(x1[count], p, log = TRUE)
          count <- count + 1
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p <- lambda * exp(alpha * d2)
          logProb <- logProb + dpois(0, p, log = TRUE)
          
        }
      }
    }
    
  
    
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dpoisLocal_normal
#' @export
rpoisLocal_normal <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  lambda = double(0, default = -999),
                  lambdaTraps = double(1),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0),
                  lengthYCombined = double(0, default = 0)
  ) {
    ## Specify return type
    returnType(double(1))
    if(detNums >= 0) stop("Random generation for the rpoisLocal_normal distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    
    #========================================================
    # RETURN TYPE DECLARATION
    if(n!=1){print("rpoisLocal_normal only allows n = 1; using n = 1")}
    # returnType(double(3))
    # len <- 2*MAX + 1
    ## GET NECESSARY INFO
    alpha <- -1.0 / (2.0 * sigma * sigma)
    # n.detectors <- dim(detector.xy)[1]
    # nMAxDetections <- length(detIndices)
    nMAxDetections <- (lengthYCombined-1)/2
    ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
    #if(indicator == 0){return(rep(0.0, 2*nMAxDetections + 1))}
    if(indicator == 0){return(rep(0.0, lengthYCombined))}
    
    ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
    
    ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
    detectOut <- rep(0, localTrapsNum[sID])
    ys <- rep(-1, nMAxDetections)
    dets <- rep(-1, nMAxDetections)
    count <- 1
    
    ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
    if(lambda==-999){# when lambda is provide through lambdaTraps
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p <- lambdaTraps[theseLocalTraps[r]] * exp(alpha * d2)
        # Draw the observation at detector j from a binomial distribution with probability p
        detectOut[r] <- rpois(1, p)
        if(detectOut[r] >0){
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                        You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          ys[count] <- detectOut[r]
          dets[count] <- theseLocalTraps[r]
          count <- count + 1
        }#if
      }#r
    }else{# when lambda is provide through lambda
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p <- lambda * exp(alpha * d2)
        # Draw the observation at detector j from a binomial distribution with probability p
        detectOut[r] <- rpois(1, p)
        if(detectOut[r] >0){
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                        You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          ys[count] <- detectOut[r]
          dets[count] <- theseLocalTraps[r]
          count <- count + 1
        }#if
      }#r
    }  
      
    count <- count - 1
    
    
    # out <- rep(-1, 2*nMAxDetections + 1)
    out <- rep(-1, lengthYCombined)
    
    out[1] <- count
    if(count >= 1){
      out[2:(count+1)] <- ys[1:count]
      out[(nMAxDetections+2):(nMAxDetections+count+1)] <- dets[1:count]
    }
    ## OUTPUT
    return(out)
  })
