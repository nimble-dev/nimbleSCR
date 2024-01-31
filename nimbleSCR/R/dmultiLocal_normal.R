#' Local evaluation of a multinomial SCR detection process 
#'
#' The \code{dmultiLocal_normal} distribution is a NIMBLE custom distribution which can be used to model 
#' and simulate multinomial observations (x) of a single individual over a set of traps defined by their coordinates \emph{trapCoords}
#' the distribution assumes that an individual’s detection probability at any trap follows a half-normal function of the distance between 
#' the individual's activity center (s) and the trap location. All coordinates (\code{s} and \code{trapCoords}) should be scaled to the habitat (see \code{\link{scaleCoordsToHabitatGrid}})
#' 
#' The \code{dmultiLocal_normal} distribution incorporates three features to increase computation efficiency (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item A local evaluation of the detection probability calculation (see Milleret et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item A sparse matrix representation (\emph{x}, \emph{detIndices} and \emph{detNums}) of the observation data to reduce the size of objects to be processed.
#' \item An indicator (\emph{indicator}) to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' The \code{dmultiLocal_normal} distribution requires x- and y- detector coordinates (\emph{trapCoords}) and activity centers coordinates (\emph{s}) to be scaled to the habitat grid (\emph{habitatGrid}) using the (\code{\link{scaleCoordsToHabitatGrid}} function.)
#'
#' When the aim is to simulate detection data: 
#' \enumerate{
#' \item \emph{x} should be provided using the \emph{yCombined} object as returned by \code{\link{getSparseY}}, 
#' \item arguments \emph{detIndices} and \emph{detNums} should not be provided, 
#' \item argument \emph{lengthYCombined} should be provided using the \emph{lengthYCombined} object as returned by  \code{\link{getSparseY}}.
#' }
#' 
#' @name dmultiLocal_normal
#'
#'
#'
#' @param x Vector of individual detection frequencies. This argument can be provided in two formats: 
#' (i) with the \emph{y} object as returned by the \code{\link{getSparseY}} function;
#' (ii) with the \emph{yCombined} object as returned by \code{\link{getSparseY}}. 
#' Note that when the random generation functionality is used (\code{rmultiLocal_normal}), only the \emph{yCombined} format can be used. 
#' The \emph{yCombined} object combines \emph{detNums}, \emph{x}, and \emph{detIndices} (in that order).  When such consolidated representation of the detection data \emph{x} is used, \emph{detIndices} and \emph{detNums} arguments should not be specified.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param detIndices Vector of indices of traps where the detections in \emph{x} were recorded; from the \emph{detIndices} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when x is provided as the  \emph{yCombined} object (returned by \code{\link{getSparseY}} ) and when detection data are simulated.
#' @param detNums umber of traps with at least one detection recorded in \emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when the \emph{yCombined} object (returned by \code{\link{getSparseY}}) is provided as \emph{x} and when detection data are simulated.
#' @param size Number of occasions.
#' @param p0 Baseline detection probability (scalar) used in the half-normal detection function. For trap-specific baseline detection probabilities use argument \emph{p0Traps} (vector) instead.
#' @param p0Traps Vector of baseline detection probabilities for each trap used in the half-normal detection function. When \emph{p0Traps} is used, \emph{p0} should not be provided. 
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Individual activity center x- and y-coordinates scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}})).
#' @param trapCoords Matrix of x- and y-coordinates of all traps scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}})).
#' @param localTrapsIndices Matrix of indices of local traps around each habitat grid cell, as returned by the \code{\link{getLocalObjects}} function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the \code{\link{getLocalObjects}} function.
#' @param resizeFactor Aggregation factor used in the \code{\link{getLocalObjects}} function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid Matrix of local habitat grid cell indices, from \emph{habitatGrid} returned by the \code{\link{getLocalObjects}} function. 
#' @param indicator Binary argument specifying whether the individual is available for detection (indicator = 1) or not (indicator = 0).
#' @param lengthYCombined The length of the x argument when the (\emph{yCombined}) format of the detection data is provided;  from the \emph{lengthYCombined} object returned by \code{\link{getSparseY}}
#' 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / 2 \sigma^2)}.
#'
#' @author Soumen Dey, Cyril Milleret
#'
#' @import nimble
#'
#' @examples
#' # I. DATA SET UP 
#'coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
#'                                    1.5, 3.5,
#'                                    2.5, 3.5,
#'                                    3.5, 3.5,
#'                                    0.5, 2.5,
#'                                    1.5, 2.5,
#'                                    2.5, 2.5,
#'                                    3.5, 2.5,
#'                                    0.5, 1.5,
#'                                    1.5, 1.5,
#'                                    2.5, 1.5,
#'                                    3.5, 1.5,
#'                                    0.5, 0.5,
#'                                    1.5, 0.5,
#'                                    2.5, 0.5,
#'                                    3.5, 0.5), ncol=2,byrow = TRUE)
#'colnames(coordsHabitatGridCenter) <- c("x","y")
#'# CREATE OBSERVATION WINDOWS
#'trapCoords <- matrix(c(1.5, 1.5, 2.5, 1.5, 1.5, 2.5, 2.5, 2.5), nrow = 4, byrow = TRUE)
#'colnames(trapCoords) <- c("x","y")
#'# PLOT CHECK
#'plot(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 
#'points(trapCoords[,"y"]~trapCoords[,"x"],col="red",pch=16) 
#'
#'# PARAMETERS
#'p0 <- 0.2
#'sigma <- 2
#'indicator <- 1 
#'# WE CONSIDER 2 INDIVIDUALS
#'
#'y <- matrix(c(2, 1, 1,-1, 1, 4, -1,#id#1 detected 2 times at detector 1 and 4 
#'              2, 2,-1,-1, 3,-1, -1),#id#2 detected 2 times at detector 3
#'               ncol=7, nrow=2, byrow = TRUE)
#'
#'
#'
#'
#'s <- matrix(c(0.5, 1,
#'              1.6, 2.3),ncol=2,nrow=2)
#'
#'# RESCALE COORDINATES 
#'ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
#'                                             coordsHabitatGridCenter = coordsHabitatGridCenter)
#'ScaledtrapCoords<- ScaledtrapCoords$coordsDataScaled
#'habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
#'
#'
#'# CREATE LOCAL OBJECTS 
#'TrapLocal <- getLocalObjects(habitatMask = habitatMask,
#'                                   coords = ScaledtrapCoords,
#'                                   dmax=2.5,
#'                                   resizeFactor = 1,
#'                                   plot.check = TRUE
#')
#'
#'# GET SPARSE MATRIX 
#'#SparseY <- getSparseY(y)
#'
#'# II. USING THE DENSITY FUNCTION 
#' # WE TAKE THE FIRST INDIVIDUAL
#'i=1
#'  # OPTION 1: USING THE RANDOM GENERATION FUNCTIONNALITY 
#'dmultiLocal_normal(x=y[i,]
#'                   ,
#'                   size=3
#'                   ,
#'                   p0 = p0
#'                   ,
#'                   sigma= sigma
#'                   , 
#'                   s=s[i,1:2]
#'                   ,
#'                   trapCoords=ScaledtrapCoords
#'                   ,
#'                   localTrapsIndices=TrapLocal$localIndices
#'                   ,
#'                   localTrapsNum=TrapLocal$numLocalIndices
#'                   ,
#'                   resizeFactor=TrapLocal$resizeFactor
#'                   ,
#'                   habitatGrid=TrapLocal$habitatGrid
#'                   ,
#'                   indicator=indicator
#'                   ,
#'                   lengthYCombined = ncol(y)
#'                   )
#'                                                                
#'  
#'
#'# III. USING THE RANDOM GENERATION FUNCTION 
#'rmultiLocal_normal(n=1,
#'                   size=3,
#'                   p0 = p0,
#'                   sigma= sigma, 
#'                   s=s[i,1:2],
#'                   trapCoords=ScaledtrapCoords,
#'                   localTrapsIndices=TrapLocal$localIndices,
#'                   localTrapsNum=TrapLocal$numLocalIndices,
#'                   resizeFactor=TrapLocal$resizeFactor,
#'                   habitatGrid=TrapLocal$habitatGrid,
#'                   indicator=indicator,
#'                   lengthYCombined = ncol(y))
#'                   
#' @export
NULL

#' @rdname dmultiLocal_normal
#' @export
dmultiLocal_normal <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(0),
                  p0 = double(0, default = -999),
                  p0Traps = double(1),
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
      if(detNums==0){
        NumOccNotDetected <- size
      }else{
        NumOccNotDetected <- size-sum(x[2:((detNums)+1)])
      }
    }else{
      x1 <- x
      detIndices1 <- detIndices
      NumOccNotDetected <- size-sum(x1)
      
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
    detIndices1 <- c(detIndices1, 0)
    count <- 1 
    
    # With dmulti, we need to recreate the entire p and y vector 
    p <- numeric(length = localTrapsNum[sID])
    y <- numeric(length = localTrapsNum[sID])
    
    
    if(p0==-999){# when p0 is provide through p0Traps
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p[r] <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          # p[r] <- -log(1-p[r])
          y[r] <- x1[count]
          count <- count + 1
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p[r] <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          # p[r] <- -log(1-p[r])
          y[r] <- 0
        }
      }
    }else{# when p0 is provide through p0
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p[r] <- p0 * exp(alpha * d2)
          # p[r] <- -log(1-p[r])
          y[r] <- x1[count]
          count <- count + 1 
        }else{
          d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          p[r] <- p0 * exp(alpha * d2)
          # p[r] <- -log(1-p[r])
          y[r] <- 0
        }
      }
    }
    y <- c(y, NumOccNotDetected)
    
    sum.p <- sum(p)
    risk <- 1 - exp(- sum.p)
    p <- risk * p / sum.p
    p <- c(p,  1 - risk)
    logProb <- dmulti(y, prob = p, size = size, log = TRUE)
    
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dmultiLocal_normal
#' @export
rmultiLocal_normal <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(0),
                  p0 = double(0, default = -999),
                  p0Traps = double(1),
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
    if(detNums >= 0) stop("Random generation for the rmultiLocal_normal distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    
    #========================================================
    # RETURN TYPE DECLARATION
    if(n!=1){print("rmultiLocal_normal only allows n = 1; using n = 1")}
    ## GET NECESSARY INFO
    alpha <- -1.0 / (2.0 * sigma * sigma)
    nMAxDetections <- (lengthYCombined-1)/2
    ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
    if(indicator == 0){
      y0 <- rep(0.0, lengthYCombined)
      #y0[2] <- size 
      return(y0)
    }
    
    ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
    sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
    theseLocalTraps <- localTrapsIndices[sID, 1:localTrapsNum[sID]]
    
    ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS
    detectOut <- rep(0, localTrapsNum[sID])
    ys <- rep(-1, nMAxDetections)
    dets <- rep(-1, nMAxDetections)
    count <- 1
    
    # With dmulti, we need to recreate the entire p and y vector 
    p <- numeric(length = localTrapsNum[sID])
    y <- numeric(length = localTrapsNum[sID])
    
    
    ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
    if(p0==-999){## when p0 is provided through p0Traps
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p[r] <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
        # p[r] <- -log(1-p[r])
      }#r 
    }else{## when p0 is provided through p0
      for(r in 1:localTrapsNum[sID]){
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p[r] <- p0 * exp(alpha * d2)
        # p[r] <- -log(1-p[r])
      }#r 
      
    }

    sum.p <- sum(p)
    risk <- 1 - exp(- sum.p)
    p <- risk * p / sum.p
    p <- c(p,  1 - risk)
    
    # p <- c(p, (prod(1-p)))
    y <- rmulti(n=1, size=size, prob=p)
    
    
    # reshape the data to get a sparse representation
    whichDetetecions <- y[1:localTrapsNum[sID]]>0
    totalNumDetections <- sum(whichDetetecions)
    
    if(nMAxDetections<totalNumDetections){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                           You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
    
    out <- rep(-1, lengthYCombined)
    out[1] <- totalNumDetections
    if(totalNumDetections >= 1){
      out[2:(totalNumDetections+1)] <- y[c(whichDetetecions, FALSE)]#[whichDetetecions]
      out[(nMAxDetections+2):(nMAxDetections+totalNumDetections+1)] <- theseLocalTraps[whichDetetecions]
    }
    ## OUTPUT
    return(out)
  })


