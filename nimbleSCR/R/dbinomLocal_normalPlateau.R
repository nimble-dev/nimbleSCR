#' Local evaluation of a binomial SCR observation process 
#'
#' The \code{dbinomLocal_normalPlateau} distribution is a NIMBLE custom distribution which can be used to model 
#' and simulate binomial observations (x) of a single individual over a set of traps defined by their coordinates \emph{trapCoords}
#' the distribution assumes that an individual’s detection probability at any trap follows a half-normal plateau function of the distance between 
#' the individual's activity center (s) and the trap location. With the half-normal plateau function, detection probability remains constant with value 
#' p0 for a plateau of width w before declining with scale sigma. 
#' 
#' All coordinates (\code{s} and \code{trapCoords}) should be scaled to the habitat (see \code{\link{scaleCoordsToHabitatGrid}})
#' 
#'
#'
#' The \code{dbinomLocal_normalPlateau} distribution incorporates three features to increase computation efficiency (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details):
#' \enumerate{
#' \item A local evaluation of the detection probability calculation (see Milleret et al., 2019 <doi:10.1002/ece3.4751> for more details)
#' \item A sparse matrix representation (\emph{x}, \emph{detIndices} and \emph{detNums}) of the observation data to reduce the size of objects to be processed.
#' \item An indicator (\emph{indicator}) to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' The \code{dbinomLocal_normalPlateau} distribution requires x- and y- detector coordinates (\emph{trapCoords}) to be scaled to the habitat grid (\emph{habitatGrid}) using the (\code{\link{scaleCoordsToHabitatGrid}} function.)
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
#' @name dbinomLocal_normalPlateau
#'
#' @param x Vector of individual detection frequencies. This argument can be provided in two formats: (i) with the \emph{y} object as returned by \code{\link{getSparseY}}; (ii) with the
#' \emph{yCombined} object as returned by \code{\link{getSparseY}} Note that when the random generation functionality is used (rbinomLocal_normal), only the yCombined format can be used. 
#' The  \emph{yCombined}  object combines \emph{detNums}, \emph{x}, and \emph{detIndices} (in that order). When such consolidated 
#' representation of the detection data x is used, \emph{detIndices} and \emph{detNums} arguments should not be specified.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param detIndices Vector of indices of traps where the detections in \emph{x} were recorded; from the \emph{detIndices} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when x is provided as the  \emph{yCombined} object (returned by \code{\link{getSparseY}} ) and when detection data are simulated.
#' @param detNums Number of traps with at least one detection recorded in \emph{x}; from the \emph{detNums} object returned by the \code{\link{getSparseY}} function. 
#' This argument should not be specified when the \emph{yCombined} object (returned by \code{\link{getSparseY}}) is provided as \emph{x} and when detection data are simulated.
#' @param size Vector of the number of trials (zero or more) for each trap (\emph{trapCoords}).
#' @param p0 Baseline detection probability (scalar) used in the half-normal detection function. For trap-specific baseline detection probabilities use argument \emph{p0Traps} (vector) instead.
#' @param p0Traps Vector of baseline detection probabilities for each trap used in the half-normal detection function. When \emph{p0Traps} is used, \emph{p0} should not be provided. 
#' @param sigma Scale parameter of the half-normal detection function.
#' @param w Length of plateau of the half-normal plateau detection function.
#' @param s Individual activity center x- and y-coordinates scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
#' @param trapCoords Matrix of x- and y-coordinates of all traps scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).
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
#'  and the half-normal plateau detection function : \eqn{p = p0} when \emph{d < w} and \eqn{p = p0 * exp(-(d-w)^2 / \sigma^2)} when \emph{d >= w}.
#'
#' @author Soumen Dey
#'
#' @references
#' 
#' Dey, S., Bischof, R., Dupont, P. P. A., & Milleret, C. (2022). Does the punishment fit the crime? Consequences and diagnosis of misspecified detection functions in Bayesian spatial capture–recapture modeling. Ecology and Evolution, 12, e8600. https://doi.org/10.1002/ece3.8600
#' 
#' @import nimble
#' @importFrom stats dbinom
#'
#' @examples
#' # A user friendly vignette is also available on github: 
#' # https://github.com/nimble-dev/nimbleSCR/blob/master/nimbleSCR/vignettes/
#' # Vignette name: Fit_with_dbinomLocal_normalPlateau_and_HomeRangeAreaComputation.rmd
#' 
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
#' p0 <- 0.25
#' sigma <- 1
#' w <- 1.5
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
#' dbinomLocal_normalPlateau(x=SparseY$y[i,,1],
#'                    detNums=SparseY$detNums[i],
#'                    detIndices=SparseY$detIndices[i,,1],
#'                    size=rep(1,4),
#'                    p0 = p0,
#'                    sigma= sigma, 
#'                    w = w,
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
#' dbinomLocal_normalPlateau(x=SparseY$yCombined[i,,1],
#'                    size=rep(1,4),
#'                    p0 = p0,
#'                    sigma= sigma, 
#'                    w = w,
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
#' rbinomLocal_normalPlateau(n=1,
#'                    size=rep(1,4),
#'                    p0 = p0,
#'                    sigma= sigma, 
#'                    w = w,
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

#' @rdname dbinomLocal_normalPlateau
#' @export
dbinomLocal_normalPlateau <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0, default = -999),
                  p0Traps = double(1),
                  sigma = double(0),
                  w = double(0, default = 2.0),
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
    detIndices1 <- c(detIndices1, 0)
    count <- 1 
    
    
    if(p0==-999){# when p0 is provided through p0Traps
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
          if(d < w) p <- p0Traps[theseLocalTraps[r]]
          if(d >=  w) p <- p0Traps[theseLocalTraps[r]] * exp(alpha * (d-w)*(d-w))
          # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          # p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          count <- count + 1
        }else{
          d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
          if(d < w) p <- p0Traps[theseLocalTraps[r]]
          if(d >=  w) p <- p0Traps[theseLocalTraps[r]] * exp(alpha * (d-w)*(d-w))
          # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          # p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
          logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          
        }
      }
    }else{# when p0 is provided through p0
      for(r in 1:localTrapsNum[sID]){
        if(theseLocalTraps[r] == detIndices1[count]){ 
          d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
          if(d < w) p <- p0
          if(d >=  w) p <- p0 * exp(alpha * (d-w)*(d-w))#/(2.0*sigma*sigma))
          # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          # p <- p0 * exp(alpha * d2)
          logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
          count <- count + 1
        }else{
          d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
          if(d < w) p <- p0
          if(d >=  w) p <- p0 * exp(alpha * (d-w)*(d-w))#/(2.0*sigma*sigma))
          # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
          # p <- p0 * exp(alpha * d2)
          logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
        }
      }
    }
    
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dbinomLocal_normalPlateau
#' @export
rbinomLocal_normalPlateau <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0, default = -999),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0, default = -999),
                  p0Traps = double(1),
                  sigma = double(0),
                  w = double(0, default = 2.0),
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
    if(detNums >= 0) stop("Random generation for the rbinomLocal_normalPlateau distribution is not currently supported without combining all individual detections information in one vector. See 'getSparseY()'")
    
    #========================================================
    # RETURN TYPE DECLARATION
    if(n!=1){print("rbinomLocal_normalPlateau only allows n = 1; using n = 1")}
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
    if(p0==-999){## when p0 is provided through p0Traps
      for(r in 1:localTrapsNum[sID]){
        d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
        if(d < w) p <- p0Traps[theseLocalTraps[r]]
        if(d >=  w) p <- p0Traps[theseLocalTraps[r]] * exp(alpha * (d-w)*(d-w))
        # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        # p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
        # Draw the observation at detector j from a binomial distribution with probability p
        detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], p)
        if(detectOut[r] >0){
          if(nMAxDetections<count){stop("Simulated individual detections occur at more traps than what can be stored within x.\n
                                          You may need to augment the size of the x object with the argument 'nMaxTraps' from the getSparseY() function")}
          ys[count] <- detectOut[r]
          dets[count] <- theseLocalTraps[r]
          count <- count + 1
        }#if
      }#r 
    }else{## when p0 is provided through p0
      for(r in 1:localTrapsNum[sID]){
        d <- pow(pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2),0.5)
        if(d < w) p <- p0
        if(d >=  w) p <- p0 * exp(alpha * (d-w)*(d-w))#/(2.0*sigma*sigma))
        # d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        # p <- p0 * exp(alpha * d2)
        # Draw the observation at detector j from a binomial distribution with probability p
        detectOut[r] <- rbinom(1, size[theseLocalTraps[r]], p)
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
