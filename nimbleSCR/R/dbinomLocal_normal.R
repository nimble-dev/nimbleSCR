#' Local evaluation of a binomial SCR observation process 
#'
#' The \code{dbinomLocal_normal} distribution is a NIMBLE custom distribution which can be used to model 
#' the binomial observations (x) of a single individual over a set of detectors defined by their 
#' coordinates (trapCoords). The distribution assumes that the detection probability at any detector 
#' follows a half-normal function of the distance between the individual's activity center (s) and the detector location.
#'
#' The \code{dbinomLocal_normal} distribution incorporates three features to increase computation efficiency:
#' \enumerate{
#' \item A local evaluation of the detection probability calculation (see Milleret et al. (2019) <doi:10.1002/ece3.4751> for more details).
#' \item It uses a sparse matrix representation (x, detIndices, detNums) of the observation data to reduce the size of objects to be processed.
#' \item It uses an indicator (indicator) to shortcut calculations for individuals unavailable for detection.
#' }
#' 
#' The \code{dbinomLocal_normal} distribution requires that the x- and y- coordinates should be scaled to the habitat (\code{\link{scaleCoordsToHabitatGrid}})
#' 
#' @name dbinomLocal_normal
#'
#' @param x Vector of individual detection frequencies, as returned by the \code{getSparseY} function 
#' (padded with -1's to maintain the square structure of the observation data).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param detIndices Vector of the detector indices where the detections in x were recorded, as returned by the \code{getSparseY} function.
#' @param detNums Number of detections recorded in x, as returned by the \code{getSparseY} function.
#' @param size Vector of the number of trials (zero or more) for each trap (trapCoords).
#' @param p0 Baseline detection probability used in the half-normal detection function.
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Individual activity center x- and y-coordinates.
#' @param trapCoords Matrix of x- and y-coordinates of all traps.
#' @param localTrapsIndices Matrix of indices of local traps around each habitat grid cell, as returned by the \code{getLocalObjects} function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the getLocalObjects function.
#' @param resizeFactor Aggregation factor used in the \code{getLocalObjects} function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid Matrix of habitat grid cells indices.
#' @param indicator Logical argument, specifying whether the individual is available for detection.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#'
#' @return The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
#'  and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / \sigma^2)}.
#'
#' @author Cyril Milleret
#'
#' @import nimble
#' @importFrom stats dbinom
#'
#' @examples
#' ## define model code
#' code <- nimbleCode({
#'     psi ~ dunif(0,1)
#'     p0 ~ dunif(0,1)
#'     sigma ~ dunif(0,100)
#'     N <- sum(z[1:M])
#'     for(i in 1:M) {
#'         s[i, 1] ~ dunif(0, 100)
#'         s[i, 2] ~ dunif(0, 100)
#'         z[i] ~ dbern(psi)
#'         y[i,1:maxDetNum] ~ dbinomLocal_normal(detNums,
#'                                                  detIndices,
#'                                                  size,
#'                                                  p0,
#'                                                  sigma,
#'                                                  s[i,1:2],
#'                                                  trapCoords,
#'                                                  localTrapsIndices,
#'                                                  localTrapsNum,
#'                                                  resizeFactor,
#'                                                  habitatGrid,
#'                                                  z[i])
#'     }
#' })
#'  
#' ## create NIMBLE model object
#' ## Rmodel <- nimbleModel(code, ...)
#'  
#' ## use model object for MCMC, etc.
#'
#' @export
NULL

#' @rdname dbinomLocal_normal
#' @export
dbinomLocal_normal <- nimbleFunction(
  run = function( x = double(1),
                  detNums = double(0),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0, default = 1.0),
                  log = integer(0, default = 0)
  ) {
    ## Specify return type
    returnType(double(0))
    
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
        if(sum(detIndices[r] == theseLocalTraps) == 0){
          if(log == 0) return(0.0)
          else return(-Inf)
        }
      }
    }
    
    ## Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logProb <- 0.0 
    detIndices1 <- c(detIndices,0)
    count <- 1 
    
    for(r in 1:localTrapsNum[sID]){
      if(theseLocalTraps[r] == detIndices1[count]){ 
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p <- p0 * exp(alpha * d2)
        logProb <- logProb + dbinom(x[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
        count <- count + 1
      }else{
        d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
        p <- p0 * exp(alpha * d2)
        logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
        
      }
    }
    
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dbinomLocal_normal
#' @export
rbinomLocal_normal <- nimbleFunction(
  run = function( n = double(0, default = 1),
                  detNums = double(0),
                  detIndices = double(1),
                  size = double(1),
                  p0 = double(0),
                  sigma = double(0),
                  s = double(1),
                  trapCoords = double(2),
                  localTrapsIndices = double(2),
                  localTrapsNum = double(1),
                  resizeFactor = double(0, default = 1),
                  habitatGrid = double(2),
                  indicator = double(0, default = 1.0)
  ) {
    ## Specify return type
    returnType(double(1))
    if(n >= 0) stop("Random generation for the dbinomLocal_normal distribution is not currently supported")
    dummyOut <- detIndices
    return(dummyOut)
  })



