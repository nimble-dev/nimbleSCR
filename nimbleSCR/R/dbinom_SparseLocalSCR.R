#' Local evaluation of a binomial SCR observation process 
#'
#' The dbinom_SparseLocalSCR distribution is a NIMBLE custom distribution which can be used to model 
#' the binomial observations (x) of a single individual over a set of detectors defined by their 
#' coordinates (trapCoords). The distribution assumes that the detection probability at any detector 
#' follows a half-normal function of the distance between the individual's activity center (s) and the detector location.
#'
#' The dbinom_SparseLocalSCR distribution incorporates three features to increase computation efficiency:
#' 1 - A local evaluation of the detection probability calculation (see Milleret et al., 2019 for more details)
#' 2 - It uses a sparse matrix representation (x, detIndices, detNum) of the observation data to reduce the size of objetcs to be processed.
#' 3 - It uses an indicator (indicator) to shortcut calculations for individuals unavailable for detection.
#'
#' @name dbinom_sparseLocalSCR
#'
#' @param x vector of individual detection frequencies, as returned by the getSparseY function 
#' (padded with -1's to maintain the square structure of the observation data).
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param detIndices vector of the detector indices where the detections in x were recorded, as returned by the getSparseY function.
#' @param detNums number of detections recorded in x, as returned by the getSparseY function.
#' @param size vector of number of trials (zero or more).
#' @param p0 Baseline detection probability used in the half-normal detection function.
#' @param sigma Scale parameter of the half-normal detection function.
#' @param s Bivariate activity center coords.
#' @param trapCoords matrix of x- and y-coordinate locations of all traps
#' @param localTrapsIndices Matrix of indices of local traps around all habitat grid cells, as returned by the getLocalTraps function.
#' @param localTrapsNum  Vector of numbers of local traps around all habitat grid cells, as returned by the getLocalTraps function.
#' @param resizeFactor aggregation factor used in the getLocalTraps function to reduce the number of habitat grid cells to retrieve local traps for.
#' @param habitatGrid matrix of habitat grid cells indices.
#' @param indicator Logical argument, specifying whether the current individual is available or not for detection.
#'
#' @return The log-likelihood value associated with the vector of SCR observations, given the current activity center (s),
#'  and the half-normal detection function : p = p0 * exp(-d^2 / sigma^2).
#'
#' @author Cyril Milleret
#'
#' @import nimble
#' @importFrom stats dbinom  
#'
#' @examples
#' \donttest{
#' ## define model code
#' # nimModel <- nimbleCode({
#' 
#' psi ~ dunif(0,1)
#' p0 ~ dunif(0,1)
#' sigma ~ dunif(0,100)
#' 
#' N <- sum(z[1:M])
#'
#'    for(i in 1:M) {
#'         s[i, 1] ~ dunif(0, 100)
#'         s[i, 2] ~ dunif(0, 100)
#'         
#'         z[i] ~ dbern(psi)
#'         
#'         y[i,1:maxDetNum] ~ dbinom_sparseLocalSCR( detNums = double(0),
#'                                                   detIndices = double(1),
#'                                                   size = double(1),
#'                                                   p0 = double(0),
#'                                                   sigma = double(0),
#'                                                   s = double(1),
#'                                                   trapCoords = double(2),
#'                                                   localTrapsIndices = double(2),
#'                                                   localTrapsNum = double(1),
#'                                                   resizeFactor = double(0, default = 1),
#'                                                   habitatGrid = double(2),
#'                                                   indicator = double(0, default = 1.0)
#' })
#'
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code)
#'
#' ## use model object for MCMC, etc.
#' }
#'
#' @export
NULL

#' @rdname dbinom_sparseLocalSCR
#' @export
dbinom_sparseLocalSCR <- nimbleFunction(
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
      if(detNum == 0){
        if(log == 0) return(1.0)
        else return(0.0)
      } else {
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
    
    ## Retrieve the index of the habitat cell where the current AC is
    sID <- habitatID[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
    
    ## Retrieve the indices of the local traps surrounding the selected habita grid cell
    theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
    
    ## Retrieve the number of local traps surrounding the selected habita grid cell
    numLocalTraps <- length(theseLocalTraps)
    
    ## Retrieve the total number of traps
    nDetectors <- length(size)
    
    ## Recreate the full detection vector (NOT NECESSARY ; COULD BE OPTIMIZED IN A LATER VERSION)
    y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
    if(detNum > 0){
      for(r in 1:detNum){
        y[detIndices[r]] <- x[r] 
        if(sum(detIndices[r] == index) == 0){
          if(log == 0) return(0.0)
          else return(-Inf)
        } 
      } 
    }
    
    ## Calculate the log-probability of the vector of detections
    alpha <- -1.0 / (2.0 * sigma * sigma)
    logProb <- 0.0 
    count <- 1 
    index1 <- c(index,0) 
    for(r in 1:nDetectors){
      if(index1[count] == r){ 
        d2 <- pow(trapCoords[r,1] - s[1], 2) + pow(trapCoords[r,2] - s[2], 2)
        p <- pZero * exp(alpha * d2)
        logProb <- logProb + dbinom(y[r], prob = p, size = trials[r], log = TRUE)
        count <- count + 1
      }
    }
    
    ## Return the probability of the vector of detections (or log-probability if required)
    if(log)return(logProb)
    return(exp(logProb))
  })


#' @rdname dbinom_sparseLocalSCR
#' @export
#' 
rbinom_sparseLocalSCR <- nimbleFunction(
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
    stop("Random generation for the dbinom_sparseLocalSCR distribution is not currently supported")
  })