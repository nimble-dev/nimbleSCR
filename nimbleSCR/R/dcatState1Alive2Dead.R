#' Density and random generation of a categorical distribution describing state transition with one alive and two dead states.
#' 
#'
#' The \code{dcatState1Alive2Dead} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' individual state transition. This function can be used to model transitions from one alive and two dead states. 
#' If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability prob1To2_t, so z_{i,t+1} ~ dcat(1- prob1To2_t, prob1To2_t, 0,0 , 0) where prob1To2_t represent the probability of an unborn individual to be recruited.
#' If z_{i,t} = 2, individual i can die from one cause of mortality (e.g. culling) and transition to z_{i,t+1}=3 with probability prob2To3, or die from another cause with probability prob2To4 z_{i,t+1}=4. 
#' If the individual does not die it can survive and remain in state 2 with probability (1-(prob2To3+prob2To4)). 
#' Individuals in dead states (z_{i,t} = 3 or 4) transition to z_{i,t+1} = 4, the absorbing state, with probability 1.
#' If transition probabilities are spatially variable, a probability vector containing the transition probability value in each habitat window can be provided using the "Hab" arguments (e.g. prob1To2Hab,prob2To3Hab).
#' 
#' @name dcatState1Alive2Dead 
#' 
#' @param x Scalar, individual state z_{i,t+1}.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param z Scalar, initial individual state z_{i,t}.
#' @param s Vector of x- and y-coordinates corresponding to the AC location of the individual. Used to extract transition spatially-explicit probabilities when they are provided.
#' @param prob1To2 scalar, probability to transition from state 1 to 2.
#' @param prob1To2Hab vector, Spatially-explicit probability to transition from state 2 to 3.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param prob2To3 scalar, probability to transition from state 2 to 3.
#' @param prob2To3Hab vector, Spatially-explicit probability to transition from state 2 to 3.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param prob2To4 scalar, probability to transition from state 2 to 4.
#' @param prob2To4Hab vector, Spatially-explicit probability to transition from state 2 to 4.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the 
#' order of habitat windows in \code{prob1To2Hab}, \code{prob2To3Hab} and in \code{lowerCoords}  and \code{upperCoords} as used in the \code{dbernppAC} function.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @return 
#' \code{dcatState1Alive2Dead} gives the (log) probability density of \code{x}. 
#' \code{rcatState1Alive2Dead} gives a randomly generated individual states conditional on the initial state \code{z}.  
#' 
#' @author Cyril Milleret
#' 
#' 
#' @examples
#' # Use the distribution in R
#' 
#' z <- 2
#' prob1To2 <- 0.2
#' prob2To3 <- 0.4
#' prob2To4 <- 0.1
#' 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(rep(1,4))
#' logSumIntensity <- log(sum(c(1:4))) 
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' s <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'                  habitatGrid, numGridRows, numGridCols)
#' 
#' ## No spatial mortality 
#' zPlusOne <- rcatState1Alive2Dead( z = z
#'                                  , prob1To2 = prob1To2
#'                                  , prob2To3 = prob2To3
#'                                  , prob2To4 = prob2To4
#'                                  , s = s
#'                                  , habitatGrid = habitatGrid)
#'    
#' dcatState1Alive2Dead(  x = zPlusOne
#'                      , z = z
#'                      , prob1To2 = prob1To2
#'                      , prob2To3 = prob2To3
#'                      , prob2To4 = prob2To4
#'                      , s = s
#'                      , habitatGrid = habitatGrid)
#'      
#' ##  With spatial mortality
#' prob2To3Hab <- c(0.10, 0.20, 0.15, 0.30)
#' prob2To4Hab <- c(0.13, 0.21, 0.12, 0.08)
#' phiSpatial <- 1-(prob2To3Hab+prob2To4Hab)
#' zPlusOne <- rcatState1Alive2Dead( z = z
#'                                  , prob1To2Hab = prob1To2Hab
#'                                  , prob2To3Hab =  prob2To3Hab
#'                                  , prob2To4Hab = prob2To4Hab
#'                                  , s = s
#'                                  , habitatGrid = habitatGrid)
#' 
#' dcatState1Alive2Dead(  x = zPlusOne
#'                      , z = z
#'                      , prob1To2Hab = prob1To2Hab
#'                      , prob2To3Hab =  prob2To3Hab
#'                      , prob2To4Hab = prob2To4Hab
#'                      , s = s
#'                      , habitatGrid = habitatGrid)
#' 
#' 
#' 
#' 
NULL
#' @rdname dcatState1Alive2Dead
#' @export



#### 1.Density function ####
dcatState1Alive2Dead  <- nimbleFunction(run = function( x = double(0)
                                                  , z = double(0)
                                                  , prob1To2 = double(0, default = -999)
                                                  , prob1To2Hab = double(1)
                                                  , prob2To3 = double(0, default = -999)
                                                  , prob2To3Hab = double(1)
                                                  , prob2To4 = double(0, default = -999)
                                                  , prob2To4Hab = double(1)
                                                  , s = double(1)
                                                  , habitatGrid = double(2)
                                                  , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    
    if(prob1To2 == -999){
      numGridRows <- nrow(habitatGrid)
      sID <- habitatGrid[numGridRows-trunc(s[2]), trunc(s[1])+1]
      indProb1To2 <- prob1To2Hab[sID]
      
    }else{
      ## EXTRACT LOCATION OF THE ID
      indProb1To2 <- prob1To2
    }
    
    logLikelihood <- dcat(x, prob = c(1 - indProb1To2, indProb1To2), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    
    numGridRows <- nrow(habitatGrid)
    sID <- habitatGrid[numGridRows-trunc(s[2]), trunc(s[1])+1]
    #prob2To3
    if(prob2To3== -999){
      indProb2To3 <- prob2To3Hab[sID]
        }else{
      indProb2To3 <- prob2To3
    }
    #w
    if(prob2To4 == -999){
      indProb2To4 <- prob2To4Hab[sID]
      }else{
        indProb2To4 <- prob2To4
    }
    #phi
    indProb2To2 <- 1 - ( indProb2To3 + indProb2To4)

    logLikelihood <- dcat(x, prob = c(0, indProb2To2, indProb2To3, indProb2To4), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3 | z == 4 ){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})


NULL
#' @rdname dcatState1Alive2Dead
#' @export
#' 
#### 2.Sampling function ####
rcatState1Alive2Dead <- nimbleFunction(run = function( n = integer(0)
                                                 , z = double(0)
                                                 , prob1To2 = double(0, default=-999)
                                                 , prob1To2Hab = double(1)
                                                 , prob2To3 = double(0, default=-999)
                                                 , prob2To3Hab = double(1)
                                                 , prob2To4 = double(0, default=-999)
                                                 , prob2To4Hab = double(1)
                                                 , s = double(1)
                                                 , habitatGrid = double(2)
){
  # Return type declaration
  returnType(double(0))

  if(z == 1){
    if(prob1To2 == -999){
      numGridRows <- dim(habitatGrid)[1]
      sID <- habitatGrid[numGridRows-trunc(s[2]), trunc(s[1])+1]
      indProb1To2 <- prob1To2Hab[sID]
    }else{
      ## EXTRACT LOCATION OF THE ID
      indProb1To2 <- prob1To2
    }
    state <- rcat(1, prob = c(1 - indProb1To2, indProb1To2))
    return(state)
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    numGridRows <- dim(habitatGrid)[1]
    sID <- habitatGrid[numGridRows-trunc(s[2]), trunc(s[1])+1]
    #prob2To3
    if(prob2To3== -999){
      indProb2To3 <- prob2To3Hab[sID]
    }else{
      indProb2To3 <- prob2To3
    }
    #w
    if(prob2To4 == -999){
      indProb2To4 <- prob2To4Hab[sID]
    }else{
      indProb2To4 <- prob2To4
    }
    #phi
    indProb2To2 <- 1 - ( indProb2To3 + indProb2To4)
    
    
    state <- rcat(1, prob = c(0, indProb2To2, indProb2To3, indProb2To4))
    return(state)
  }
  
  if(z == 3 | z == 4 ){
    state <- 4
    return(state)
  }
  
})




