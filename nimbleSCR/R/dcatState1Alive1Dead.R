#' Density and random generation of a categorical distribution describing state transition with one alive and one dead states.
#' 
#'
#' The \code{dcatState1Alive1Dead} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' individual state transition. This function can be used to model transitions from one alive and one dead state. 
#' If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability prob1To2_t, so z_{i,t+1} ~ dcat(1- prob1To2_t, prob1To2_t, 0,0 , 0) where prob1To2_t represent the probability of an unborn individual to be recruited.
#' If z_{i,t} = 2, individual i can die and transition to z_{i,t+1}=3 with probability prob2To3, or survive with probability 1-prob2To3
#' Individuals in dead states (z_{i,t} = 3 ) remain in that state with probability 1, the absorbing state.
#' If transition probabilities are spatially variable, a probability vector containing the transition probability value in each habitat window can be provided using the "Hab" arguments (e.g. prob1To2Hab,prob2To3Hab).
#' 
#' @name dcatState1Alive1Dead 
#' 
#' @param x Scalar, individual state z_{i,t+1}.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param z Scalar, initial individual state z_{i,t}.
#' @param s Vector of x- and y-coordinates corresponding to the AC location of the individual. Used to extract transition spatially-explicit probabilities when they are provided.
#' @param prob1To2 scalar, probability to transition from state 1 to 2.
#' @param prob1To2Hab vector, Spatially-explicit probability to transition from state 2 to 3.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param prob2To3 scalar, probability to transition from state 2 to 3.
#' @param prob2To3Hab vector, Spatially-explicit probability to transition from state 2 to 3.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the 
#' order of habitat windows in \code{prob1To2Hab}, \code{prob2To3Hab} and in \code{lowerCoords}  and \code{upperCoords} as used in the \code{dbernppAC} function.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' \code{dcatState1Alive1Dead} gives the (log) probability density of \code{x}. 
#' \code{rcatState1Alive2Dead} gives a randomly generated individual states conditional on the initial state \code{z}.  
#' 
#' @author Cyril Milleret
#' 
#' @import nimble
#' @examples
#' # Use the distribution in R
#' z <- 2
#' prob1To2 <- 0.2
#' prob2To3 <- 0.7
#' 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(rep(1,4))
#' logSumIntensity <- log(sum(c(1:4))) 
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' s <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'                habitatGrid, numGridRows, numGridCols)
#' 
#' ## No spatial mortality 
#' zPlusOne <- rcatState1Alive1Dead( z = z
#'                                  , prob1To2 = prob1To2
#'                                  , prob2To3 = prob2To3
#'                                  , s = s
#'                                  , habitatGrid = habitatGrid)
#' zPlusOne     
#' 
#' dcatState1Alive1Dead(  x = zPlusOne
#'                      , z = z
#'                      , prob1To2 = prob1To2
#'                      , prob2To3 = prob2To3
#'                      , s = s
#'                      , habitatGrid = habitatGrid)
#'     
#' ##  With spatial mortality
#' prob2To3Hab <- c(0.60, 0.70, 0.74, 0.65)
#' prob1To2Hab <- c(0.4,0.5,0.1,0.3)
#' zPlusOne <- rcatState1Alive1Dead( z = z
#'                                 , prob1To2Hab = prob1To2Hab
#'                                 , prob2To3Hab = prob2To3Hab
#'                                 , s = s
#'                                 , habitatGrid = habitatGrid)
#' zPlusOne    
#' dcatState1Alive1Dead(  x = zPlusOne
#'                      , z = z
#'                      , prob1To2Hab = prob1To2Hab
#'                      , prob2To3Hab = prob2To3Hab
#'                      , s = s
#'                      , habitatGrid = habitatGrid)
#' 
#' 
#' 
#' 
#' 
NULL
#' @rdname dcatState1Alive1Dead
#' @export



#### 1.Density function ####
dcatState1Alive1Dead  <- nimbleFunction(run = function( x = double(0)
                                                  , z = double(0)
                                                  , prob1To2 = double(0, default = -999)
                                                  , prob1To2Hab = double(1)
                                                  , prob2To3 = double(0, default = -999)
                                                  , prob2To3Hab = double(1)
                                                  , s = double(1)
                                                  , habitatGrid = double(2)
                                                  , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    
    if(prob1To2 == -999){
      sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
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
    
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
   
    #prob2To3
    if(prob2To3 == -999){
      indProb2To3 <- prob2To3Hab[sID]
    }else{
      indProb2To3 <- prob2To3
    }
    
    logLikelihood <- dcat(x, prob = c(0, 1-indProb2To3, indProb2To3), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3){
    logLikelihood <- dcat(x, prob = c(0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})


NULL
#' @rdname dcatState1Alive1Dead
#' @export
#' 
#### 2.Sampling function ####
rcatState1Alive1Dead <- nimbleFunction(run = function( n = integer(0)
                                                     , z = double(0)
                                                     , prob1To2 = double(0, default = -999)
                                                     , prob1To2Hab = double(1)
                                                     , prob2To3 = double(0, default = -999)
                                                     , prob2To3Hab = double(1)
                                                     , s = double(1)
                                                     , habitatGrid = double(2)
){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    if(prob1To2 == -999){
      sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
      indProb1To2 <- prob1To2Hab[sID]
    }else{
      ## EXTRACT LOCATION OF THE ID
      indProb1To2 <- prob1To2
    }
    state <- rcat(1, prob = c(1 - indProb1To2, indProb1To2, 0))
    return(state)
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    #prob2To3
    if(prob2To3 == -999){
      indProb2To3 <- prob2To3Hab[sID]
    }else{
      indProb2To3 <- prob2To3
    }
    
    state <- rcat(1, prob = c(0, 1- indProb2To3, indProb2To3))
    return(state)
  }
  
  if(z == 3 ){
    state <- 3
    return(state)
  }
  
})




