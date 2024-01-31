#' Density and random generation of a categorical distribution describing state transition with two alive and two dead states.
#' 
#'
#' The \code{dcatState2Alive2Dead} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' individual state transition. This function can be used to model transitions from two alive and two dead states. 
#' If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability prob1To2_t, so z_{i,t+1} ~ dcat(1- prob1To2_t, prob1To2_t, 0,0 , 0) where prob1To2_t represent the probability of an unborn individual to be recruited.
#' If z_{i,t} = 2, individual i can die from one cause of mortality (e.g. culling) and transition to z_{i,t+1}=4 with probability prob2To4, or die from another cause with probability prob2To5 z_{i,t+1}=5. 
#' If the individual does not die (1-(prob2To4+prob2To5)), it can either transition to the second state alive (z_{i,t+1}=3) with probability prob2To3 or remain in the first state alive (z_{i,t+1}=2) with probability (1-prob2To3). 
#' If z_{i,t} = 3, individual i can die from one cause of mortality (e.g. culling) and transition to z_{i,t+1}=4 with probability prob3To4, or die from another cause with probability prob3To5 z_{i,t+1}=5. 
#' if the individual does not die (1-(prob3To4+prob3To5)), the individual remain in state 3. 
#' Individuals in dead states (z_{i,t} = 4 or 5) transition to z_{i,t+1} = 5, the absorbing state, with probability 1.
#' If transition probabilities are spatially variable, a probability vector containing the transition probability value in each habitat window can be provided using the "Hab" arguments (e.g. prob1To2Hab,prob2To3Hab).
#' 
#' 
#' @name dcatState2Alive2Dead 
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
#' @param prob2To5 scalar, probability to transition from state 2 to 5.
#' @param prob2To5Hab vector, Spatially-explicit probability to transition from state 2 to 5.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param prob3To4 scalar, probability to transition from state 3 to 4.
#' @param prob3To4Hab vector, Spatially-explicit probability to transition from state 3 to 4.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param prob3To5 scalar, probability to transition from state 3 to 5.
#' @param prob3To5Hab vector, Spatially-explicit probability to transition from state 3 to 5.  The length of the vector should be equal the number of habitat windows in \code{habitatGrid}.
#' @param habitatGrid Matrix of habitat window indices. Cell values should correspond to the 
#' order of habitat windows in \code{prob1To2Hab}, \code{prob2To3Hab} and in \code{lowerCoords}  and \code{upperCoords} as used in the \code{dbernppAC} function.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @return 
#' \code{dcatState2Alive2Dead} gives the (log) probability density of \code{x}. 
#' \code{dcatState2Alive2Dead} gives a randomly generated individual states conditional on the initial state \code{z}.  
#' 
#' @author Cyril Milleret
#' 
#' 
#' @import nimble
#' @examples
#' # Use the distribution in R
#' 
#'z <- 3
#'prob1To2 <- 0.2
#'prob2To3 <- 0.4
#'prob2To4 <- 0.1
#'prob2To5 <- 0.1
#'
#'prob3To4 <- 0.2
#'prob3To5 <- 0.1
#'
#'
#'lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#'upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#'logIntensities <- log(rep(1,4))
#'logSumIntensity <- log(sum(c(1:4))) 
#'habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#'numGridRows <- nrow(habitatGrid)
#'numGridCols <- ncol(habitatGrid)
#'s <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'               habitatGrid, numGridRows, numGridCols)
#'
#'## No spatial mortality 
#'zPlusOne <- rcatState2Alive2Dead( z = z
#'                                  , prob1To2 = prob1To2
#'                                  , prob2To3 = prob2To3
#'                                  , prob2To4 = prob2To4
#'                                  , prob2To5 = prob2To5
#'                                  , prob3To4 = prob3To4
#'                                  , prob3To5 = prob3To5
#'                                  , s = s
#'                                  , habitatGrid = habitatGrid)
#'
#'dcatState2Alive2Dead(  x = zPlusOne
#'                       , z = z
#'                       , prob1To2 = prob1To2
#'                       , prob2To3 = prob2To3
#'                       , prob2To4 = prob2To4
#'                       , prob2To5 = prob2To5
#'                       , prob3To4 = prob3To4
#'                       , prob3To5 = prob3To5
#'                       , s = s
#'                       , habitatGrid = habitatGrid)
#'
#'##  With spatial mortality
#'prob2To3Hab <- runif(length(habitatGrid),0,0.1)
#'prob2To4Hab <- runif(length(habitatGrid),0,0.1)
#'prob2To5Hab <- runif(length(habitatGrid),0,0.1)
#'prob3To4Hab <- runif(length(habitatGrid),0,0.1)
#'prob3To5Hab <- runif(length(habitatGrid),0,0.1)
#'
#'
#'
#'zPlusOne <- rcatState2Alive2Dead( z = z
#'                                  , prob1To2 = prob1To2
#'                                  , prob2To3Hab = prob2To3Hab
#'                                  , prob2To4Hab = prob2To4Hab
#'                                  , prob2To5Hab = prob2To5Hab
#'                                  , prob3To4Hab = prob3To4Hab
#'                                  , prob3To5Hab = prob3To5Hab
#'                                  , s = s
#'                                  , habitatGrid = habitatGrid)
#'
#'dcatState2Alive2Dead(  x = zPlusOne
#'                       , z = z
#'                       , prob1To2 = prob1To2
#'                       , prob2To3Hab = prob2To3Hab
#'                       , prob2To4Hab = prob2To4Hab
#'                       , prob2To5Hab = prob2To5Hab
#'                       , prob3To4Hab = prob3To4Hab
#'                       , prob3To5Hab = prob3To5Hab
#'                       , s = s
#'                       , habitatGrid = habitatGrid)
#' 
#' 
NULL
#' @rdname dcatState2Alive2Dead
#' @export



#### 1.Density function ####
dcatState2Alive2Dead  <- nimbleFunction(run = function( x = double(0)
                                                        , z = double(0)
                                                        , prob1To2 = double(0, default = -999)
                                                        , prob1To2Hab = double(1)
                                                        , prob2To3 = double(0, default = -999)
                                                        , prob2To3Hab = double(1)
                                                        , prob2To4 = double(0, default = -999)
                                                        , prob2To4Hab = double(1)
                                                        , prob3To4 = double(0, default = -999)
                                                        , prob3To4Hab = double(1)
                                                        , prob2To5 = double(0, default = -999)
                                                        , prob2To5Hab = double(1)
                                                        , prob3To5 = double(0, default = -999)
                                                        , prob3To5Hab = double(1)
                                                        , s = double(1)
                                                        , habitatGrid = double(2)
                                                        , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  # z=1
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
  # z=2
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    #prob2To3
    if(prob2To3== -999){
      indProb2To3 <- prob2To3Hab[sID]
    }else{
      indProb2To3 <- prob2To3
    }
    #prob2To4
    if(prob2To4 == -999){
      indProb2To4 <- prob2To4Hab[sID]
    }else{
      indProb2To4 <- prob2To4
    }
    #prob2To5
    if(prob2To5 == -999){
      indProb2To5 <- prob2To5Hab[sID]
    }else{
      indProb2To5 <- prob2To5
    }
    
    probAlive <- 1-(indProb2To4+indProb2To5)
    
    indProb2To2 <- (1-indProb2To3) * probAlive 
    indProb2To3 <- indProb2To3 * probAlive 
    logLikelihood <- dcat(x, prob = c(0, indProb2To2, indProb2To3, indProb2To4, indProb2To5), log=1)
    
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  # z=3
  if(z == 3){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
   
    #prob3To4
    if(prob3To4 == -999){
      indProb3To4 <- prob3To4Hab[sID]
    }else{
      indProb3To4 <- prob3To4
    }
    #prob3To5
    if(prob3To5 == -999){
      indProb3To5 <- prob3To5Hab[sID]
    }else{
      indProb3To5 <- prob3To5
    }
    
    indProb3To3 <- 1-(indProb3To4+indProb3To5)
    
    logLikelihood <- dcat(x, prob = c(0, 0, indProb3To3, indProb3To4, indProb3To5), log=1)
    
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
  # z=4|5
  if(z == 4 | z == 5 ){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})


NULL
#' @rdname dcatState2Alive2Dead
#' @export
#' 
#### 2.Sampling function ####
rcatState2Alive2Dead <- nimbleFunction(run = function( n = integer(0)
                                                       , z = double(0)
                                                       , prob1To2 = double(0, default = -999)
                                                       , prob1To2Hab = double(1)
                                                       , prob2To3 = double(0, default = -999)
                                                       , prob2To3Hab = double(1)
                                                       , prob2To4 = double(0, default = -999)
                                                       , prob2To4Hab = double(1)
                                                       , prob3To4 = double(0, default = -999)
                                                       , prob3To4Hab = double(1)
                                                       , prob2To5 = double(0, default = -999)
                                                       , prob2To5Hab = double(1)
                                                       , prob3To5 = double(0, default = -999)
                                                       , prob3To5Hab = double(1)
                                                       , s = double(1)
                                                       , habitatGrid = double(2)
){
  # Return type declaration
  returnType(double(0))
  # z=1 
  if(z == 1){
    if(prob1To2 == -999){
      sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
      indProb1To2 <- prob1To2Hab[sID]
    }else{
      ## EXTRACT LOCATION OF THE ID
      indProb1To2 <- prob1To2
    }
    state <- rcat(1, prob = c(1 - indProb1To2, indProb1To2))
    return(state)
  }
  
  # z=2 
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    #prob2To3
    if(prob2To3== -999){
      indProb2To3 <- prob2To3Hab[sID]
    }else{
      indProb2To3 <- prob2To3
    }
    #prob2To4
    if(prob2To4 == -999){
      indProb2To4 <- prob2To4Hab[sID]
    }else{
      indProb2To4 <- prob2To4
    }
    #prob2To5
    if(prob2To5 == -999){
      indProb2To5 <- prob2To5Hab[sID]
    }else{
      indProb2To5 <- prob2To5
    }
    
    probAlive <- 1-(indProb2To4+indProb2To5)
    
    indProb2To2 <- (1-indProb2To3) * probAlive 
    indProb2To3 <- indProb2To3 * probAlive 
    
    
    state <- rcat(1, prob = c(0, indProb2To2, indProb2To3, indProb2To4, indProb2To5))
    return(state)
  }
  
  # z=3 
  if(z == 3){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    
    #prob3To4
    if(prob3To4 == -999){
      indProb3To4 <- prob3To4Hab[sID]
    }else{
      indProb3To4 <- prob3To4
    }
    #prob3To5
    if(prob3To5 == -999){
      indProb3To5 <- prob3To5Hab[sID]
    }else{
      indProb3To5 <- prob3To5
    }
    
    indProb3To3 <- 1-(indProb3To4+indProb3To5)
    
    state <- rcat(1, prob = c(0, 0, indProb3To3, indProb3To4, indProb3To5))
    return(state)
    
  }
  # z=4|5
  if(z == 4 | z == 5 ){
    state <- 5
    return(state)
  }
  
})




