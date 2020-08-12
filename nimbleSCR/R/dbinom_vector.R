#' Vectorized binomial distribution
#'
#' The \code{dbinom_vector} distribution is a vectorized version of the binomial distribution.
#' It can be used to model a vector of binomial realizations. NB: using the vectorized version 
#' is beneficial only when the entire joint likelihood of the vector of binomial realizations (x)
#' is calculated simultaneously.
#'
#' @name dbinom_vector
#'
#' @param x Vector of quantiles.
#' @param prob Vector of success probabilities on each trial
#' @param size Vector of number of trials (zero or more).
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @param n Number of observations. Only n = 1 is supported.
#'
#' @return The log-likelihood value associated with the vector of binomial observations.
#'
#' @author Pierre Dupont
#'
#' @import nimble
#' @importFrom stats dbinom rbinom
#'
#' @examples
#' ## define vectorized model code
#' code <- nimbleCode({
#'     p ~ dunif(0,1)
#'     p_vector[1:J] <- p
#'     y[1:J] ~ dbinom_vector(size = trials[1:J],
#'                            prob = p_vector[1:J])
#' })
#'  
#' ## simulate binomial data
#' J <- 1000
#' trials <- sample(x = 10, size = J, replace = TRUE)
#' y <- rbinom_vector(J, size = trials, prob = 0.21)
#'  
#' constants <- list(J = J, trials = trials)
#'  
#' data <- list(y = y)
#'  
#' inits <- list(p = 0.5)
#'  
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code, constants, data, inits)
#'  
#' ## use model object for MCMC, etc.
#'
#' @export
NULL

#' @rdname dbinom_vector
#' @export
dbinom_vector <- nimbleFunction(
  run = function( x = double(1),
                  size = double(1),
                  prob = double(1), 
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })

#' @rdname dbinom_vector
#' @export
rbinom_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  size = double(1),
                  prob = double(1)
  ) {
    returnType(double(1))
    return(rbinom(length(size), prob = prob, size = size))
  })



