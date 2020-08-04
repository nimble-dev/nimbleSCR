#' Vectorized binomial distribution
#'
#' The \code{dbinom_vector} distribution is a vectorized version of the binomial distribution.
#' It can be used to model a vector of binomial realizations. NB: using the vectorized version 
#' is beneficial only when the entire joint likelihood of the vector of binomial realizations (x)
#' can be calculated simultaneously.
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
#' \donttest{
#' # simulate binomial data
#' p <- rep(0.21, 1000)
#' trials <- sample(x = 10, size = 1000, replace = T)
#' y <- rbinom_vector(n = 1000, size = trials, prob = p)
#' 
#' ## define model code
#' code <- nimbleCode({
#'     p ~ dunif(0,1)
#'     for(j in 1:J) {
#'       p_vector[i] <- p
#'       y[j] ~ dbinom( size = trials[i],
#'                      prob = p_vector[i])
#'     }
#' })
#'
#'## define vectorized model code
#' code <- nimbleCode({
#'     p ~ dunif(0,1)
#'     p_vector[1:J] <- p
#'     y[1:J] ~ dbinom_vector( size = trials[1:J],
#'                             prob = p_vector[1:J])
#' })
#'
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code,data=list(y=y,trials=trials),constants = list(J=length(y)))
#'
#' ## use model object for MCMC, etc.
#' }
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



