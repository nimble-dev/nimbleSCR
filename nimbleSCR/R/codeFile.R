
#' One phrase description of function.
#' 
#' \code{exampleFunction} returns the sum of the two arguments.  This paragraph, and subsequent paragraphs, contain the longer description of the function.
#' 
#' @param a Description of first argument here
#' @param b Description of second argument here
#' 
#' @return Returns the sum of first two arguments
#' 
#' @examples
#' x <- exampleFunction(1, 2)
#'
#' @export
#' 
exampleFunction <- function(a, b) {
    c <- a + b
    return(c)
}


#' Bivariate dispersal distribution for activity centers
#'
#' The dDisperal distribution is a bivariate distribution which can be used to model the latent bivariate activity centers (ACs) of individuals in a population.  This distribution models the situation when individual AC dispersal is uniform in direction (that is, dispersal occurs in a direction theta, where theta is uniformly distributed on [-pi, pi]), and with an exponential distribution for the radial dispersal distance.
#'
#' The dDispersal distribution models the location of an AC at time (t+1), conditional on the previous AC location at time (t), and also the rate parameter lambda of the exponential distribution for dispersal distance.
#'
#' @name dDispersal
#'
#' @param x Bivariate activity center coordinates (at time t+1).
#' @param S Current location of the bivariate activity center (at time t).
#' @param lam Rate parameter of the exponential distribution for dispersal distance.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#'
#' @return The log-likelihood value associated with the bivariate activity center location x, given the current activity center S, and the lambda parameter of the exponential dispersal distance distribution.
#'
#' @author Daniel Turek
#'
#' @import nimble
#' @importFrom stats dexp rexp runif
#'
#' @examples
#' \donttest{
#' ## define model code
#' code <- nimbleCode({
#'     lambda ~ dgamma(0.001, 0.001)
#'     for(i in 1:N) {
#'         AC[i, 1, 1] ~ dunif(0, 100)
#'         AC[i, 2, 1] ~ dunif(0, 100)
#'         for(t in 2:T) {
#'             AC[i, 1:2, t+1] ~ dDispersal(S = AC[i, 1:2, t], lam = lambda)
#'         }
#'     }
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

#' @rdname dDispersal
#' @export
dDispersal <- nimbleFunction(
    run = function(x = double(1), S = double(1), lam = double(), log = double()) {
        dist <- sqrt(sum((x-S)^2))
        lp <- dexp(dist, rate = lam, log = TRUE) - log(dist)
        returnType(double())
        if(log) return(lp) else return(exp(lp))
    }
)

#' @rdname dDispersal
#' @export
rDispersal <- nimbleFunction(
    run = function(n = integer(), S = double(1), lam = double()) {
        if(n != 1) stop('rDispersal only implemented for n = 1')
        theta <- runif(1, 0, 2*pi)
        d <- rexp(1, rate = lam)
        newS <- numeric(2)
        newS[1] <- S[1] + d * cos(theta)
        newS[2] <- S[2] + d * sin(theta)
        returnType(double(1))
        return(newS)
    }
)




