#' Bivariate exponential dispersal distribution for activity centers
#'
#' This function is deprecated, and it will be removed from a future release.
#' 
#' The dDispersal_exp distribution is a bivariate distribution which can be used to model the latent bivariate activity centers (ACs) of individuals in a population.  This distribution models the situation when individual AC dispersal is uniform in direction (that is, dispersal occurs in a direction theta, where theta is uniformly distributed on [-pi, pi]), and with an exponential distribution for the radial dispersal distance.
#'
#' The dDispersal_exp distribution models the location of an AC at time (t+1), conditional on the previous AC location at time (t) and the rate parameter (rate) of the exponential distribution for dispersal distance.
#'
#' @name dDispersal_exp
#'
#' @param x Bivariate activity center coordinates (at time t+1).
#' @param s Current location of the bivariate activity center (at time t).
#' @param rate Rate parameter of the exponential distribution for dispersal distance.
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#'
#' @return The log-probability value associated with the bivariate activity center location x, given the current activity center s, and the rate parameter of the exponential dispersal distance distribution.
#'
#' @author Daniel Turek
#'
#' @import nimble
#' @importFrom stats dexp rexp runif
#'
#' @examples
#'
#' \dontrun{
#' 
#' ## define model code
#' code <- nimbleCode({
#'     lambda ~ dgamma(0.001, 0.001)
#'     for(i in 1:N) {
#'         AC[i, 1, 1] ~ dunif(0, 100)
#'         AC[i, 2, 1] ~ dunif(0, 100)
#'         for(t in 2:T) {
#'             AC[i, 1:2, t+1] ~ dDispersal_exp(s = AC[i, 1:2, t], rate = lambda)
#'         }
#'     }
#' })
#'
#' constants <- list(N = 10, T = 6)
#'
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code, constants)
#'
#' ## use model object for MCMC, etc.
#'
#' }
#'
#' @export
NULL

#' @rdname dDispersal_exp
#' @export
dDispersal_exp <- nimbleFunction(
    run = function( x = double(1),
                    s = double(1),
                    rate = double(),
                    log = double()) {
        stop('This function is deprecated, and will be removed from a future release.')
        dist <- sqrt(sum((x-s)^2))
        lp <- dexp(dist, rate = rate, log = TRUE) - log(dist)
        returnType(double())
        if(log) return(lp) else return(exp(lp))
    })

#' @rdname dDispersal_exp
#' @export
rDispersal_exp <- nimbleFunction(
    run = function( n = integer(),
                    s = double(1),
                    rate = double()) {
        stop('This function is deprecated, and will be removed from a future release.')
        if(n != 1) stop('rDispersal_exp only implemented for n = 1')
        theta <- runif(1, 0, 2*pi)
        d <- rexp(1, rate = rate)
        newS <- numeric(2)
        newS[1] <- s[1] + d * cos(theta)
        newS[2] <- s[2] + d * sin(theta)
        returnType(double(1))
        return(newS)
    })
