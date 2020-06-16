
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
#' Longer description here.  Describe how to use it, etc.
#'
#' @name dDispersal
#'
#' @param x DESCRIPTION OF x PARAMETER
#' @param S DESCRIPTION OF S PARAMETER
#' @param lam DESCRIPTION OF lam PARAMETER
#' @param log DESCRIPTION OF log PARAMETER
#' @param n DESCRIPTION OF n PARAMETER
#'
#' @return returns the log-likelihood value associated with the bivariate activity center location x, given the current activity center S, and the lambda parameter of the exponential dispersal distance distribution
#'
#' @author Daniel Turek
#'
#' @import nimble
#' @importFrom stats dexp rexp runif
#'
#' @examples
#' ## R CODE EXAMPLE HERE USING dDispersal DISTRIBUTION
#' x[1:2] ~ dDispersal()
#' ### etc ...
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




