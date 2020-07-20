
#' Zero's trick distribution for irregular habitat shapes
#'
#' The dHabitatMask distribution checks that the proposed activity center location (s) falls
#' within the suitable habitat (defined in the binary matrix habitatMask).
#'
#' The rHabitatMask function returns the value of the habitat mask cell (0 or 1) where the proposed activity center falls. 
#' 
#' @name dHabitatMask
#'
#' @param x Zero's trick data (==1).
#' @param s Bivariate activity center coordinates.
#' @param xmin Minimal value among all trap location x-coordinates
#' @param xmax Maximal value among all trap location x-coordinates
#' @param ymin Minimal value among all trap location y-coordinates
#' @param ymax Maximal value among all trap location y-coordinates
#' @param habitatMask two-dimensional binary habitat mask matrix 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#'
#' @return The log-likelihood value associated with the bivariate activity center location s being in the suitable habitat
#' (i.e. 0 if it falls within the habitat mask and -Inf otherwise)
#'
#' @author Pierre Dupont
#'
#' @import nimble
#'
#' @examples
#' \donttest{
#' 
#' OK <- rep(1, N)
#'
#' habitatMask <- matrix(rbinom(10000,1,0.75), nrow = 100)
#'
#' ## define model code
#' code <- nimbleCode({
#'     for(i in 1:N) {
#'         s[i, 1] ~ dunif(0, 100)
#'         s[i, 2] ~ dunif(0, 100)
#'         
#'         OK[i] ~ dHabitatMask( s = s[i,1:2],
#'                               xmax = 100,
#'                               xmin = 0,
#'                               ymax = 100,
#'                               ymin = 0,
#'                               habitatMask = habMask[1:100,1:100])
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

#' @rdname dHabitatMask
#' @export
dHabitatMask <- nimbleFunction(
  run = function( x = double(0),
                  s = double(1),
                  xmax = double(0),
                  xmin = double(0),
                  ymax = double(0),
                  ymin = double(0),
                  habitatMask = double(2),
                  log = integer(0, default = 0)
  ) {
    returnType(double(0))
    
    if(s[1] < xmin) return(-Inf)   # min x-coordinates
    if(s[1] > xmax) return(-Inf)   # max x-coordinates
    if(s[2] < ymin) return(-Inf)   # min y-coordinates
    if(s[2] > ymax) return(-Inf)   # max y-coordinates
    
    test <- 1-(habitatMask[trunc(s[2])+1, trunc(s[1])+1] == 0)
    
    if(log) return(log(test)) else return(test)
  })

#' @rdname dHabitatMask
#' @export
rHabitatMask <- nimbleFunction(
  run = function( n = integer(),
                  s = double(1),
                  xmax = double(0),
                  xmin = double(0),
                  ymax = double(0),
                  ymin = double(0),
                  habitatMask = double(2)
  ) {
    returnType(double(0))
    
    if(s[1] < xmin) return(0)   # min x-coordinates
    if(s[1] > xmax) return(0)   # max x-coordinates
    if(s[2] < ymin) return(0)   # min y-coordinates
    if(s[2] > ymax) return(0)   # max y-coordinates
    
    if(habitatMask[trunc(s[2])+1, trunc(s[1])+1] == 0) return(0) else return(1)
})


