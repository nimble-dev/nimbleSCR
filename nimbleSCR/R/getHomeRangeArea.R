###############################################################
################ HOME RANGE RADIUS AND AREA ###################
###############################################################
#' @title Computation of home range radius and area
#'
#' @description
#' \code{getHomeRangeArea} returns the estimates of home range radius and area for a given set of parameters with respect to a specified detection function using bisection algorithm. The following circular detection functions are available to use in nimbleSCR: half-normal (detFun = 0), half-normal plateau (detFun = 1), exponential (detFun = 2), asymmetric logistic (detFun = 3), bimodal (detFun = 4) and donut (detFun = 5).
#' 
#' @name getHomeRangeArea
#' 
#' @param x \code{Vector} or \code{matrix} (parameters in columns) of values for different parameters corresponding to the specified detection function.
#' @param detFun \code{Numeric} variable denoting the type of detection function. 0 = Half-normal, 1 = Half-normal plateau, 2 = Exponential, 3 = Asymmetric logistic, 4 = Bimodal, 5 = Donut.
#' @param prob \code{Numeric}  variable denoting the quantile probability to compute the home range radius.
#' @param d \code{Numeric} variable giving an initial value of the radius.
#' @param xlim \code{Vector} of length 2 giving the range along x-axis.
#' @param ylim \code{Vector} of length 2 giving the range along y-axis.
#' @param nBreaks \code{Numeric} variable denoting the number of breaks along an axis.
#' @param tol \code{Numeric} variable denoting the allowed tolerance in the radius estimate.
#' @param nIter \code{Numeric} variable giving the maximum number of iterations in bisection algorithm.
#'
#' @author Soumen Dey
#' 
#' @references
#' 
#' Dey, S., Bischof, R., Dupont, P. P. A., & Milleret, C. (2022). Does the punishment fit the crime? Consequences and diagnosis of misspecified detection functions in Bayesian spatial capture–recapture modeling. Ecology and Evolution, 12, e8600. https://doi.org/10.1002/ece3.8600
#' 
#' @import nimble
#' 
#' 
#' @examples
#'
#' \dontrun{
#' 
#' # A user friendly vignette is also available on github: 
#' # https://github.com/nimble-dev/nimbleSCR/blob/master/nimbleSCR/vignettes/
#' # Vignette name: Fit_with_dbinomLocal_normalPlateau_and_HomeRangeAreaComputation.rmd
#' 
#' # HALF-NORMAL PLATEAU FUNCTION (detFun = 1) 
#' habitatMask <- matrix(1, nrow = 30, ncol= 30, byrow = TRUE)
#' 
#' prob <- 0.95
#' paramnames.hr <- c("HRradius", "HRarea")
#' sigma <- 1
#' w <- 1.5
#' params <- c(sigma, w)
#' names(params) <- c("sigma", "w")
#' HRAnim <- getHomeRangeArea( x = params, detFun = 1, prob = prob, d = 6,
#'                       xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                       nBreaks = 800, tol = 1E-5, nIter = 2000)
#' 
#' # Different values of argument "detFun"
#' # 0 = Half-normal, 1 = Half-normal plateau, 2 = Exponential,
#' # 3 = Aysmmetric logistic, 4 = Bimodal, 5 = Donut.
#' HR.hnp <- c(HRAnim$run())
#' names(HR.hnp) <- paramnames.hr
#' print(HR.hnp)
#' # FASTER HRA COMPUTATION USING NIMBLE
#' samples <- cbind(rgamma(n = 500, shape = 1, rate = 1), rgamma(n = 500, shape = 1.5, rate = 1))
#' colnames(samples) <- c("sigma", "w")
#' HRAnim.mat <- getHomeRangeArea(x = samples, detFun = 1, prob = prob, d = 6, 
#'                          xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                          nBreaks = 800, tol = 1E-5, nIter = 2000)
#' 
#' cHRAnim.arr <- compileNimble(HRAnim.mat, resetFunctions = TRUE)
#' 
#' HRA.Runtime <- system.time(
#'   HR.chain <- cHRAnim.arr$run()
#' )
#' print(HRA.Runtime)
#' dimnames(HR.chain)[[2]] <- paramnames.hr
#' HRest <- do.call(rbind, lapply(c(1:2), function(j){
#'   c(mean(HR.chain[,j], na.rm = TRUE), sd(HR.chain[,j], na.rm = TRUE))
#' }))
#' dimnames(HRest) <- list(paramnames.hr, c("MEAN", "SD"))
#' 
#' cat("Numerical estimates using MCMC samples: \n", sep = "")
#' print(HRest)
#' 
#' # HALF-NORMAL FUNCTION (detFun = 0) 
#' sigma = 2
#' params <- c(sigma)
#' names(params) <- c("sigma")
#' 
#' HRAnim <- getHomeRangeArea(x = params, detFun = 0, prob = prob, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      nBreaks = 800, tol = 1E-5, nIter = 2000)
#' 
#' HR.hn <- c(HRAnim$run())
#' names(HR.hn) <- paramnames.hr
#' print(HR.hn)
#' 
#' # Exponential (detFun = 2) 
#' 
#' rate = 1/2
#' params <- c(rate)
#' names(params) <- c("rate")
#' HRAnim <- getHomeRangeArea(x = params, detFun = 2, prob = prob, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      nBreaks = 800, tol = 1E-5, nIter = 2000)
#' HR.exp <- c(HRAnim$run())
#' names(HR.exp) <- paramnames.hr
#' print(HR.exp)
#' 
#' # Asymmetric logistic (detFun = 3) 
#' 
#' sigma = 2
#' alpha.a = 5 
#' alpha.b = 1
#' params <- c(sigma, alpha.a, alpha.b)
#' names(params) <- c("sigma", "alpha.a", "alpha.b")
#' HRAnim <- getHomeRangeArea(x = params, detFun = 3, prob = prob, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      nBreaks = 800, tol = 1E-5, nIter = 2000)
#' HR.al <- c(HRAnim$run())
#' names(HR.al) <- paramnames.hr
#' print(HR.al)
#' 
#' 
#' # Bimodal (detFun = 4) 
#' 
#' p0.a = 0.25
#' sigma.a = 0.5
#' p0.b = 0.15 
#' sigma.b = 1
#' w = 2
#' params <- c(sigma.a, sigma.b, p0.a, p0.b, w)
#' names(params) <- c("sigma.a", "sigma.b", "p0.a", "p0.b", "w")
#' HRAnim <- getHomeRangeArea(x = params, detFun = 4, prob = prob, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      nBreaks = 800, tol = 1E-5, nIter = 2000)
#' HR.bi <- c(HRAnim$run())
#' names(HR.bi) <- paramnames.hr
#' print(HR.bi)
#' 
#' # Donut (detFun = 5) 
#' 
#' sigma.a = 1.5 
#' sigma.b = 1 
#' w = 1 
#' params <- c(sigma.a, sigma.b, w)
#' names(params) <- c("sigma.a", "sigma.b", "w")
#' HRAnim <- getHomeRangeArea(x = params, detFun = 5, prob = prob, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      nBreaks = 800, tol = 1E-5, nIter = 2000)
#' HR.dn <- c(HRAnim$run())
#' names(HR.dn) <- paramnames.hr
#' print(HR.dn)
#'
#' }
#' 
#' @export
NULL

#' @rdname getHomeRangeArea
#' @export
#' 
#Bisection algorithm
getHomeRangeArea <- nimbleFunction(
  setup = function(x=c(2), detFun=0, prob = 0.95, d = 6,
                   xlim=c(0,30), ylim=c(0,30), 
                   nBreaks=800, tol=1E-5, nIter=2000){
    s <- c(sum(xlim)/2, sum(ylim)/2) # center of the rectangular grid
    x1 <- rep(seq(xlim[1], xlim[2], length.out = nBreaks), nBreaks) # vector (nBreaks^2x 1)
    x2.temp <- seq(ylim[1], ylim[2], length.out = nBreaks)
    x2 <- nimNumeric(value = 0, length = nBreaks^2, init = TRUE) 
    for(i in 1:nBreaks){
      x2[((i-1)*nBreaks+1):(i*nBreaks)] <- rep(x2.temp[i], nBreaks) # NIMBLE
    }
    n <- length(x1) # nBreaks^2
    D <- sqrt((s[1] - x1[1:n])^2 + (s[2] - x2[1:n])^2 ) # vector (n x 1) # NIMBLE
    resolution <- min(diff(x1[1:10])) # Resolution of the grid
    
    if(!is.matrix(x)){ x <- t(as.matrix(x))}
    
  }, #setup
  run = function(){
    
    nrows <- dim(x)[1]
    # nrows <- length(p0Vec)
    out <- matrix(NA, nrows, 2)
    
    for(this.row in 1:nrows ){
      if (this.row%%1000 == 0)  cat('..... row #', this.row, '\n')
      
      if(detFun == 0){ # HALF-NORMAL
        # paramnames = c('sigma') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- x[this.row, 1] # Scale parameter 
        p <- exp(-D[1:n]*D[1:n]/(2*sigma*sigma)) # vector (n x 1)
      }
      if(detFun == 1){ #  HALF-NORMAL PLATEAU
        # paramnames = c('sigma', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- x[this.row, 1] # Scale parameter 
        w <- x[this.row, 2] # Length of plateau
        p <- nimNumeric(length = n, value = 1, init = TRUE)
        for(i in 1:n){
          if(D[i] > w) {p[i] <- exp(-(D[i] - w)*(D[i] - w)/(2*sigma*sigma))}
        }
      }
      if(detFun == 2){ # EXPONENTIAL
        # paramnames = c('rate') # PARAMETER NAMES IN THEIR ORDER OF USE
        rate <- x[this.row, 1] # Rate parameter 
        p <- exp(-rate * D[1:n]) # vector (n x 1)
      }
      if(detFun == 3){ # ASYMMETRIC LOGISTIC
        # paramnames = c('sigma', 'alpha.a', 'alpha.b') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- x[this.row, 1] # Scale parameter 
        alpha.a <- x[this.row, 2]
        alpha.b <- x[this.row, 3]
        Anuslope <- (2*abs(alpha.a)*alpha.b)/(1+alpha.b)
        fx <- 1/ (1+(D[1:n]/sigma)^Anuslope)
        dens <- 1+fx[1:n]*((D[1:n]/sigma)^(alpha.a))+(1-fx[1:n])*((D[1:n]/sigma)^(alpha.a*alpha.b))
        p <- 1/dens[1:n]
      }
      if(detFun == 4){ # BIMODAL
        # paramnames = c('sigma', 'sigma.b', 'p0', 'p0.b', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma.a <- x[this.row, 1] # Scale parameter of the first peak 
        sigma.b <- x[this.row, 2] # Scale parameter of the second peak 
        p0.a <- x[this.row, 3] # Baseline detection probability of the first peak 
        p0.b <- x[this.row, 4] # Baseline detection probability of the second peak 
        w <- x[this.row, 5] # Distance between the home-range center and the "doughnut" center
        densa <- p0.a *  exp(-D[1:n]*D[1:n]/(2.0*sigma.a*sigma.a)) 
        densb <-  p0.b * exp(-(D[1:n]-w)*(D[1:n]-w)/(2.0*sigma.b*sigma.b))
        p <- densa[1:n] + densb[1:n]
        
      }
      if(detFun == 5){ # DONUT
        # paramnames = c('sigma', 'sigma.b', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma.a <- x[this.row, 1] # Scale parameter of the left tail
        sigma.b <- x[this.row, 2] # Scale parameter of the right tail
        w <- x[this.row, 3] # Distance between the activity centre and the "donut" centre
        p <- nimNumeric(length = n, value = 0, init = TRUE)
        for(i in 1:n){
          if(D[i] <= w) {p[i] <- exp(-(D[i] - w)*(D[i] - w)/(2*sigma.a*sigma.a))}
          if(D[i] > w) {p[i] <- exp(-(D[i] - w)*(D[i] - w)/(2*sigma.b*sigma.b))}
        }
      }#if
      
      psi <- p[1:n]/sum(p[1:n]) # vector (n x 1)
      
      d.lower <- min(D)
      d.upper <- max(D)
      
      temp <- 0
      for(i in 1:n){
        if(D[i] <= d) temp <- temp + psi[i]
      }
      iter <- 1 
      while(iter <= nIter &
            abs(d.upper-d.lower) > d*tol &
            abs(d.upper-d.lower) > resolution){
        
        if(temp >= prob){
          d <<- (d.lower+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < prob) d.lower <- d
          if(temp >= prob) d.upper <- d ## NEW
        }
        if(temp < prob){
          d <<- (d.upper+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < prob) d.lower <- d ## NEW
          if(temp >= prob) d.upper <- d
        }
        iter <- iter + 1
        
      }#while
      if(temp < prob) {d <<- d*(1+tol) }
      radius <- d
      area <- pi * radius^2
      out[this.row, 1:2] <- c(radius, area)
      
    }#this.row
    
    returnType(double(2))
    
    return(out) 
  }) #run
