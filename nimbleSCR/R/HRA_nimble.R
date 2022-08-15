###############################################################
################ HOME RANGE RADIUS AND AREA ###################
###############################################################
#' @title Computation of home range radius and area
#'
#' @description
#' \code{HRA_nimble} returns the estimates of home range radius and area for a given set of parameters with respect to a specified detection function using bisection algorithm. The following circular detection functions are available to use in nimbleSCR: half-normal(HN, detfun = 0), half-normal plateau (HNP, detfun = 1), exponential (EX, detfun = 2), asymmetric logistic (AL, detfun = 3), bimodal (BI, detfun = 4) and donut (DN, detfun = 5).
#' 
#' @name HRA_nimble
#' 
#' @param param \code{Vector} or \code{matrix} (parameters in columns) of values for different parameters corresponding to the specified detection function.
#' @param detfun \code{Numeric} variable denoting the type of detection function. 0 = Half-normal (HN), 1 = Half-normal plateau (HNP), 2 = Exponential (EX), 3 = Asymmetric logistic (AL), 4 = Bimodal (BI), 5 = Donut (DN).
#' @param pr \code{Numeric}  variable denoting the quantile probability to compute the home range radius.
#' @param d \code{Numeric} variable giving an initial value of the radius.
#' @param xlim \code{Vector} of length 2 giving the range along x-axis.
#' @param ylim \code{Vector} of length 2 giving the range along y-axis.
#' @param ng \code{Numeric} variable denoting the number of breaks along an axis.
#' @param tol \code{Numeric} variable denoting the allowed tolerance in the radius estimate.
#' @param niter \code{Numeric} variable giving the maximum number of iterations in bisection algorithm.
#'
#' @author Soumen Dey
#' 
#' @import nimble
#' 
#' 
#' @examples
#' # A user friendly vignette is also available on github: 
#' # https://github.com/nimble-dev/nimbleSCR/blob/master/nimbleSCR/vignettes/
#' # Vignette name: Fit_SCR_models_with_dbinomLocal_HNP_and_HomeRangeRadiusComputation.rmd
#' 
#' # HALF-NORMAL PLATEAU FUNCTION (HNP, detfun = 1) 
#' habitatMask <- matrix(1, nrow = 30, ncol= 30, byrow = TRUE)
#' 
#' pr <- 0.95
#' paramnames.hr <- c("HRradius", "HRarea")
#' sigma <- 1
#' w <- 1.5
#' params <- c(sigma, w)
#' names(params) <- c("sigma", "w")
#' HRAnim <- HRA_nimble( param = params, detfun = 1, pr = pr, d = 6,
#'                       xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                       ng = 800, tol = 1E-5, niter = 2000)
#' 
#' # Different values of argument "detfun"
#' # 0 = Half-normal (HN), 1 = Half-normal plateau (HNP), 2 = Exponential (EX),
#' # 3 = Aysmmetric logistic (AL), 4 = Bimodal (BI), 5 = Donut (DN).
#' HR.hnp <- c(HRAnim$run())
#' names(HR.hnp) <- paramnames.hr
#' print(HR.hnp)
#' # FASTER HRA COMPUTATION USING NIMBLE
#' samples <- cbind(rgamma(n = 500, shape = 1, rate = 1), rgamma(n = 500, shape = 1.5, rate = 1))
#' colnames(samples) <- c("sigma", "w")
#' HRAnim.mat <- HRA_nimble(param = samples, detfun = 1, pr = pr, d = 6, 
#'                          xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                          ng = 800, tol = 1E-5, niter = 2000)
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
#' # HALF-NORMAL FUNCTION (HN, detfun = 0) 
#' sigma = 2
#' params <- c(sigma)
#' names(params) <- c("sigma")
#' 
#' HRAnim <- HRA_nimble(param = params, detfun = 0, pr = pr, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      ng = 800, tol = 1E-5, niter = 2000)
#' 
#' HR.hn <- c(HRAnim$run())
#' names(HR.hn) <- paramnames.hr
#' print(HR.hn)
#' 
#' # Exponential (EX, detfun = 2) 
#' 
#' sigma = 2
#' params <- c(sigma)
#' names(params) <- c("sigma")
#' HRAnim <- HRA_nimble(param = params, detfun = 2, pr = pr, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      ng = 800, tol = 1E-5, niter = 2000)
#' HR.ex <- c(HRAnim$run())
#' names(HR.ex) <- paramnames.hr
#' print(HR.ex)
#' 
#' # Asymmetric logistic (AL, detfun = 3) 
#' 
#' sigma = 2
#' alpha.a = 5 
#' alpha.b = 1
#' params <- c(sigma, alpha.a, alpha.b)
#' names(params) <- c("sigma", "alpha.a", "alpha.b")
#' HRAnim <- HRA_nimble(param = params, detfun = 3, pr = pr, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      ng = 800, tol = 1E-5, niter = 2000)
#' HR.al <- c(HRAnim$run())
#' names(HR.al) <- paramnames.hr
#' print(HR.al)
#' 
#' 
#' # Bimodal (BI, detfun = 4) 
#' 
#' p0.a = 0.25
#' sigma.a = 0.5
#' p0.b = 0.15 
#' sigma.b = 1
#' w = 2
#' params <- c(sigma.a, sigma.b, p0.a, p0.b, w)
#' names(params) <- c("sigma.a", "sigma.b", "p0.a", "p0.b", "w")
#' HRAnim <- HRA_nimble(param = params, detfun = 4, pr = pr, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      ng = 800, tol = 1E-5, niter = 2000)
#' HR.bi <- c(HRAnim$run())
#' names(HR.bi) <- paramnames.hr
#' print(HR.bi)
#' 
#' # Donut (DN, detfun = 5) 
#' 
#' sigma.a = 1.5 
#' sigma.b = 1 
#' w = 1 
#' params <- c(sigma.a, sigma.b, w)
#' names(params) <- c("sigma.a", "sigma.b", "w")
#' HRAnim <- HRA_nimble(param = params, detfun = 5, pr = pr, d = 6, 
#'                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#'                      ng = 800, tol = 1E-5, niter = 2000)
#' HR.dn <- c(HRAnim$run())
#' names(HR.dn) <- paramnames.hr
#' print(HR.dn)
#' 
#' @export
NULL

#' @rdname HRA_nimble
#' @export
#' 
#Bisection algorithm
HRA_nimble <- nimbleFunction(
  setup = function(param=c(2), detfun=0, pr = 0.95, d = 6,
                   xlim=c(0,30), ylim=c(0,30), 
                   ng=800, tol=1E-5, niter=2000){
    s <- c(sum(xlim)/2, sum(ylim)/2) # center of the rectangular grid
    x1 <- rep(seq(xlim[1], xlim[2], length.out = ng), ng) # vector (ng^2x 1)
    x2.temp <- seq(ylim[1], ylim[2], length.out = ng)
    x2 <- nimNumeric(value = 0, length = ng^2, init = TRUE) 
    for(i in 1:ng){
      x2[((i-1)*ng+1):(i*ng)] <- rep(x2.temp[i], ng) # NIMBLE
    }
    n <- length(x1) # ng^2
    D <- sqrt((s[1] - x1[1:n])^2 + (s[2] - x2[1:n])^2 ) # vector (n x 1) # NIMBLE
    resolution <- min(diff(x1[1:10])) # Resolution of the grid
    
    if(!is.matrix(param)){ param <- t(as.matrix(param))}
    
  }, #setup
  run = function(){
    
    nrows <- dim(param)[1]
    # nrows <- length(p0Vec)
    out <- matrix(NA, nrows, 2)
    
    for(this.row in 1:nrows ){
      if (this.row%%1000 == 0)  cat('..... row #', this.row, '\n')
      
      if(detfun == 0){ #'HN'
        # paramnames = c('sigma') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- param[this.row, 1] # Scale parameter 
        p <- exp(-D[1:n]*D[1:n]/(2*sigma*sigma)) # vector (n x 1)
      }
      if(detfun == 1){ #  'HNP'
        # paramnames = c('sigma', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- param[this.row, 1] # Scale parameter 
        w <- param[this.row, 2] # Length of plateau
        p <- nimNumeric(length = n, value = 1, init = TRUE)
        for(i in 1:n){
          if(D[i] > w) {p[i] <- exp(-(D[i] - w)*(D[i] - w)/(2*sigma*sigma))}
        }
      }
      if(detfun == 2){ # 'EX'
        # paramnames = c('sigma') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- param[this.row, 1] # Scale parameter 
        p <- exp(-D[1:n]/sigma) # vector (n x 1)
      }
      if(detfun == 3){ #'AL'
        # paramnames = c('sigma', 'alpha.a', 'alpha.b') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma <- param[this.row, 1] # Scale parameter 
        alpha.a <- param[this.row, 2]
        alpha.b <- param[this.row, 3]
        Anuslope <- (2*abs(alpha.a)*alpha.b)/(1+alpha.b)
        fx <- 1/ (1+(D[1:n]/sigma)^Anuslope)
        dens <- 1+fx[1:n]*((D[1:n]/sigma)^(alpha.a))+(1-fx[1:n])*((D[1:n]/sigma)^(alpha.a*alpha.b))
        p <- 1/dens[1:n]
      }
      if(detfun == 4){ #'BI'
        # paramnames = c('sigma', 'sigma.b', 'p0', 'p0.b', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma.a <- param[this.row, 1] # Scale parameter of the first peak 
        sigma.b <- param[this.row, 2] # Scale parameter of the second peak 
        p0.a <- param[this.row, 3] # Baseline detection probability of the first peak 
        p0.b <- param[this.row, 4] # Baseline detection probability of the second peak 
        w <- param[this.row, 5] # Distance between the home-range center and the "doughnut" center
        densa <- p0.a *  exp(-D[1:n]*D[1:n]/(2.0*sigma.a*sigma.a)) 
        densb <-  p0.b * exp(-(D[1:n]-w)*(D[1:n]-w)/(2.0*sigma.b*sigma.b))
        p <- densa[1:n] + densb[1:n]
        
      }
      if(detfun == 5){ #'DN'
        # paramnames = c('sigma', 'sigma.b', 'w') # PARAMETER NAMES IN THEIR ORDER OF USE
        sigma.a <- param[this.row, 1] # Scale parameter of the left tail
        sigma.b <- param[this.row, 2] # Scale parameter of the right tail
        w <- param[this.row, 3] # Distance between the activity centre and the "donut" centre
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
      while(iter <= niter &
            abs(d.upper-d.lower) > d*tol &
            abs(d.upper-d.lower) > resolution){
        
        if(temp >= pr){
          d <<- (d.lower+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < pr) d.lower <- d
          if(temp >= pr) d.upper <- d ## NEW
        }
        if(temp < pr){
          d <<- (d.upper+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < pr) d.lower <- d ## NEW
          if(temp >= pr) d.upper <- d
        }
        iter <- iter + 1
        
      }#while
      if(temp < pr) {d <<- d*(1+tol) }
      radius <- d
      area <- pi * radius^2
      out[this.row, 1:2] <- c(radius, area)
      
    }#this.row
    
    returnType(double(2))
    
    return(out) 
  }) #run
