
nimbleOptions(verbose = FALSE)
# tempFileName <- 'distributionsTestLog.Rout'
# generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForDistributionsTesting'))
# outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForDistributionsTesting'), goldFileName) else tempFileName
# 
# sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)



#---- 1. SIMULATE AND FIT MODELS ----
#----   1.1 dpoisLocal_normal() WITH sxy~dunif() ----
test_that("Simulate and fit dpoisLocal_normal ",{
  
  coordsHabitatGridCenter <- cbind(rep(seq(11.5, 0.5, by=-1), 10),
                                   sort(rep(seq(0.5, 11.5, by=1), 10)))
  colnames(coordsHabitatGridCenter) <- c("x", "y")
  
  # CREATE A 8*8 TRAP GRID CENTERED ON THE HABITAT GRID
  trapCoords <- cbind(rep(seq(2.5, 9.5,by=1),8),
                      sort(rep(seq(2.5, 9.5,by=1),8)))
  colnames(trapCoords) <- c("x","y")
  
  # PLOT CHECK
  # plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch=16)
  # points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16)
  
  ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                               coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
  habitatMask <- matrix(1, nrow = 12, ncol= 12, byrow = TRUE)
  
  # CREATE LOCAL OBJECTS TO APPLY LOCAL EVALUATION
  TrapLocal <- getLocalObjects(habitatMask = habitatMask,
                               coords = ScaledtrapCoords,
                               dmax = 7,
                               resizeFactor = 1,
                               plot.check = F
  )
  modelCode <- nimbleCode({
    ## priors
    psi ~ dunif(0, 1)
    sigma ~ dunif(0, 50)
    lambda ~ dunif(0, 10)
    ## loop over individuals
    for(i in 1:n.individuals) {
      ## AC coordinates
      sxy[i,1] ~ dunif(0, x.max)
      sxy[i,2] ~ dunif(0, y.max)
      ## habitat constraint
      ones[i] ~ dHabitatMask( s = sxy[i,1:2],
                              xmin = lowerCoords[1],
                              xmax = upperCoords[1],
                              ymin = lowerCoords[2],
                              ymax = upperCoords[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
      ## latent dead/alive indicators
      z[i] ~ dbern(psi)
      ## likelihood
      y[i, 1:lengthYCombined] ~ dpoisLocal_normal(# detNums = nbDetections[i],
        # detIndices = yDets[i,1:nMaxDetectors],
        lambda = lambda,
        s = sxy[i,1:2],
        sigma = sigma,
        trapCoords = trapCoords[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
        localTrapsNum = nDetectors[1:n.cells],
        resizeFactor = ResizeFactor,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = z[i],
        lengthYCombined = lengthYCombined)
    }
    ## derived quantity: total population size
    N <- sum(z[1:n.individuals])
  })
  
  lambda <- 0.2
  sigma <- 2
  psi <- 0.5
  n.individuals <- 200
  lengthYCombined <- TrapLocal$numLocalIndicesMax*2+1
  
  nimConstants <- list(n.individuals = n.individuals,
                       n.detectors = dim(ScaledtrapCoords)[1],
                       y.max = dim(habitatMask)[1],
                       x.max = dim(habitatMask)[2],
                       y.maxDet = dim(TrapLocal$habitatGrid)[1],
                       x.maxDet = dim(TrapLocal$habitatGrid)[2],
                       ResizeFactor = TrapLocal$resizeFactor,
                       n.cells = dim(TrapLocal$localIndices)[1],
                       maxNBDets = TrapLocal$numLocalIndicesMax,
                       detectorIndex = TrapLocal$localIndices,
                       nDetectors = TrapLocal$numLocalIndices,
                       habitatIDDet = TrapLocal$habitatGrid,
                       lengthYCombined = lengthYCombined)
  
  
  nimData <- list(trapCoords = ScaledtrapCoords,
                  habitat.mx = habitatMask,
                  ones = rep(1, nimConstants$n.individuals),
                  lowerCoords = c(min(coordsHabitatGridCenter[,1]) - 0.5, min(coordsHabitatGridCenter[,2]) - 0.5),
                  upperCoords = c(max(coordsHabitatGridCenter[,1]) + 0.5, max(coordsHabitatGridCenter[,2]) + 0.5))
                  

  # We set the parameter values as inits
  nimInits <- list(lambda = lambda,
                   psi = psi,
                   sigma = sigma)
  
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData,
                        inits = nimInits,
                        check = F,
                        calculate = F)
  # FIRST WE GET THE NODES TO SIMULATE
  nodesToSim <- model$getDependencies(c("sxy", "z"), self=T)
  # THEN WE SIMULATE THOSE NODES
  set.seed(100)
  model$simulate(nodesToSim, includeData = FALSE)
  N <- sum(model$z)
  nimData1 <- nimData
  nimData1$y <- model$y
  nimInits1 <- nimInits
  nimInits1$z <- model$z
  nimInits1$sxy <- model$sxy
  
  
  # CREATE AND COMPILE THE NIMBLE MODEL
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData1,
                        inits = nimInits1,
                        check = F,
                        calculate = F)
  logProb <- model$calculate()#"-2009.3471687696"
  
  ## test R logProp
  expect_match(sprintf("%.10f", logProb),## USE THIS TO THE GET DIGITS VALUES
               "-2009.3471687696",
               info = paste0("Incorrect dpoisLocal_normal LogProb in R"))
  
  cmodel <- compileNimble(model)
  ClogProb <- cmodel$calculate()#-1963.177
  ## test c++ logProp
  expect_match(sprintf("%.10f", ClogProb),
               "-2009.3471687696",
               info = paste0("Incorrect dpoisLocal_normal LogProb in C++"))
  
  
  MCMCconf <- configureMCMC(model = model,
                            monitors = c("N", "sigma", "lambda","psi"),
                            control = list(reflective = TRUE),
                            thin = 1)
  
  ## Add block sampling of sxy coordinates
  MCMCconf$removeSamplers("sxy")
  ACnodes <- paste0("sxy[", 1:nimConstants$n.individuals, ", 1:2]")
  for(node in ACnodes) {
    MCMCconf$addSampler(target = node,
                        type = "RW_block",
                        control = list(adaptScaleOnly = TRUE),
                        silent = TRUE)
  }
  
  MCMC <- buildMCMC(MCMCconf)
  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
  # RUN THE MCMC
  set.seed(100)
  myNimbleOutput <- runMCMC( mcmc = cMCMC,
                             nburnin = 1000,
                             niter = 5000,
                             nchains = 3,
                             samplesAsCodaMCMC = TRUE)
  sum <- summary(myNimbleOutput)
  ## COMPARE N
  expect_equal(sum$statistics["N","Mean"] ,
               N,
               tolerance=1,
               info = paste0("Non-comparable mean estimated N from simulation with dpoisLocal_normal"))
  
  #COMPARE "p0","psi","sigma"
  True <- sum$statistics[c("lambda","psi","sigma"),"Mean"]
  True[] <- c(lambda, psi, sigma)
  
  expect_equal(as.numeric(sum$statistics[c("lambda","psi","sigma"),"Mean"]) ,
               as.numeric(True),
               tolerance = 0.07,
               info = paste0("Non-comparable mean estimated lambda, psi, sigma from simulation with dpoisLocal_normal"))
})



#----   2.1 dpoisLocal_normal() WITH sxy~dunif() and trap covariates  ----
test_that("Simulate and fit dpoisLocal_normal with trap covariates",{
  
  coordsHabitatGridCenter <- cbind(rep(seq(11.5, 0.5, by=-1), 10),
                                   sort(rep(seq(0.5, 11.5, by=1), 10)))
  colnames(coordsHabitatGridCenter) <- c("x","y")
  
  # CREATE A 8*8 TRAP GRID CENTERED ON THE HABITAT GRID
  trapCoords <- cbind(rep(seq(2.5, 9.5,by=1),8),
                      sort(rep(seq(2.5, 9.5,by=1),8)))
  colnames(trapCoords) <- c("x","y")
  
  # PLOT CHECK
  # plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch=16)
  # points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16)
  
  ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                               coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
  habitatMask <- matrix(1, nrow = 12, ncol= 12, byrow = TRUE)
  
  # CREATE LOCAL OBJECTS TO APPLY LOCAL EVALUATION
  TrapLocal <- getLocalObjects(habitatMask = habitatMask,
                               coords = ScaledtrapCoords,
                               dmax = 7,
                               resizeFactor = 1,
                               plot.check = F
  )
  
  set.seed(50)
  trapCov <- cbind(runif(dim(ScaledtrapCoords)[1],-1, 1),
                   runif(dim(ScaledtrapCoords)[1],-1, 1))
  
  
  
  modelCode <- nimbleCode({
    ## priors
    psi ~ dunif(0, 1)
    sigma ~ dunif(0, 50)
    lambda ~ dunif(0, 10)
    
    
    betaTrap[1] ~ dunif(-5,5)
    betaTrap[2] ~ dunif(-5,5)
    
    
    lambdaTraps[1:n.detectors] <- exp(lambda + betaTrap[1] * trapCov[1:n.detectors, 1] +
                                       betaTrap[2] * trapCov[1:n.detectors, 2])
    
    
    
    
    ## loop over individuals
    for(i in 1:n.individuals) {
      ## AC coordinates
      sxy[i,1] ~ dunif(0, x.max)
      sxy[i,2] ~ dunif(0, y.max)
      ## habitat constraint
      ones[i] ~ dHabitatMask( s = sxy[i,1:2],
                              xmin = lowerCoords[1],
                              xmax = upperCoords[1],
                              ymin = lowerCoords[2],
                              ymax = upperCoords[2],
                              habitat = habitat.mx[1:y.max,1:x.max])
      ## latent dead/alive indicators
      z[i] ~ dbern(psi)
      ## likelihood
      y[i, 1:lengthYCombined] ~ dpoisLocal_normal(# detNums = nbDetections[i],
        # detIndices = yDets[i,1:nMaxDetectors],
        lambdaTraps = lambdaTraps[1:n.detectors],
        s = sxy[i,1:2],
        sigma = sigma,
        trapCoords = trapCoords[1:n.detectors,1:2],
        localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
        localTrapsNum = nDetectors[1:n.cells],
        resizeFactor = ResizeFactor,
        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
        indicator = z[i],
        lengthYCombined = lengthYCombined)
    }
    ## derived quantity: total population size
    N <- sum(z[1:n.individuals])
  })
  
  lambda <- 0.2
  sigma <- 2
  psi <- 0.5
  betaTrap <- c(-2, 2)
  n.individuals <- 200
  lengthYCombined <- TrapLocal$numLocalIndicesMax*2+1
  
  nimConstants <- list(n.individuals = n.individuals,
                       n.detectors = dim(ScaledtrapCoords)[1],
                       y.max = dim(habitatMask)[1],
                       x.max = dim(habitatMask)[2],
                       y.maxDet = dim(TrapLocal$habitatGrid)[1],
                       x.maxDet = dim(TrapLocal$habitatGrid)[2],
                       ResizeFactor = TrapLocal$resizeFactor,
                       n.cells = dim(TrapLocal$localIndices)[1],
                       maxNBDets = TrapLocal$numLocalIndicesMax,
                       detectorIndex = TrapLocal$localIndices,
                       nDetectors = TrapLocal$numLocalIndices,
                       habitatIDDet = TrapLocal$habitatGrid,
                       lengthYCombined = lengthYCombined,
                       trapCov = trapCov)
  
  
  nimData <- list(trapCoords = ScaledtrapCoords,
                  habitat.mx = habitatMask,
                  ones = rep(1, nimConstants$n.individuals),
                  lowerCoords = c(min(coordsHabitatGridCenter[,1]) - 0.5, min(coordsHabitatGridCenter[,2]) - 0.5),
                  upperCoords = c(max(coordsHabitatGridCenter[,1]) + 0.5, max(coordsHabitatGridCenter[,2]) + 0.5))
                  
  
  # We set the parameter values as inits
  nimInits <- list(lambda = lambda,
                   psi = psi,
                   sigma = sigma,
                   betaTrap = betaTrap)
  
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData,
                        inits = nimInits,
                        check = F,
                        calculate = T)
  # FIRST WE GET THE NODES TO SIMULATE
  nodesToSim <- model$getDependencies(c("sxy", "z"), self=T)
  # THEN WE SIMULATE THOSE NODES
  set.seed(100)
  model$simulate(nodesToSim, includeData = FALSE)
  N <- sum(model$z)
  nimData1 <- nimData
  nimData1$y <- model$y
  nimInits1 <- nimInits
  nimInits1$z <- model$z
  nimInits1$sxy <- model$sxy
  
  
  # CREATE AND COMPILE THE NIMBLE MODEL
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData1,
                        inits = nimInits1,
                        check = F,
                        calculate = F)
  logProb <- model$calculate()#-1963.177
  
  ## test R logProp
  expect_match(sprintf("%.10f", logProb),## USE THIS TO THE GET DIGITS VALUES
               "-4575.2386357219",
               info = paste0("Incorrect dpoisLocal_normal LogProb in R"))
  
  cmodel <- compileNimble(model)
  ClogProb <- cmodel$calculate()#-1963.177
  ## test c++ logProp
  expect_match(sprintf("%.10f", ClogProb),
               "-4575.2386357219",
               info = paste0("Incorrect dpoisLocal_normal LogProb in C++"))
  
  
  MCMCconf <- configureMCMC(model = model,
                            monitors = c("N", "sigma", "lambda","psi", "betaTrap"),
                            control = list(reflective = TRUE),
                            thin = 1)
  
  ## Add block sampling of sxy coordinates
  MCMCconf$removeSamplers("sxy")
  ACnodes <- paste0("sxy[", 1:nimConstants$n.individuals, ", 1:2]")
  for(node in ACnodes) {
    MCMCconf$addSampler(target = node,
                        type = "RW_block",
                        control = list(adaptScaleOnly = TRUE),
                        silent = TRUE)
  }
  
  MCMC <- buildMCMC(MCMCconf)
  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
  # RUN THE MCMC
  set.seed(100)
  myNimbleOutput <- runMCMC( mcmc = cMCMC,
                             nburnin = 1000,
                             niter = 5000,
                             nchains = 3,
                             samplesAsCodaMCMC = TRUE)
  sum <- summary(myNimbleOutput)
  ## COMPARE N
  expect_equal(sum$statistics["N","Mean"] ,
               N,
               tolerance=1,
               info = paste0("Non-comparable mean estimated N from simulation with dpoisLocal_normal"))
  
  #COMPARE "lambda","psi","sigma"
  True <- sum$statistics[c("lambda","psi","sigma","betaTrap[1]","betaTrap[2]"),"Mean"]
  True[] <- c(lambda, psi, sigma, betaTrap)
  
  expect_equal(as.numeric(sum$statistics[c("lambda","psi","sigma"),"Mean"]),
               as.numeric(True[c("lambda","psi","sigma")]),
               tolerance = 0.07,
               info = paste0("Non-comparable mean estimated lambda, psi, sigma from simulation with dpoisLocal_normal"))
  
  
  expect_equal(as.numeric(sum$statistics[c("betaTrap[1]","betaTrap[2]"),"Mean"]) ,
               as.numeric(True[c("betaTrap[1]","betaTrap[2]")]),
               tolerance = 0.07,
               info = paste0("Non-comparable mean estimated betaTrap from simulation with dpoisLocal_normal"))
})


