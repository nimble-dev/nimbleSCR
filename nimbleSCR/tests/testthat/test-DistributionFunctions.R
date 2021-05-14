library(nimble)
nimbleOptions(verbose = FALSE)
# tempFileName <- 'distributionsTestLog.Rout'
# generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForDistributionsTesting'))
# outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForDistributionsTesting'), goldFileName) else tempFileName
# 
# sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)


if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else stop('unknown user')

#---- 1. REPRODUCE RESULTS FROM PUBLISHED PAPERS ----
#----   1.1 Turek et al 2021 Ecosphere ----
## REPRODUCE WOLVERINE LIKELIHODD CALCULATION FROM THE WOLVERINE EXAMPLE IN DANIEL'S PAPER


# IN R
test_that("Wolverine LogProb is calculated correctly in R",{
  truthWolverineLikelihood <- "-15726.3284159870"# as calculated in the paper
  #Code as in the vignette https://nimble-dev.github.io/nimbleSCR/wolverine_example.html
  load(file.path(baseDir,"nimbleSCR/tests/testthat/WolverineData.RData"))
  code <- nimbleCode({
    ## priors
    psi ~ dunif(0, 1)
    sigma ~ dunif(0, 50)
    p0 ~ dunif(0, 1)
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
      y[i, 1:nMaxDetectors] ~ dbinomLocal_normal( detNums = nbDetections[i],
                                                  detIndices = yDets[i,1:nMaxDetectors],
                                                  size = trials[1:n.detectors],
                                                  p0 = p0,
                                                  s = sxy[i,1:2],
                                                  sigma = sigma,
                                                  trapCoords = detector.xy[1:n.detectors,1:2],
                                                  localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                                  localTrapsNum = nDetectors[1:n.cells],
                                                  resizeFactor = ResizeFactor,
                                                  habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                  indicator = z[i])
    }
    ## derived quantity: total population size
    N <- sum(z[1:n.individuals])
  })

  data <- list(y = my.jags.input$y,
               z = my.jags.input$z,
               detector.xy = my.jags.input$detector.xy,
               habitat.mx = my.jags.input$habitat.mx,
               ones = my.jags.input$OK,
               lowerCoords = c(0,0),
               upperCoords = c(
                 dim(my.jags.input$habitat.mx)[2],
                 dim(my.jags.input$habitat.mx)[1]),
               trials = rep(1, dim(my.jags.input$detector.xy)[1]))
  constants <- list(n.individuals = my.jags.input$n.individuals,
                    n.detectors = dim(my.jags.input$detector.xy)[1],
                    y.max = dim(my.jags.input$habitat.mx)[1],
                    x.max = dim(my.jags.input$habitat.mx)[2])
  inits <- list(sxy = inits.1$sxy,
                z = inits.1$z,
                p0 = 0.05,
                psi = 0.5,
                sigma = 6)
  DetectorIndex <- getLocalObjects(habitatMask = data$habitat.mx,
                                   coords = data$detector.xy,
                                   dmax = 38,
                                   resizeFactor = 24,
                                   plot.check = F)
  constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
  constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
  constants$ResizeFactor <- DetectorIndex$resizeFactor
  constants$n.cells <- dim(DetectorIndex$localIndices)[1]
  constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
  data$detectorIndex <- DetectorIndex$localIndices
  data$nDetectors <- DetectorIndex$numLocalIndices
  data$habitatIDDet <- DetectorIndex$habitatGrid

  ySparse <- getSparseY(x = my.jags.input$y)
  data$y <- ySparse$y[,,1]
  data$yDets <- ySparse$detIndices[,,1]
  data$nbDetections <- ySparse$detNums[,1]
  constants$nMaxDetectors <- ySparse$maxDetNums
  # R  
  Rmodel <- nimbleModel(code, constants, data, inits, calculate = F)
  logProb <- Rmodel$calculate()
 
  expect_match(sprintf("%.10f", logProb),# to get the right digits
               truthWolverineLikelihood,
               info = paste0("Incorrect wolverine LogProb in R"))
  # COMPILED
   cRmodel <- compileNimble(Rmodel)
   clogProb <- cRmodel$calculate()
  
   expect_equal(sprintf("%.10f", logProb),# to get the right digits
                 truthWolverineLikelihood,
                 info = paste0("Incorrect wolverine LogProb in C++"))
  })


#---- 2. SIMULATE AND FIT MODELS ----
#----   2.1 dbinomLocal_normal() WITH sxy~dunif() ----
test_that("Simulate and fit dbinomLocal_normal ",{
  
coordsHabitatGridCenter <- cbind(rep(seq(11.5, 0.5, by=-1), 10),
                                 sort(rep(seq(0.5, 11.5, by=1), 10)))
colnames(coordsHabitatGridCenter) <- c("x","y")

# CREATE A 8*8 TRAP GRID CENTERED ON THE HABITAT GRID
trapCoords <- cbind(rep(seq(2.5, 9.5,by=1),8),
                    sort(rep(seq(2.5, 9.5,by=1),8)))
colnames(trapCoords) <- c("x","y")

# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch=16)
points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16)

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
  p0 ~ dunif(0, 1)
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
    y[i, 1:lengthYCombined] ~ dbinomLocal_normal(# detNums = nbDetections[i],
      # detIndices = yDets[i,1:nMaxDetectors],
      size = trials[1:n.detectors],
      p0 = p0,
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

p0 <- 0.2
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
                upperCoords = c(max(coordsHabitatGridCenter[,1]) + 0.5, max(coordsHabitatGridCenter[,2]) + 0.5),
                trials = rep(1, dim(ScaledtrapCoords)[1])

)
# We set the parameter values as inits
nimInits <- list(p0 = p0,
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
logProb <- model$calculate()#-1963.177

## test R logProp
expect_match(sprintf("%.10f", logProb),## USE THIS TO THE GET DIGITS VALUES
             "-1963.1772292583",
             info = paste0("Incorrect dbinomLocal_normal LogProb in R"))

cmodel <- compileNimble(model)
ClogProb <- cmodel$calculate()#-1963.177
## test c++ logProp
expect_match(sprintf("%.10f", ClogProb),
             "-1963.1772292583",
             info = paste0("Incorrect dbinomLocal_normal LogProb in C++"))


MCMCconf <- configureMCMC(model = model,
                          monitors = c("N", "sigma", "p0","psi"),
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
             info = paste0("Non-comparable mean estimated N from simulation with dbinomLocal_normal"))

#COMPARE "p0","psi","sigma"
True <- sum$statistics[c("p0","psi","sigma"),"Mean"]
True[] <- c(p0, psi, sigma)

expect_equal(sum$statistics[c("p0","psi","sigma"),"Mean"] ,
             True,
             tolerance = 0.07,
             info = paste0("Non-comparable mean estimated p0, psi, sigma from simulation with dbinomLocal_normal"))
})




