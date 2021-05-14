library(nimble)
nimbleOptions(verbose = FALSE)

#---- 1.TEST UTILITY FUNCTIONS  ----
#----   1.1 getSparseY() ----

test_that("Number of spatial recaptures correclty retrieved by getSparseY",{
  set.seed(10)
  y.full <- matrix(rbinom(5000, 5, 0.02), ncol = 100)
  sum(y.full)
  y <- getSparseY(y.full)
  
  expect_equal(sum(y.full>0),
               sum(y$detNums),
               info = paste0("Incorrect number of spatial recaptures retrieved by getSparseY"))
  
  expect_equal(which(y.full[10,]>0),
               y$detIndices[10,1:y$detNums[10,1],1],
               info = paste0("Incorrect setector ids retrieved by getSparseY"))              
               })
#----   1.2 scaleCoordsToHabitatGrid() ----
#----     1.2.1 with dim(coordsData) == 2  ----



#
test_that("scaleCoordsToHabitatGrid correctly scale and rescale coordinates",{
  coordsGridCenter <- expand.grid(list(x = seq(50.5, 100, 1),
                                       y = seq(100.5, 150, 1)))
  coordsData <- expand.grid(list(x = seq(60, 90, 1),
                                 y = seq(110, 140, 1)))
  scaled <- scaleCoordsToHabitatGrid(coordsData = coordsData
                                     , coordsHabitatGridCenter = coordsGridCenter)
  
  reScaled <- scaleCoordsToHabitatGrid(coordsData = scaled$coordsDataScaled
                                       , coordsHabitatGridCenter = coordsGridCenter,
                                       scaleToGrid = F)
  
  expect_equal(sum(reScaled$coordsDataScaled==coordsData),
               prod(dim(reScaled$coordsDataScaled)),
               info = paste0("scaleCoordsToHabitatGrid incorrectly scale and rescale coordinates"))})

#----     1.2.2 with dim(coordsData) == 3  ----


test_that("scaleCoordsToHabitatGrid correctly scale and rescale coordinates with >2 dimensions",{
  coordsGridCenter <- expand.grid(list(x = seq(50.5, 100, 1),
                                       y = seq(100.5, 150, 1)))
  coordsData <- expand.grid(list(x = seq(60, 90, 1),
                                 y = seq(110, 140, 1)))
  coordsData.ar <- array(0, c(dim(coordsData),3))
  dimnames(coordsData.ar)[2] <- list(c("x","y"))
  coordsData.ar[,,1] <- as.matrix(coordsData)
  coordsData.ar[,,2] <- as.matrix(coordsData)
  coordsData.ar[,,3] <- as.matrix(coordsData)
  
  scaled.ar <- scaleCoordsToHabitatGrid(coordsData = coordsData.ar
                                        , coordsHabitatGridCenter = coordsGridCenter)
  
  reScaled.ar <- scaleCoordsToHabitatGrid(coordsData = scaled.ar$coordsDataScaled
                                          , coordsHabitatGridCenter = coordsGridCenter
                                          , scaleToGrid = F)
  
  expect_equal(sum(reScaled.ar$coordsDataScaled==coordsData.ar),
               prod(dim(reScaled.ar$coordsDataScaled)),
               info = paste0("scaleCoordsToHabitatGrid incorrectly scale and rescale coordinates with >2 dimensions"))})

#----   1.3 getLocalObjects() ----
#----     1.3.1 check that a list of local detetors is provided for each habitat==1 ----

test_that("getLocalObjects correctly retrieved local indices",{
  set.seed(5)
  colNum <- sample(20:100,1)
  set.seed(5)
  rowNum <- sample(20:100,1)
  coords <- expand.grid(list(x = seq(0.5, colNum, 1),
                             y = seq(0.5, rowNum, 1)))
  set.seed(5)
  habitatMask <- matrix(rbinom(colNum*rowNum, 1, 0.8), ncol = colNum, nrow = rowNum)
  dmax <- 7
  localObject.list <- getLocalObjects(habitatMask, coords,  dmax = dmax, resizeFactor = 1, plot.check = F)
  expect_equal(sum(habitatMask),
               nrow(localObject.list$localIndices),
               info = paste0("getLocalObjects correctly retrieved local indices"))
  
  #----     1.3.2 check that local detectors are located in distance < dmax----
  s <- coords[40,]
  cellID <- localObject.list$habitatGrid[trunc(s[1,2])+1, trunc(s[1,1])+1 ]
  localTraps <- localObject.list$localIndices[cellID,]
  coordsLocalTraps <- coords[localTraps,]
  D <- matrix(unlist(lapply(1:nrow(coordsLocalTraps),
                            function(dd){sqrt((s[1] - coordsLocalTraps[dd,1])^2 + (s[2]  - coordsLocalTraps[dd,2])^2) })),
              nrow(coordsLocalTraps), nrow(coordsLocalTraps), byrow = TRUE)
  
  expect_equal(sum(D < dmax),
               prod(dim(D)),
               info = paste0("Local objects correctly retrieved < dmax"))
  
  
  
  })



# test_that("Local objects correctly retrieved < dmax",{
#   s <- coords[40,]
#   cellID <- localObject.list$habitatGrid[trunc(s[1,2])+1, trunc(s[1,1])+1 ]
#   localTraps <- localObject.list$localIndices[cellID,]
#   coordsLocalTraps <- coords[localTraps,]
#   D <- matrix(unlist(lapply(1:nrow(coordsLocalTraps),
#                             function(dd){sqrt((s[1] - coordsLocalTraps[dd,1])^2 + (s[2]  - coordsLocalTraps[dd,2])^2) })),
#               nrow(coordsLocalTraps), nrow(coordsLocalTraps), byrow = TRUE)
#   
#   expect_equal(sum(D < dmax),
#                prod(dim(D)),
#                info = paste0("Local objects correctly retrieved < dmax"))})






