#' Local Trap Identification
#'
#' R utility function to find the traps/detector locations within a given radius of all the cells in a habitat mask.
#'
#' The \code{getLocalTraps} function is used in advance of model building.
#'
#' A visualization option (plot.check) is also provided, which displays which traps are considered "local traps"
#' around one habitat grid cell selected at random.
#'
#' @param habitatMask a binary matrix object indicating which cells are considered as suitable habitat.
#' @param trapCoords A matrix giving the x- and y-coordinate locations of all traps.
#' @param dmax The maximal radius from a cell center for performing local trap calculations
#' @param resizeFactor An aggregation factor to reduce the number of habitat cells to retrieve local traps for. 
#' @param plot.check A visualization option (if TRUE), which displays which traps are considered "local traps"
#' around one habitat grid cell selected at random.
#'
#' @return It returns a list of objects:
#' - localTrapIndices : a matrix with number of rows equal to the reduced number of habitat grid cells.
#' Each row gives the id numbers of the local traps for this grid cell.
#' - habitatGrid : a matrix of habitat grid cells ID corresponding to the row indices in localTrapIndices. 
#' - numLocalTraps: a vector of the number of local traps for each habitat grid cell in habitatGrid.
#' - numLocalTrapsMax: the maximum number of local traps for any habitat grid cell ; corresponds to the number of columns in habitatGrid.
#' - resizeFactor: the aggregation fatcor used to reduces the number of habitat grid cells ; default = 1 means no aggregation.
#' 
#' @author Cyril Milleret
#'
#' @importFrom stats aggregate
#' @importFrom graphics plot points
#'
#' @examples
#' \donttest{
#' 
#' habitatMask <- matrix(rbinom(10000,1,0.75), nrow = 100)
#' trapCoords <- expand.grid(list(1:100,1:100))
#' localTraps.list <- getLocalTraps(habitatMask, trapCoords, dmax = 6)
#' }
#'
#' @export
getLocalTraps <- function( habitatMask,
                           trapCoords,
                           dmax,
                           resizeFactor = 1,
                           plot.check = TRUE
){
  ## ==== 1. CREATE THE XY COORDINATES OF THE HABITAT ====
  habitatID <- habCoordsx <- habCoordsy <- habitatMask 
  if(resizeFactor == 1){
    from <- 0.5
    resizeFactor <- 1 
  } else {  
    from <- (resizeFactor/2) + 0.5
  }
  # GET DIMENSIONS MATRIX
  dimCoords <- c(ceiling(dim(habCoordsx)[1]/resizeFactor), ceiling(dim(habCoordsx)[2]/resizeFactor))
  CoordsMax <- dimCoords * resizeFactor
  
  habCoordsx <- matrix(rep(seq(from, CoordsMax[2], by = resizeFactor ), dimCoords[1]),
                       nrow = dimCoords[1],
                       ncol = dimCoords[2],
                       byrow = T)
  
  habCoordsy <- matrix(rep(seq(from, CoordsMax[1], by = resizeFactor ), dimCoords[2]),
                       nrow = dimCoords[1],
                       ncol = dimCoords[2],
                       byrow = F)
  habCoordsy <- as.vector(habCoordsy)
  habCoordsx <- as.vector(habCoordsx)
  habCoordsxy <- cbind(habCoordsx, habCoordsy)
  
  
  ## ==== 2. RESCALE THE HABITAT ====
  r <- raster(habitatMask)
  if(resizeFactor > 1){ r <- aggregate(r, fact = resizeFactor)}
  r[r > 0] <-1
  habitatMask1 <- as.matrix(r)
  
  # CREATE A HABITAT ID MATRIX 
  habitatID <- habitatMask1
  habitatID[] <- as.character(habitatMask1)
  
  # ONLY GIVE AN ID TO THE HABITAT CELLS == 1 
  habitatID[habitatMask1 == "1"] <- 1:sum(habitatMask1 == "1")
  m <- sapply(habitatID, FUN = as.numeric)
  habitatID <- matrix(m, nrow = dim(habitatID)[1], ncol = dim(habitatID)[2], byrow = F)
  
  # REMOVE HABITAT CELL COORDINATES THAT ARE NOT HABITAT  
  habCoordsxy  <- habCoordsxy[as.character(as.vector(habitatMask1))=="1",]
  
  ## ==== 3. FIND DETECTORS THAT ARE WITHIN dmax OF EACH HABITAT CELL  ====
  # Determine detector within radius distance from the center of each habitat cell
  localTrapsIndices  <- apply(habCoordsxy, 1, function(x){
    D <- sqrt((x[1] - trapCoords[,1])^2 + (x[2] - trapCoords[,2])^2) 
    which(D < dmax)
  })
  
  # Make sure it always returns a list. 
  if(class(localTrapsIndices) == "matrix"){
    localTrapsIndices <- lapply(1:dim(localTrapsIndices)[2], function(x) localTrapsIndices[ ,x])
  }
  
  ## ==== 4. STORE LOCAL DETECTOR INDICES IN A MATRIX ====
  # get number of detectors within dmax of each cell
  numLocalTraps <- unlist(lapply(localTrapsIndices, function(x) length(x)))
  maxLocalTraps <- max(numLocalTraps)
  # store detector index (colums) for each habitat cell (rows)
  detectorIndex <- matrix(0, nrow = length(localTrapsIndices), ncol = maxNBDets)
  for(j in 1:length(localTrapsIndices)){
    if(length(localTrapsIndices[[j]])!=0){
      detectorIndex[j, 1:numLocalTraps[j]] <- localTrapsIndices[[j]]
    }
  }
  
  ## PLOT CHECK 
  if(plot.check){
    SXY <- habCoordsxy[sample(1:dim(habCoordsxy)[1], size = 1), ]
    sxyID <- habitatID[trunc(SXY[2]/resizeFactor)+1, trunc(SXY[1]/resizeFactor)+1]

    index <- detectorIndex[sxyID,1:numLocalTraps[sxyID]]
    plot(habCoordsxy[,2] ~ habCoordsxy[,1], pch=16, cex=0.1)
    points(habCoordsxy[sxyID,2]~habCoordsxy[sxyID,1], pch=16, cex=0.4, col="orange")
    
    points(trapCoords[,2]~trapCoords[,1], pch=16, cex=0.2, col="red")
    points(trapCoords[index,2]~trapCoords[index,1], pch=16, cex=0.4, col="blue")
    points(SXY[2]~SXY[1], bg="red", pch=21, cex=1.2)
    
  }
  
  output <-list( habitatGrid = habitatID,
                 localTrapsIndices = detectorIndex,
                 numLocalTraps = numLocalTraps,
                 numLocalTrapsMax = maxNBDets,
                 resizeFactor = resizeFactor)
  
  return(output)
}
