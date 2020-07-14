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
#' @param plot.check A visualization option (if TRUE); displays which traps are considered "local traps"
#' around one randomly selected habitat grid cell.
#'
#' @return It returns a list of objects:
#' - localTrapIndices : a matrix with number of rows equal to the reduced number of habitat grid cells.
#' Each row gives the id numbers of the local traps for this grid cell.
#' - habitatGrid : a matrix of habitat grid cells ID corresponding to the row indices in localTrapIndices. 
#' - numLocalTraps: a vector of the number of local traps for each habitat grid cell in habitatGrid.
#' - numLocalTrapsMax: the maximum number of local traps for any habitat grid cell ; corresponds to the number of columns in habitatGrid.
#' - resizeFactor: the aggregation fatcor used to reduces the number of habitat grid cells ; default = 1 means no aggregation.
#' 
#' @author Pierre Dupont
#'
#' @importFrom graphics plot points
#'
#' @examples
#' \donttest{
#' colNum <- sample(20:100,1)
#' rowNum <- sample(20:100,1)
#' trapCoords <- expand.grid(list(x = seq(0.5, colNum, 1),
#'                                y = seq(0.5, rowNum, 1)))
#' 
#' habitatMask <- matrix(rbinom(colNum*rowNum, 1, 0.8), ncol = colNum, nrow = rowNum)
#' 
#' localTraps.list <- getLocalTraps(habitatMask, trapCoords, resizeFactor = 1, dmax = 7)
#' }
#'
#' @export
getLocalTraps <- function( habitatMask,
                           trapCoords,
                           dmax,
                           resizeFactor = 1,
                           plot.check = TRUE
){
  ## STORE THE COORDINATES OF THE ORIGINAL HABITAT CELLS
  oldCoords <- which(habitatMask == 1, arr.ind = T)-0.5
  
  ## GET ORIGIN FOR THE NEW (RESIZED) HABITAT MATRIX
  origin <- (resizeFactor/2) 

  ## GET DIMENSIONS FOR THE NEW HABITAT MATRIX
  xMax <- ceiling(dim(habitatMask)[2]/resizeFactor) * resizeFactor
  yMax <- ceiling(dim(habitatMask)[1]/resizeFactor) * resizeFactor

  ## GET COORDINATES FOR THE NEW HABITAT CELLS
  xCoords <- seq(origin, xMax, by = resizeFactor)
  yCoords <- seq(origin, yMax, by = resizeFactor)
  habitatCoords <- expand.grid(list(x = xCoords, y = yCoords))
  habitatCoords <- habitatCoords[order(habitatCoords[ ,1],habitatCoords[ ,2]), ]
  
  ## GET UPPER AND LOWER COORDS FOR THE NEW HABITAT CELLS
  habUpCoords <- habitatCoords + resizeFactor/2
  habLoCoords <- habitatCoords - resizeFactor/2
  
  # check which "new cells" contains at least one "old cell"
  # Assess the intersection between the coordinates and the observation windows
  # Initialise an output vector with each element set to 1
  isIn <- NULL
  for(c in 1:dim(habUpCoords)[1]){
    isIn[c] <- sum((habLoCoords[c,1] <= oldCoords[ ,2]) * (habLoCoords[c,2] <= oldCoords[ ,1]) *
      (habUpCoords[c, 1] > oldCoords[ ,2]) * (habUpCoords[c,2] > oldCoords[ ,1])) > 0
      }#c
  
  ## REMOVE NEW HABITAT CELL COORDINATES THAT ARE NOT HABITAT  
  habitatCoords <- habitatCoords[isIn, ]
  
  ## CREATE AN EMPTY MATRIX OF NEW HABITAT CELL IDs
  habitatID <- matrix(0, nrow = length(yCoords), ncol = length(xCoords))
  for(c in 1:dim(habitatCoords)[1]){
    habitatID[trunc(habitatCoords[c,2]/resizeFactor)+1, trunc(habitatCoords[c,1]/resizeFactor)+1] <- c
  }
  
  ## DETERMINE WHICH DETECTORS ARE WITHIN dmax OF THE CENTER OF EACH NEW HABITAT CELL
  localTrapsIndices  <- apply(habitatCoords, 1, function(x){
    D <- sqrt((x[1] - trapCoords[,1])^2 + (x[2] - trapCoords[,2])^2) 
    which(D < dmax)
  })
  
  ## MAKE SURE IT ALWAYS RETURN A LIST
  if(class(localTrapsIndices) == "matrix"){
    localTrapsIndices <- lapply(1:dim(localTrapsIndices)[2], function(x) localTrapsIndices[ ,x])
  }
  
  ## GET NUMBER OF DETECTORS WITHIN dmax OF EACH NEW HABITAT CELL
  numLocalTraps <- unlist(lapply(localTrapsIndices, function(x) length(x)))
  maxLocalTraps <- max(numLocalTraps)
  
  ## STORE LOCAL DETECTOR INDICES IN A MATRIX 
  detectorIndex <- matrix(0, nrow = length(localTrapsIndices), ncol = maxLocalTraps)
  for(j in 1:length(localTrapsIndices)){
    if(length(localTrapsIndices[[j]])!=0){
      detectorIndex[j, 1:numLocalTraps[j]] <- localTrapsIndices[[j]]
    }
  }
  
  ## PLOT CHECK 
  if(plot.check){
    SXY <- as.numeric(habitatCoords[sample(1:dim(habitatCoords)[1], size = 1), ])
    sxyID <- habitatID[as.numeric(trunc(SXY[2]/resizeFactor)+1),
                       as.numeric(trunc(SXY[1]/resizeFactor)+1)]
    
    index <- detectorIndex[sxyID,1:numLocalTraps[sxyID]]
    plot(habitatCoords[,2] ~ habitatCoords[,1], pch=16, cex=0.1)
    points(habitatCoords[sxyID,2] ~ habitatCoords[sxyID,1], pch=16, cex=0.4, col="orange")
    
    points(trapCoords[ ,2] ~ trapCoords[,1], pch=16, cex=0.2, col="red")
    points(trapCoords[index,2]~trapCoords[index,1], pch=16, cex=0.4, col="blue")
    points(SXY[2]~SXY[1], bg="red", pch=21, cex=1.2)
    
  }
  
  output <-list( habitatGrid = habitatID,
                 localTrapsIndices = detectorIndex,
                 numLocalTraps = numLocalTraps,
                 numLocalTrapsMax = maxLocalTraps,
                 resizeFactor = resizeFactor)
  
  return(output)
}