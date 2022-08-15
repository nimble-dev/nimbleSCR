#' Local Objects Identification
#'
#' R utility function to identify all objects (e.g. traps) within a given radius dmax of each cell in a habitat mask.
#' Used in the implementation of the local evaluation approach in SCR models (\code{\link{dbinomLocal_normal}};\code{\link{dpoisLocal_normal}}).
#' The distance to the activity center and the detection probability are then calculated for local objects only (i.e. the detection probability is assumed to be 0 
#' for all other objects as they are far enough from the activity center).
#'
#' The \code{getLocalObjects} function is used in advance of model building.
#'
#' @param habitatMask a binary matrix object indicating which cells are considered as suitable habitat.
#' @param coords A matrix giving the x- and y-coordinate of each object (i.e. trap). x- and y- coordinates should be scaled to the habitat (\code{\link{scaleCoordsToHabitatGrid}}). 
#' @param dmax The maximal radius from a habitat cell center within which detection probability is evaluated locally for each trap.
#' @param resizeFactor An aggregation factor to reduce the number of habitat cells to retrieve local objects for. Defaults to 1; no aggregation.
#' @param plot.check A visualization option (if TRUE); displays which objects are considered "local objects" for a randomly chosen habitat cell.
#'
#' @return This function returns a list of objects:
#' \itemize{
#' \item localIndices: a matrix with number of rows equal to the reduced number of habitat grid cells (following aggregation).
#' Each row gives the id numbers of the local objects associated with this grid cell.
#' \item habitatGrid: a matrix of habitat grid cells ID corresponding to the row indices in localIndices. 
#' \item numLocalIndices: a vector of the number of local objects for each habitat grid cell in habitatGrid.
#' \item numLocalIndicesMax: the maximum number of local objects for any habitat grid cell ; corresponds to the number of columns in habitatGrid.
#' \item resizeFactor: the aggregation factor used to reduce the number of habitat grid cells.
#' }
#'
#' @author Cyril Milleret and Pierre Dupont
#'
#' @importFrom graphics plot points
#'
#' @examples
#' colNum <- sample(20:100,1)
#' rowNum <- sample(20:100,1)
#' coords <- expand.grid(list(x = seq(0.5, colNum, 1),
#'                                y = seq(0.5, rowNum, 1)))
#' 
#' habitatMask <- matrix(rbinom(colNum*rowNum, 1, 0.8), ncol = colNum, nrow = rowNum)
#' 
#' localObject.list <- getLocalObjects(habitatMask, coords,  dmax = 7,resizeFactor = 1)
#'
#' @export
getLocalObjects <- function( habitatMask,
                             coords,
                             dmax,
                             resizeFactor = 1,
                             plot.check = TRUE
){
  ## STORE THE COORDINATES OF THE ORIGINAL HABITAT CELLS
  oldCoords <- which(habitatMask > 0, arr.ind = T) - 0.5
  oldCoords <- cbind(oldCoords[,2], oldCoords[,1])
  names(oldCoords) <- c("x", "y")
  
  ## GET ORIGIN FOR THE NEW (RESIZED) HABITAT MATRIX 
  origin <- (resizeFactor/2) 
  
  ## GET DIMENSIONS FOR THE NEW HABITAT MATRIX
  xMax <- ceiling(dim(habitatMask)[2]/resizeFactor) * resizeFactor
  yMax <- ceiling(dim(habitatMask)[1]/resizeFactor) * resizeFactor
  
  ## GET COORDINATES FOR THE NEW HABITAT CELLS
  xCoords <- seq(origin, xMax, by = resizeFactor)
  yCoords <- seq(origin, yMax, by = resizeFactor)
  habitatCoords <- expand.grid(list(x = xCoords, y = yCoords))
  habitatCoords <- habitatCoords[order(habitatCoords[ ,1], habitatCoords[ ,2]), ]
  
  ## GET UPPER AND LOWER COORDINATES FOR THE NEW HABITAT CELLS
  habUpCoords <- habitatCoords + resizeFactor/2
  habLoCoords <- habitatCoords - resizeFactor/2
  
  ## CHECK WHICH "NEW CELLS" CONTAIN AT LEAST ONE "OLD CELL"
  isIn <- unlist(lapply(1:dim(habUpCoords)[1], function(c){
    sum((habLoCoords[c,1] <= oldCoords[ ,1]) *
          (habLoCoords[c,2] <= oldCoords[ ,2]) *
          (habUpCoords[c, 1] > oldCoords[ ,1]) *
          (habUpCoords[c,2] > oldCoords[ ,2])) > 0
  })) 
  
  ## REMOVE NEW HABITAT CELL COORDINATES THAT ARE NOT HABITAT  
  habitatCoords <- habitatCoords[isIn, ]
  
  ## CREATE AN EMPTY MATRIX OF NEW HABITAT CELL IDs
  habitatID <- matrix(0, nrow = length(yCoords), ncol = length(xCoords))
  for(c in 1:dim(habitatCoords)[1]){
    habitatID[trunc(habitatCoords[c,2]/resizeFactor)+1, trunc(habitatCoords[c,1]/resizeFactor)+1] <- c
  }
  
  ## DETERMINE WHICH POINTS ARE WITHIN dmax OF THE CENTER OF EACH NEW HABITAT CELL
  localIndices  <- apply(habitatCoords, 1, function(x){
    D <- sqrt((x[1] - coords[,1])^2 + (x[2] - coords[,2])^2) 
    which(D < dmax)
  })
  
  ## MAKE SURE IT ALWAYS RETURN A LIST
  if(inherits(localIndices,"matrix")){
    localIndices <- lapply(1:dim(localIndices)[2], function(x) localIndices[ ,x])
  }
  
  ## GET NUMBER OF DETECTORS WITHIN dmax OF EACH NEW HABITAT CELL
  numLocalIndices <- unlist(lapply(localIndices, function(x) length(x)))
  maxLocalIndices <- max(numLocalIndices)
  
  ## FOR ALL HABITAT GRIDS, THE LOCAL EVALUATION SHOULD BE LARGE ENOUGH TO OVERLAP WITH >0 TRAP
  if(any(numLocalIndices %in% 0 )){
    stop("dmax value too small. All habitat grid cells should have at least one local object within a radius of dmax.")
  }
  ## STORE LOCAL  INDICES IN A MATRIX 
  Index <- matrix(0, nrow = length(localIndices), ncol = maxLocalIndices)
  for(j in 1:length(localIndices)){
    if(length(localIndices[[j]])!=0){
      Index[j, 1:numLocalIndices[j]] <- localIndices[[j]]
    }
  }
  
  
  ## PLOT CHECK 
  if(plot.check){
    SXY <- as.numeric(habitatCoords[sample(1:dim(habitatCoords)[1], size = 1), ])
    sxyID <- habitatID[trunc(SXY[2]/resizeFactor)+1, trunc(SXY[1]/resizeFactor)+1]
    index <- Index[sxyID, 1:numLocalIndices[sxyID]]
    
    plot(habitatCoords[ ,2] ~ habitatCoords[ ,1], pch = 16, cex = 0.1)
    points(habitatCoords[sxyID,2] ~ habitatCoords[sxyID,1], pch = 16, cex = 0.4, col = "orange")
    points(coords[ ,2] ~ coords[ ,1], pch = 16, cex = 0.2, col = "red")
    points(coords[index,2] ~ coords[index,1], pch = 16, cex = 0.4, col = "blue")
    points(SXY[2] ~ SXY[1], bg = "red", pch = 21, cex = 1.2)
  }
  
  
  
  ## OUTPUT LIST
  output <- list( habitatGrid = habitatID,
                  localIndices = Index,
                  numLocalIndices = numLocalIndices,
                  numLocalIndicesMax = maxLocalIndices,
                  resizeFactor = resizeFactor)
  
  return(output)
}
