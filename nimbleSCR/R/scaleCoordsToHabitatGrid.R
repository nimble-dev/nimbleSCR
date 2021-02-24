#' Scale x- and y-coordinates to grid cells coordinates.
#'
#' R utility function to scale x- and y- coordinates to the habitat grid.
#' Scaling the coordinates to the habitat grid allows implementation of the fast look-up approach to identify the habitat grid cell in which a point is located. 
#' This technique was first applied by Mike Meredith in SCR (\href{https://mmeredith.net/blog/1309_SECR_in_JAGS_patchy_habitat.htm}{https://mmeredith.net/blog/1309_SECR_in_JAGS_patchy_habitat.htm}).
#' Re-scaling the entire coordinate system of the data input is a requirement to run SCR models with the local evaluation approach. 
#' This function requires square grid cells and coordinates using projection with units in meters or km (e.g., UTM but not latitude/longitude)
#'
#' @param coordsData A matrix or array of x- and y-coordinates to be scaled to the habitat grid. x- and y- coordinates must be identified using "x" and "y" dimnames.
#' @param coordsHabitatGridCenter A matrix of x- and y-coordinates for each habitat grid cell center. 
#' @param scaleToGrid Defaults to TRUE. If FALSE, coordsData are already scaled and will be rescaled to its original coordinates.
#'
#' @return This function returns a list of objects:
#' \itemize{
#' \item coordsDataScaled: A matrix or array of scaled (rescaled if scaleToGrid==FALSE) x- and y-coordinates for coordsData.
#' \item coordsHabitatGridCenterScaled: A matrix of scaled x- and y-cell coordinates for coordsHabitatGridCenter.
#' }
#'
#' @author Richard Bischof, Cyril Milleret
#'
#' @importFrom graphics plot points
#'
#' @examples
#'coordsGridCenter <- expand.grid(list(x = seq(50.5, 100, 1),
#'                                     y = seq(100.5, 150, 1)))
#'coordsData <- expand.grid(list(x = seq(60, 90, 1),
#'                              y = seq(110, 140, 1)))
#'plot(coordsGridCenter[,1]~coordsGridCenter[,2])
#'points(coordsData[,1]~coordsData[,2], col="red")
#'scaled <- scaleCoordsToHabitatGrid(coordsData = coordsData
#'                                   , coordsHabitatGridCenter = coordsGridCenter)
#'plot(scaled$coordsHabitatGridCenterScaled[,1]~scaled$coordsHabitatGridCenterScaled[,2])
#'points(scaled$coordsDataScaled[,1]~scaled$coordsDataScaled[,2], col="red")
#'
#' @export
scaleCoordsToHabitatGrid <- function(coordsData = coordsData,
                                     coordsHabitatGridCenter = coordsHabitatGridCenter,
                                     scaleToGrid = TRUE
){
  
  ## check that arrays have x and y dimnames 
  whereXYdata <- unlist(lapply(dimnames(coordsData), function(x) any(x %in% c("x","y"))))
  if(sum(whereXYdata)==0){stop("Provide the dimnames 'x' and 'y' for coordsData")}
  whereXYGrid <- unlist(lapply(dimnames(coordsHabitatGridCenter), function(x) any(x %in% c("x","y"))))
  if(sum(whereXYGrid)==0){stop("Provide the dimnames 'x' and 'y' for coordsHabitatGridCenter")}
  
  ## Identify the dimension coding for x and y (coordsData may be vastly different depending on the application)
  whereYdataList <- lapply(1:length(dimnames(coordsData)), function(x){
    where <- which(dimnames(coordsData)[[x]] == c("y"))
    if(length(where)==0){where <- 1:dim(coordsData)[x]}
    return(where)
  })
  
  whereXdataList <- lapply(1:length(dimnames(coordsData)), function(x){
    where <- which(dimnames(coordsData)[[x]] == c("x"))
    if(length(where)==0){where <- 1:dim(coordsData)[x]}
    return(where)
  })
  
  ## Find habitat resolution 
  resolution <- min(diff(unique(sort(coordsHabitatGridCenter[ ,"x"]))))#---assumes square grid cells and utm projection (units in meters or km; not latlong!)
  
  ## Obtain x and y min
  start0.y <- max(coordsHabitatGridCenter[ ,"y"]) + resolution/2 #---because we are moving from top to bottom
  start0.x <- min(coordsHabitatGridCenter[ ,"x"]) - resolution/2 #---because we are moving from left to right
  
  ## Re-projectthe grid cell centers
  coordsHabitatGridCenterScaled <- coordsHabitatGridCenter
  coordsHabitatGridCenterScaled[ ,"y"] <- (start0.y - coordsHabitatGridCenter[ ,"y"])/resolution
  coordsHabitatGridCenterScaled[ ,"x"] <- (coordsHabitatGridCenter[ ,"x"] - start0.x)/resolution
  
  ## Re-project data coordinates 
  coordsDataScaled <- coordsData
  ## If coordsData has 2 dimensions    
  if(length(dim(coordsData)) == 2){
    if(scaleToGrid){
      coordsDataScaled[ ,"y"] <- (start0.y - coordsDataScaled[ ,"y"])/resolution
      coordsDataScaled[ ,"x"] <- (coordsDataScaled[ ,"x"] - start0.x)/resolution 
    } else {
      coordsDataScaled[ ,"y"] <- start0.y - coordsData[ ,"y"] * resolution 
      coordsDataScaled[ ,"x"] <- start0.x + coordsData[ ,"x"] * resolution 
    }
    
  }
  
  ## If data has 3 dimensions    
  if(length(dim(coordsData)) == 3){
    Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] 
    X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]]  
    if(scaleToGrid){
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] <- (start0.y - Y)/resolution
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]] <- (X - start0.x)/resolution 
    } else {
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] <- start0.y - Y * resolution 
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]] <- start0.x + X * resolution 
    }
  }
  
  ## If data has 4 dimensions    
  if(length(dim(coordsData)) == 4){
    Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] 
    X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]]  
    if(scaleToGrid){
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] <- (start0.y - Y)/resolution
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]] <- (X - start0.x)/resolution 
    } else {
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] <- start0.y - Y * resolution 
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]] <- start0.x + X * resolution 
    }
  }
  
  ## If data has 5 dimensions    
  if(length(dim(coordsData)) == 5){
    Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] 
    X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]]  
    if(scaleToGrid){
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] <- (start0.y - Y)/resolution
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]] <- (X - start0.x)/resolution 
    } else {
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] <- start0.y - Y * resolution 
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]] <- start0.x + X * resolution 
    }
  }
  
  if(length(dim(coordsData)) >= 6){
    stop("The function can only handle arrays of up to 5 dimensions for now...sorry")
  }
  
  return(list(coordsDataScaled = coordsDataScaled,
              coordsHabitatGridCenterScaled = coordsHabitatGridCenterScaled))
}