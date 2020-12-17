#' Scale x and y coordinates to the habitat grid  
#'
#' R utility function to scale x and y coordinates to the habitat grid.
#' Rescaling the coordinates to the habitat allows a fast look up of a matrix of habitat cellID from x and y coordinates. 
#' This technique was first applied by Mike Meredith in SCR (\href{https://mmeredith.net/blog/1309_SECR_in_JAGS_patchy_habitat.htm}).
#' Rescaling the entire coordinate system of the data input is a requirement to run SCR models with the (\code{\link{nimbleSCR}}) package.  
#' This square grid cells and utm projection (units in meters or km; not latlong!)
#' 
#' Using \code{scaleCoordsToHabitatGrid} to rescale the entire coordinate system of the data input is a requirement to run SCR models with the (\code{\link{nimbleSCR}}) package.  
#'
#' @param coordsData A matrix or an array of the x and y data to scale to the habitat grid. Is is necessary to provide the dimnames of the x and y data such as c("x", "y")
#' @param coordsHabitatGridCenter A matrix with the x and y coordinates of each habitat  grid cell center. 
#' @param scaleToGrid If FALSE, the coordsData are already scaled and will be rescaled to its original coordinates.
#'
#' @return This function returns a list of objects:
#' \itemize{
#' \item coordsDataScaled: A matrix or an array of the scaled (rescaled if scaleToGrid==FALSE) x and y coordinates coordsData
#' \item coordsHabitatGridCenterScaled: A matrix of the scaled x and y cell coordinate centers coordsHabitatGridCenter
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
scaleCoordsToHabitatGrid <- function(coordsData = coordsData
                             ,
                              coordsHabitatGridCenter = coordsHabitatGridCenter
                             ,
                              scaleToGrid = TRUE
                             
){
  
  ## check that arrays have a x and y dimnames 
  whereXYdata <- unlist(lapply(dimnames(coordsData), function(x) any(x %in% c("x","y"))))
  if(sum(whereXYdata)==0){stop("Provide the dimnames 'x' and 'y' for coordsData")}
  whereXYGrid <- unlist(lapply(dimnames(coordsHabitatGridCenter), function(x) any(x %in% c("x","y"))))
  if(sum(whereXYGrid)==0){stop("Provide the dimnames 'x' and 'y' for coordsHabitatGridCenter")}
  
  
  ## Find automatically on which dimension are x and y
  ## when there are more than 2 dimensions, it may not be obvious (using posterior sxy coordinates...)
  whereYdataList <- lapply(1:length(dimnames(coordsData)), function(x){
    where <- dimnames(coordsData)[[x]] == c("y")
    if(length(where)==0){where <- 1:dim(coordsData)[x]}
    return(where)
  }
    )
  whereXdataList <- lapply(1:length(dimnames(coordsData)), function(x){
    where <- dimnames(coordsData)[[x]] == c("X")
    if(length(where)==0){where <- 1:dim(coordsData)[x]}
    return(where)
  }
  )
  
  ## find habitat resolution 
  resolution <- min(diff(unique(sort(coordsHabitatGridCenter[ ,"x"]))))#---assumes square grid cells and utm projection (units in meters or km; not latlong!)
  
  ## obtain x and y min
  start0.y <- max(coordsHabitatGridCenter[ ,"y"]) + resolution/2 #---because we are moving from top to bottom
  start0.x <- min(coordsHabitatGridCenter[ ,"x"]) - resolution/2 #---because we are moving from left to right
  
  ## re-projecting the grid cell centers
  coordsHabitatGridCenterScaled <- coordsHabitatGridCenter
  coordsHabitatGridCenterScaled[ ,"y"] <- (start0.y - coordsHabitatGridCenter[ ,"y"])/resolution
  coordsHabitatGridCenterScaled[ ,"x"] <- (coordsHabitatGridCenter[ ,"x"] - start0.x)/resolution
  
  
  ## reprojecting sxy 
  if(length(dim(coordsData))==2){
    if(scaleToGrid==T){
      coordsDataScaled <- coordsData
      coordsDataScaled[ ,"y"] <- (start0.y - coordsDataScaled[ ,"y"])/resolution
      coordsDataScaled[ ,"x"] <- (coordsDataScaled[ ,"x"] - start0.x)/resolution 
    }else{
      coordsDataScaled[ ,"y"] <- start0.y - coordsData[ ,"y"] * resolution 
      coordsDataScaled[ ,"x"] <- start0.x + coordsData[ ,"x"] * resolution 
    }
    
  }
  
  ## If coordsData== 3 dimensions    
  if(length(dim(coordsData))==3){
      if(scaleToGrid==T){
        Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] 
        X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]]  
        
        coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] <- (start0.y - Y)/resolution
        coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]] <- (X - start0.x)/resolution 
      }else{
        coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]]] <- start0.y - Y * resolution 
        coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]]] <- start0.x + X * resolution 
      }
  }
  
  ## If coordsData== 4 dimensions    
  if(length(dim(coordsData))==4){
    if(scaleToGrid==T){
      Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] 
      X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]]  
      
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] <- (start0.y - Y)/resolution
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]] <- (X - start0.x)/resolution 
    }else{
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]]] <- start0.y - Y * resolution 
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]]] <- start0.x + X * resolution 
    }
  }
  ## If coordsData== 5 dimensions    
  if(length(dim(coordsData))==5){
    if(scaleToGrid==T){
      Y <- coordsData[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] 
      X <- coordsData[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]]  
      
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] <- (start0.y - Y)/resolution
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]] <- (X - start0.x)/resolution 
    }else{
      coordsDataScaled[whereYdataList[[1]], whereYdataList[[2]], whereYdataList[[3]], whereYdataList[[4]], whereYdataList[[5]]] <- start0.y - Y * resolution 
      coordsDataScaled[whereXdataList[[1]], whereXdataList[[2]], whereXdataList[[3]], whereXdataList[[4]], whereYdataList[[5]]] <- start0.x + X * resolution 
    }
  }
  

  
  return(list(coordsDataScaled = coordsDataScaled,
              coordsHabitatGridCenterScaled = coordsHabitatGridCenterScaled))
  
  
}