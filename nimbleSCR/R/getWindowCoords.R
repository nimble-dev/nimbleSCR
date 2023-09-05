#' Get lower and upper windows coordinates 
#'
#' The \code{getWindowCoords} is an R utility function to create lower and upper habitat and observation windows coordinates, as well an habitat grid with cell ids. 
#' Those objects are necessary to run all point process (pp) functions. All input data should be scaled to the habitat grid using \code{\link{scaleCoordsToHabitatGrid}}. 
#' Note that we assume homogeneous window sizes.
#'
#' @param scaledHabGridCenter A matrix with the scaled "x" and "y" habitat windows grid cell centers (after using \code{\link{scaleCoordsToHabitatGrid}}).
#' @param scaledObsGridCenter A matrix with the scaled "x" and "y" observation windows grid cell centers (afer using \code{\link{scaleCoordsToHabitatGrid}}). 
#' This is an optional argument and only necessary when modelling detection as a point process (e.g. \code{\link{dpoisppDetection_normal}}). 
#' @param plot.check A visualization option (if TRUE); displays habitat and detection windows.
#' 
#'
#' 
#' @return A list of objects :
#'  
#' \itemize{
#' \item \emph{lowerHabCoords} A matrix with the "x" and "y" lower habitat window coordinates.
#' \item \emph{upperHabCoords} A matrix with the "x" and "y" upper habitat window coordinates.
#' \item \emph{habitatGrid} A matrix of habitat cell ID that can be used to lookup efficiently the cell ID from a coordinate scaled to the habitat grid:
#' habitatGrid[trunc(scaledHabGridCenter[1,"y"]) + 1, trunc(scaledHabGridCenter[1,"x"]) + 1]. See \code{\link{scaleCoordsToHabitatGrid}} for more details.
#' \item \emph{lowerObsCoords} A matrix with the "x" and "y" lower observation window coordinates. Only returned when \emph{scaledObsGridCenter} is provided.
#' \item \emph{upperObsCoords} A matrix with the "x" and "y" upper observation window coordinates. Only returned when \emph{scaledObsGridCenter} is provided.
#'
#'
#' }
#'
#' @author Cyril Milleret
#' @importFrom grDevices adjustcolor
#' @importFrom graphics rect
#' @examples
#' coordsGridCenter <- expand.grid(list(x = seq(50.5, 100, 1),
#'                                      y = seq(100.5, 150, 1)))
#' coordsData <- expand.grid(list(x = seq(60, 90, 1),
#'                                y = seq(110, 140, 1)))
#' 
#' plot(coordsGridCenter[,2] ~ coordsGridCenter[,1])
#' points(coordsData[,2] ~ coordsData[,1], col="red")
#' scaled <- scaleCoordsToHabitatGrid(coordsData = coordsData
#'                                    , coordsHabitatGridCenter = coordsGridCenter)
#' plot(scaled$coordsHabitatGridCenterScaled[,2] ~ scaled$coordsHabitatGridCenterScaled[,1])
#' points(scaled$coordsDataScaled[,2] ~ scaled$coordsDataScaled[,1], col="red")
#' 
#' LowerAndUpperCoords <- getWindowCoords(scaledHabGridCenter = scaled$coordsHabitatGridCenterScaled,
#'                                          scaledObsGridCenter = scaled$coordsDataScaled)
#' 
#' # Plot habitat window cell centers and lower/upper coordinates 
#' plot(scaled$coordsHabitatGridCenterScaled[,2] ~ 
#'      scaled$coordsHabitatGridCenterScaled[,1], 
#'      pch=16, cex=0.3, col=grey(0.5))
#' points(LowerAndUpperCoords$lowerHabCoords[,2] ~ 
#'        LowerAndUpperCoords$lowerHabCoords[,1],
#'        pch=16, cex=0.3, col=grey(0.1))
#' points(LowerAndUpperCoords$upperHabCoords[,2] ~
#'        LowerAndUpperCoords$upperHabCoords[,1],
#'        pch=16, cex=0.3, col=grey(0.1))
#' 
#' # Plot observation window cells center and lower/upper coordinates
#' points(scaled$coordsDataScaled[,2]~scaled$coordsDataScaled[,1], pch=16,
#'  cex=0.3, col = adjustcolor("red",alpha.f = 0.8))
#' points(LowerAndUpperCoords$lowerObsCoords[,2] ~ 
#'         LowerAndUpperCoords$lowerObsCoords[,1],
#'        pch=16, cex=0.3, col = adjustcolor("red", alpha.f = 0.8))
#' points(LowerAndUpperCoords$upperObsCoords[,2] ~
#'        LowerAndUpperCoords$upperObsCoords[,1],
#'        pch=16, cex=0.3, col = adjustcolor("red", alpha.f = 0.8))
#'        
#' @export
getWindowCoords <- function(scaledHabGridCenter = scaledHabGridCenter,
                              scaledObsGridCenter = NULL,
                              plot.check = TRUE
                              
){
  ## GET LOWER AND UPPER COORDS 
  lowerHabCoords <- scaledHabGridCenter - 0.5
  upperHabCoords <- scaledHabGridCenter + 0.5
  
  ## CONSTRUCT GRID CELL
  habitatGrid <- matrix(0, nrow = max(scaledHabGridCenter[,"y"])+1,
                        ncol = max(scaledHabGridCenter[,"x"])+1 )
  
  for(i in 1:nrow(scaledHabGridCenter)){
    habitatGrid[dim(habitatGrid)[1] -trunc(scaledHabGridCenter[i,"y"]) , trunc(scaledHabGridCenter[i,"x"]) + 1] <- i
  }
  
  ## CONSTRUCT LOWER AND UPPER OBSGRID 
  ## ASSUME THAT OBS WINDOW SIZE IS HOMOGENEOUS
  if(!is.null(scaledObsGridCenter)){
    # GET RESOLUTION HABITAT GRID 
    minResX <- abs(diff(scaledObsGridCenter[,"x"]))
    minResX[minResX==0] <- NA
    minResY <- abs(diff(scaledObsGridCenter[,"y"]))
    minResY[minResY==0] <- NA
    resY <- min(minResY, na.rm = T)
    resX <- min(minResX, na.rm = T)
    
    # CONSTRUCT LOWER AND OBS COORDS 
    lowerObsCoords <- scaledObsGridCenter - resX/2
    upperObsCoords <- scaledObsGridCenter + resY/2
    
    if(plot.check){
      xrange <- range(c(lowerHabCoords[,1], upperHabCoords[,1]))
      yrange <- range(c(lowerHabCoords[,2], upperHabCoords[,2]))
      # just add some extent for plotting
      xrange[1] <- xrange[1] - 1
      xrange[2] <- xrange[2] + 1
      yrange[1] <- yrange[1] - 1
      yrange[2] <- yrange[2] + 1
      
      #habitat windows 
      plot(lowerHabCoords[,2] ~ lowerHabCoords[,1], xlim=xrange, ylim=yrange, pch=21, cex=0.6, bg="red",
           ylab= "y", xlab= "x")
      rect(xleft = lowerHabCoords[,1] ,
           ybottom = lowerHabCoords[,2] ,
           xright = upperHabCoords[,1],ytop = upperHabCoords[,2],
           col=adjustcolor("red",alpha.f = 0.4),
           border="red")
      points(lowerHabCoords[,2] ~ lowerHabCoords[,1], xlim=xrange, ylim=yrange, pch=21,cex=0.6, bg="red")
      points(upperHabCoords[,2] ~ upperHabCoords[,1], xlim=xrange, ylim=yrange, pch=21,cex=0.6, bg="red")
      points(scaledHabGridCenter[,2] ~ scaledHabGridCenter[,1], xlim=xrange, ylim=yrange,
             pch=16, cex=0.3, col="red")
      
      
      #Observation windows 
      points(lowerObsCoords[,2] ~ lowerObsCoords[,1], xlim=xrange, ylim=yrange,
             pch=21, cex=0.6, bg="blue")
      rect(xleft = lowerObsCoords[,1] ,
           ybottom = lowerObsCoords[,2] ,
           xright = upperObsCoords[,1],ytop = upperObsCoords[,2],
           col=adjustcolor("blue",alpha.f = 0.1),
           border="blue")
      points(lowerObsCoords[,2] ~ lowerObsCoords[,1], xlim=xrange, ylim=yrange,
             pch=21,cex=0.6, bg="blue")
      points(upperObsCoords[,2] ~ upperObsCoords[,1], xlim=xrange, ylim=yrange,
             pch=21,cex=0.6, bg="blue")
      points(scaledObsGridCenter[,2] ~ scaledObsGridCenter[,1], xlim=xrange, ylim=yrange,
             pch=16, cex=0.3, col="blue")
      
      
    }
    
    return(list(lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid,
                lowerObsCoords=lowerObsCoords,
                upperObsCoords=upperObsCoords
    ))
    
  }else{
    
    
    if(plot.check){
      xrange <- range(c(lowerHabCoords[,1], upperHabCoords[,1]))
      yrange <- range(c(lowerHabCoords[,2], upperHabCoords[,2]))
      # just add some extent for plotting
      xrange[1] <- xrange[1] - 1
      xrange[2] <- xrange[2] + 1
      yrange[1] <- yrange[1] - 1
      yrange[2] <- yrange[2] + 1
      
      #habitat windows 
      plot(lowerHabCoords[,2] ~ lowerHabCoords[,1], xlim=xrange, ylim=yrange, pch=21, cex=0.6, bg="red",
           ylab= "y", xlab= "x")
      rect(xleft = lowerHabCoords[,1] ,
           ybottom = lowerHabCoords[,2] ,
           xright = upperHabCoords[,1],ytop = upperHabCoords[,2],
           col=adjustcolor("red",alpha.f = 0.4),
           border="red")
      points(lowerHabCoords[,2] ~ lowerHabCoords[,1], xlim=xrange, ylim=yrange, pch=21,cex=0.6, bg="red")
      points(upperHabCoords[,2] ~ upperHabCoords[,1], xlim=xrange, ylim=yrange, pch=21,cex=0.6, bg="red")
      points(scaledHabGridCenter[,2] ~ scaledHabGridCenter[,1], xlim=xrange, ylim=yrange,
             pch=16, cex=0.3, col="red")
      
      rect(xleft = lowerHabCoords[,1] ,
           ybottom = lowerHabCoords[,2] ,
           xright = upperHabCoords[,1],ytop = upperHabCoords[,2],
           col=adjustcolor("red",alpha.f = 0.4),
           border="red")
      
    }
    
    return(list(lowerHabCoords = lowerHabCoords,
                upperHabCoords = upperHabCoords,
                habitatGrid = habitatGrid
    ))
    
  }
  
  
  
  
}
