#' Sparse Matrix Preparation
#'
#' R utility function to turn a two or three-dimensional detection array into a sparse matrix representation.
#' Used in the implementation of the \code{\link{dbinom_sparseLocalSCR}} function.
#' 
#' The \code{getSparseY} function is used in advance of model building to create a sparse matrix representation of the observation data. 
#' It creates and returns a list of objects:
#' 
#'
#' @param x A two- or three-dimensional observation data array.
#' @param noDetections The value indicating no detection. Defaults to -1.
#' 
#' @return a list of objects which constitute a sparse representation of the observation data:
#' 
#' \itemize{
#' \item \emph{detNums} a vector of the number of traps at which each individual was detected.
#' \item \emph{maxDetNums} the maximum number of traps at which any individual was detected (i.e, the maximum of detNums).
#' \item \emph{detIndices} a matrix of dimensions n.individuals * maxDetNums which contains the IDs of the traps. 
#' where each individuals was detected.
#' \item \emph{y} a matrix of dimensions n.individuals * maxDetNums which contains the number of observations of each individual
#' at the traps it was detected at.
#' }
#'
#' @author Cyril Milleret
#'
#' @examples
#' y.full <- matrix(rbinom(5000, 5, 0.02), ncol = 100)
#' 
#' y <- getSparseY(y.full)
#' 
#' @export
getSparseY <- function( x,
                        noDetections = -1
){
  # IF 2D OBSERVATION MATRIX, CONVERT TO 3D ARRAY
  if(length(dim(x))==2){
    Y <- array(x, c(dim(x),1))
  }else{Y <- x}
  
  # RETRIEVE THE NUMBER OF DETECTIONS FOR EACH ID
  detNums <- apply(Y, c(1,3), function(x) length(which(x>0)))
  
  # INITIALIAZE EMPTY ARRAYS FOR THE DETECTIONS AND DETECTOR INDICES
  detIndices <- array(-1, c(dim(Y)[1],max(detNums), dim(Y)[3]))
  ySparse <- array(-1, c(dim(Y)[1],max(detNums), dim(Y)[3]))
  
  # FILL IN THE ARRAYS
  for(t in 1:dim(Y)[3]){
    for(i in 1:dim(Y)[1]){
      if(detNums[i,t] > 0){
        # GET WHERE (DETECTOR ID) DETECTIONS OCCUR
        detIndices[i,1:detNums[i,t],t] <- which(Y[i, ,t] > 0)
        # GET NUMBER OF DETECTIONS 
        ySparse[i,1:detNums[i,t],t] <- Y[i,which(Y[i, ,t] > 0),t]
      }
    }
  }   
  
  return(list( y = ySparse,
               detIndices = detIndices,
               detNums = detNums,
               maxDetNums = max(detNums)))
}
