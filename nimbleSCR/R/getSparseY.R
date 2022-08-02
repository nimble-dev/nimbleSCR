#' Sparse Matrix Preparation
#'
#' R utility function to turn a two or three-dimensional detection array into a sparse matrix representation (see Turek et al., 2021 <doi.org/10.1002/ecs2.3385>  for more details).
#' Used in the implementation of the \code{\link{dbinomLocal_normal}} and \code{\link{dpoisLocal_normal}} functions.
#' 
#' The \code{getSparseY} function is used in advance of model building to create a sparse matrix representation of the observation data. 
#' It creates and returns a list of objects:
#' 
#'
#' @param x A two- or three-dimensional observation data array with dimensions : number of  individuals, number of traps, (and number of detection occasions/sessions).
#' @param noDetections The value indicating no detection. Defaults to -1.
#' @param nMaxTraps The maximum number of traps at which detections can occur. 
#' It is necessary to artificially augment the sparse detection array when using the random generation functionality of the \link{dbinomLocal_normal} or \link{dpoisLocal_normal} functions.
#' When simulating detection data, augmenting the size of the detection array is necessary to avoids artificially limiting the number of detectors at which individuals can be detected.
#' Default value is maxDetNums * 2, which doubles the maximum number of traps at which an individual can be detected. 
#' We generally recommend using \emph{numLocalIndicesMax} obtained from \code{\link{getLocalObjects}} when aiming at randomly generating detections from \link{dbinomLocal_normal} or \link{dpoisLocal_normal}.

#' 
#' @return A list of objects which constitute a sparse representation of the observation data:
#'  
#' \itemize{
#' \item \emph{detNums} A matrix with number of traps at which each individual (in rows) was detected at each occasions/sessions (in columns).
#' \item \emph{maxDetNums} The maximum number of traps at which an individual was detected (i.e., the maximum of \emph{detNums}).
#' \item \emph{detIndices} An array of dimensions n.individuals, maxDetNums, and number of occasions/sessions, which contains the IDs of the traps where each individual was detected.
#' \item \emph{y} An array of dimensions n.individuals, maxDetNums, and occasions/sessions, which contains the number of observations of each individual at the traps it was detected at.
#' \item \emph{yCombined} An array that combines \emph{detNums}, \emph{y}, and \emph{detIndices} by columns (in that specific order). 
#' Note that  \emph{y}, and \emph{detIndices} are augmented before combining, such that the maximum number of detectors at which an individual can be detected is equal to \emph{nMaxTraps}
#' Consequently, the number of columns of \emph{lengthYCombined} is 2*nMaxTraps + 1.
#' \item \emph{lengthYCombined} Dimension of the augmented lengthYCombined object to be specified as the argument \emph{lengthYCombined} of the \code{\link{dbinomLocal_normal}}  or \code{\link{dpoisLocal_normal}} functions when simulating detection data.
#' 
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
                        noDetections = -1,
                        nMaxTraps = NULL 
                        
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
  
  # IF DETECTIONS SHOULD BE BINDED TO ALLOW NUMBER LARGE NUMBER OF DETECTORS AT WHICH DETECTIONS CAN OCCUR
  if(!is.null(nMaxTraps)){
    if(max(detNums) >= nMaxTraps){
      
      nMaxTraps <-max(detNums)*2
      print("Warnings! nMaxTraps was less than or equal to maxDetNums (the maximum number of spatial recaptures). We have set nMaxTraps to 2*maxDetNums")
    }
  }else{# if nMaxTraps is not provided, the default value is max(detNums)*2
    nMaxTraps <- max(detNums)*2
  }
  
  increaseSize <- nMaxTraps 
  yCombined <- array(NA, c(dim(ySparse)[1], 1 + (increaseSize)*2, dim(ySparse)[3])) 
  
  if(dim(ySparse)[2]==1){#Deal with cases where there is a maximum of one detection per individuals
    for(t in 1:dim(yCombined)[3]){
      yCombined[,,t] <- cbind(detNums[,t],
                              ySparse[,,t],
                              matrix(-1, nrow=length(ySparse[,,t])[1], ncol=increaseSize-max(detNums)),
                              detIndices[,,t],
                              matrix(-1, nrow=length(ySparse[,,t])[1], ncol=increaseSize-max(detNums))
      )
    }
  }else{
    for(t in 1:dim(yCombined)[3]){
      yCombined[,,t] <- cbind(detNums[,t],
                              ySparse[,,t],
                              matrix(-1, nrow=dim(ySparse[,,t])[1], ncol=increaseSize-max(detNums)),
                              detIndices[,,t],
                              matrix(-1, nrow=dim(ySparse[,,t])[1], ncol=increaseSize-max(detNums))
      )
    }
  }
  
  return(list( y = ySparse,
               detIndices = detIndices,
               detNums = detNums,
               maxDetNums = max(detNums),
               yCombined = yCombined,
               lengthYCombined = dim(yCombined)[2]))
  
  
}
