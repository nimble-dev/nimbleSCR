% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLocalTraps.R
\name{getLocalTraps}
\alias{getLocalTraps}
\title{Local Trap Identification. This function is deprecated, use \code{getLocalObjects} instead.}
\usage{
getLocalTraps(
  habitatMask,
  trapCoords,
  dmax,
  resizeFactor = 1,
  plot.check = TRUE
)
}
\arguments{
\item{habitatMask}{a binary matrix object indicating which cells are considered as suitable habitat.}

\item{trapCoords}{A matrix giving the x- and y-coordinate of each trap.}

\item{dmax}{The maximal radius from a habitat cell center within which detection probability is evaluated locally for each trap.}

\item{resizeFactor}{An aggregation factor to reduce the number of habitat cells to retrieve local traps for. Defaults to 1; no aggregation.}

\item{plot.check}{A visualization option (if TRUE); displays which traps are considered "local traps" for a randomly chosen habitat cell.}
}
\value{
This function returns a list of objects:
\itemize{
\item localTrapIndices: a matrix with number of rows equal to the reduced number of habitat grid cells (following aggregation).
Each row gives the id numbers of the local traps associated with this grid cell.
\item habitatGrid: a matrix of habitat grid cells ID corresponding to the row indices in localTrapIndices. 
\item numLocalTraps: a vector of the number of local traps for each habitat grid cell in habitatGrid.
\item numLocalTrapsMax: the maximum number of local traps for any habitat grid cell ; corresponds to the number of columns in habitatGrid.
\item resizeFactor: the aggregation factor used to reduce the number of habitat grid cells.
}
}
\description{
R utility function to identify all traps within a given radius dmax of each cell in a habitat mask.
Used in the implementation of the local evaluation approach in SCR models (\code{\link{dbinomLocal_normal}}). The distance to the activity center and the 
detection probability are then calculated for these local traps only (i.e. the detection probability is assumed to be 0 
for all other traps as they are far enough from the activity center).
}
\details{
The \code{getLocalTraps} function is used in advance of model building.
}
\examples{
colNum <- sample(20:100,1)
rowNum <- sample(20:100,1)
trapCoords <- expand.grid(list(x = seq(0.5, colNum, 1),
                               y = seq(0.5, rowNum, 1)))

habitatMask <- matrix(rbinom(colNum*rowNum, 1, 0.8), ncol = colNum, nrow = rowNum)

localTraps.list <- getLocalTraps(habitatMask, trapCoords, resizeFactor = 1, dmax = 7)

}
\author{
Cyril Milleret and Pierre Dupont
}
