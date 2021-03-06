% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpoisppDetection_normal.R
\name{dpoisppDetection_normal}
\alias{dpoisppDetection_normal}
\alias{rpoisppDetection_normal}
\title{Poisson point process detection model}
\usage{
dpoisppDetection_normal(
  x,
  lowerCoords,
  upperCoords,
  s,
  sd,
  baseIntensities,
  windowIndices,
  numPoints,
  numWindows,
  indicator,
  log = 0
)

rpoisppDetection_normal(
  n,
  lowerCoords,
  upperCoords,
  s,
  sd,
  baseIntensities,
  windowIndices,
  numPoints,
  numWindows,
  indicator
)
}
\arguments{
\item{x}{Matrix of x- and y-coordinates of a set of spatial points (detection locations). One row corresponds to one point.}

\item{lowerCoords, upperCoords}{Matrices of lower and upper x- and y-coordinates of all detection windows. One row for each window.}

\item{s}{Vector of x- and y-coordinates of the isotropic multivariate normal distribution mean (the AC location).}

\item{sd}{Standard deviation of the isotropic multivariate normal distribution.}

\item{baseIntensities}{Vector of baseline detection intensities for all detection windows.}

\item{windowIndices}{Vector of indices of detection windows where each point of \code{x} falls.}

\item{numPoints}{Number of points that should be considered. For \code{dpoisppDetection_normal}, this should be a non-negative integer and is used to truncate \code{x} 
so that extra rows beyond \code{numPoints} are ignored. For \code{rpoisppDetection_normal}, this can be any integer. 
Usually it is set to be a negative value, which indicates that the number of points is generated from a Poisson model. 
Otherwise, the number of points generated will be equal to the specified value.}

\item{numWindows}{Number of detection windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
so that extra rows beyond \code{numWindows} are ignored.}

\item{indicator}{Binary variable (0 or 1) used for data augmentation. \code{indicator = 0} means the individual does not exist 
and thus the probability of no detection is 1.}

\item{log}{Logical argument, specifying whether to return the log-probability of the distribution.}

\item{n}{Integer specifying the number of realisations to generate.  Only n = 1 is supported.}
}
\value{
\code{dpoisppDetection_normal} gives the (log) probability density of the observation matrix \code{x}.
\code{rpoisppDetection_normal} gives coordinates of a set of randomly generated spatial points.
}
\description{
Density and random generation functions of the Poisson point process for detection. 
An isotropic multivariate normal distribution is used as the decay kernel.
}
\examples{
x <- matrix(c(0.5, 0.2, 1.1, 0.4, 1.4, 0.3, 0.1, 1.3, 1, 1.5), nrow = 5, byrow = TRUE)
lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
s <- c(1, 1)
sd <- 0.1
baseIntensities <- c(1:4)
windowIndices <- c(1, 2, 2, 3, 4)
numPoints <- 5
numWindows <- 4
indicator <- 1
dpoisppDetection_normal(x, lowerCoords, upperCoords, s, sd, baseIntensities,
                        windowIndices, numPoints, numWindows, indicator, log = TRUE)
}
\references{
W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035
}
\author{
Wei Zhang
}
