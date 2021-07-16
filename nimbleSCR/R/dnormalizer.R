#' Normalizing constant generator
#'
#' A normalizer used for normalizing \code{nimble} distributions. 
#' It is particularly useful for fitting \code{dpoisppDetection_normal} and \code{dpoisppLocalDetection_normal} models 
#' using the semi-complete data likelihood approach.
#' 
#' @name dnormalizer
#'
#' @param x Input data, which can be any scalar and will not influence the return value.
#' @param n Integer specifying the number of realisations to generate.  Only n = 1 is supported.
#' @param logNormConstant Normalizing constant on a log scale.
#' @param log Logical. If \code{TRUE} return the log normalizing constant. Otherwise return the normalizing constant.
#'
#' @return The normalizing constant.
#' @author Wei Zhang
#' 
#' @examples 
#' dnormalizer(1, log(0.5), log = TRUE)
#' dnormalizer(0, log(0.5), log = FALSE)
#' 
NULL
#' @rdname dnormalizer
#' @export
dnormalizer <- nimbleFunction(
  run = function(
    x               = double(0),
    logNormConstant = double(0),
    log             = integer(0, default = 0)
  ) {
    if(log) return(logNormConstant)
    else return(exp(logNormConstant))
    returnType(double(0))
  }
)

NULL
#' @rdname dnormalizer
#' @export
rnormalizer <- nimbleFunction(
  run = function(
    n               = integer(0),
    logNormConstant = double(0)
  ) {
    print("rnormalizer is not implemented.")
    returnType(double(0))
    return(logNormConstant)
  }
)

