#' \code{nimble} constant function generator
#' 
#' Create global constants that can be used in \code{nimble} functions.
#' As of \code{nimble} 0.6-12 there does not appear to be official support for global constants. 
#' This function creates functions that can be used in place of constants.
#' 
#' @param inValue Scalar to be set as a constant
#' @param isDouble Logical. If \code{TRUE} the constant is a numeric scalar, and if \code{FALSE} the constant is an integer.
#' 
#' @return An object of the type returned by \code{nimbleFunction} that can be called inside
#' other \code{nimble} code to retrieve the value of the global constant.
#' 
#' @examples 
#' maxMNormTrunc <- makeConstantNimbleFunction(0.0001, TRUE)   
#' 
#' @author Joseph D. Chipperfield
#' @keywords internal
#' @export
makeConstantNimbleFunction <- function(inValue, isDouble = TRUE) {
  ## Define whether the constant is an integer or a double
  retType <- ifelse(as.logical(isDouble)[1], "double", "integer")
  ## Create a string containing the call to create a nimble function that returns the constant
  ## This will also resolve 'inValue'
  functionText <- paste(
    "nimbleFunction(run = function() {",
    paste("\treturnType(", retType, "(0))", sep = ""),
    paste("\treturn(", inValue, ")", sep = ""),
    "})", sep = "\n")
  ## Parse the string and create the nimble function
  eval(parse(text = functionText))
}
