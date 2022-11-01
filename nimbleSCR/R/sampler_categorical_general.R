#' @title \code{nimble} MCMC sampler function for general categorial distributions
#'
#' @description
#' The \code{categorical_general} sampler operates within \code{nimble}'s MCMC engine to perform Gibbs sampling for a single node, which must in essence follow a categorical distribution.  However, the prior distribution need not be \code{nimble}'s \code{dcat} distribution, but rather can be any (potentially user-defined) distribution which has the same support as a standard categorical (\code{dcat}) distribution.  Specifically: the distribution must define a discrete random variable, which can only attain values from the set {1, 2, 3, ..., \code{numCategories}}.
#'
#' The \code{categorical_general} sampler requires one control list argument, named \code{numCategories}, which specifies the fixed upper-bound for the range of the random variable.
#'
#' The \code{categorical_general} sampler is designed to be used in \code{nimble}'s MCMC engine, and can be added to an MCMC configuration object using the \code{addSampler} method.  See \code{help(configureMCMC)} for more information about MCMC configuration objects and adding custom samplers.
#' 
#' @param model (uncompiled) model on which the MCMC is to be run. 
#' @param mvSaved \code{modelValues} object to be used to store MCMC samples. 
#' @param target node on which the sampler will operate.
#' @param control named list containing an elemented named \code{numCategories}, which specifies the upper-bound for the range of the random variable.
#'
#' @author Daniel Turek
#'
#' @examples
#'
#' \dontrun{
#' ## define custom dmy_categorical distribution as a nimbleFunction
#' dmy_categorical <- nimbleFunction(...)
#'
#' ## nimble model code, using custom-written dmy_categorical distribution
#' code <- nimbleCode({
#'   x ~ dmy_categorical(...)
#' })
#'
#' ## create NIMBLE model object
#' Rmodel <- nimbleModel(code)
#'
#' ## create MCMC configuration object with no samplers
#' conf <- configureMCMC(Rmodel, nodes = NULL)
#'
#' ## add categorical_general sampler to MCMC configuration
#' conf$addSampler(target = 'x', type = 'categorical_general', control = list(numCategories = 10))
#'
#' ## build MCMC algorithm
#' Rmcmc <- buildMCMC(conf)
#'
#' ## compile model and MCMC, run MCMC algorithm
#' }
#'
#' @export
sampler_categorical_general <- nimbleFunction(
  name = 'sampler_categorical_general',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## control list extraction
    k <- extractControlElement(control, 'numCategories', error = 'categorical_general sampler missing required control argument: numCategories')
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
    calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    probs <- numeric(k)
    logProbs <- numeric(k)
    ## checks
    if(length(targetAsScalar) > 1)   stop('cannot use categorical_general sampler on more than one target node')
  },
  run = function() {
    currentValue <- model[[target]]
    logProbs[currentValue] <<- getLogProb(model, calcNodes)
    for(i in 1:k) {
      if(i != currentValue) {
        model[[target]] <<- i
        logProbPrior <- calculate(model, target)
        if(logProbPrior == -Inf) {
          logProbs[i] <<- -Inf
        } else {
          if(is.nan(logProbPrior)) {
            logProbs[i] <<- -Inf
          } else {
            logProbs[i] <<- logProbPrior + calculate(model, calcNodesNoSelf)
            if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
          }
        }
      }
    }
    logProbs <<- logProbs - max(logProbs)
    probs <<- exp(logProbs)
    newValue <- rcat(1, probs)   ## rcat normalizes the probabilitiess internally
    if(newValue != currentValue) {
      model[[target]] <<- newValue
      model$calculate(calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    } else {
      nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
      nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
    }
  },
  methods = list(
    reset = function() { }
  )
)


