.onAttach <- function(libname, pkgname) {
    
    packageStartupMessage("Warning message:\n'getLocalTraps' and 'dbinom_sparseLocalSCR' are deprecated.\nUse 'getLocalObjects' and 'dbinomLocal_normal' instead.")
    
    suppressMessages({
 
        # dDispersal_exp
        registerDistributions(
            list(
                dDispersal_exp = list(
                    BUGSdist = 'dDispersal_exp(s, rate)',
                    types = c('value = double(1)', 's = double(1)', 'rate = double()'),
                    discrete = FALSE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)

        # dHabitatMask
        registerDistributions(
            list(
                dHabitatMask = list(
                    BUGSdist = 'dHabitatMask(s, xmax, xmin, ymax, ymin, habitatMask)',
                    types = c('s = double(1)', 'habitatMask = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)

        # dbinom_sparseLocalSCR
        registerDistributions(
            list(
                dbinom_sparseLocalSCR = list(
                    BUGSdist = 'dbinom_sparseLocalSCR(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor, habitatGrid, indicator)',
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)
        
        # dbinomLocal_normal
        registerDistributions(
            list(
                dbinomLocal_normal = list(
                    BUGSdist ='dbinomLocal_normal(detNums       , detIndices    , size, p0       , p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                    Rdist = c('dbinomLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)'
                    ),
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0Traps = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = F)
        
        # dpoisLocal_normal
        registerDistributions(
            list(
                dpoisLocal_normal = list(
                    BUGSdist ='dpoisLocal_normal(detNums       , detIndices    , lambda       , lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                    Rdist = c('dpoisLocal_normal(detNums       , detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda = -999, lambdaTraps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              
                              'dpoisLocal_normal(detNums       , detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums = -999, detIndices = s, lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dpoisLocal_normal(detNums = -999, detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dpoisLocal_normal(detNums       , detIndices    , lambda       , lambdaTraps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)'
                    ),
                    types = c('value = double(1)', 'detIndices = double(1)', 'lambdaTraps = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = F)
        
        # dbinom_vector
        registerDistributions(
            list(
                dbinom_vector = list(
                    BUGSdist = 'dbinom_vector(size, prob)',
                    types = c('value = double(1)', 'size = double(1)', 'prob = double(1)'),
                    discrete = TRUE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)
        
        # dbernppAC
        registerDistributions(
            list(
                dbernppAC = list(
                    BUGSdist = "dbernppAC(lowerCoords, upperCoords, logIntensities, logSumIntensity, habitatGrid, numGridRows, numGridCols)",
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)", "logIntensities = double(1)", 
                              "logSumIntensity = double(0)", "habitatGrid = double(2)", "numGridRows = double(0)", "numGridCols = double(0)"),
                    
                    discrete = TRUE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)
        
        # dbernppACmovement_normal
        registerDistributions(
            list(
            dbernppACmovement_normal= list(
                BUGSdist = "dbernppACmovement_normal (lowerCoords, upperCoords, s, sd, baseIntensities, habitatGrid,
                                                      numGridRows, numGridCols, numWindows)",
                types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                          "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                          "habitatGrid = double(2)", "numGridRows = double(0)", "numGridCols = double(0)", "numWindows = double(0)" ),
                pqAvail = FALSE,
                mixedSizes = FALSE   
            )
          ),
        verbose = FALSE)
        
        # dbernppDetection_normal
        registerDistributions(
            list(
                dbernppDetection_normal = list(
                    BUGSdist = "dbernppDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, windowIndex,
                                                        numPoints, numWindows, indicator)",
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "windowIndex = double(0)", "numPoints = double(0)", "numWindows = double(0)", "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        # dbernppLocalACmovement_normal
        registerDistributions(
            list(
                dbernppLocalACmovement_normal = list(
                    BUGSdist = "dbernppLocalACmovement_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGrid, habitatGridLocal,
                                                               resizeFactor, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "habitatGrid = double(2)", "habitatGridLocal = double(2)",
                              "resizeFactor = double(0)","localHabWindowIndices  = double(2)", "numLocalHabWindows = double(1)",
                              "numGridRows = double(0)", "numGridCols = double(0)", "numWindows = double(0)" ),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        # dbernppLocalDetection_normal
        registerDistributions(
            list(
                dbernppLocalDetection_normal = list(
                    BUGSdist = "dbernppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, windowIndex, 
                                                             habitatGridLocal, resizeFactor, localObsWindowIndices, numLocalObsWindows, numPoints, numWindows, indicator)",
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)", "windowIndex = double(0)",
                              "habitatGridLocal = double(2)", "resizeFactor = double(0)",
                              "localObsWindowIndices = double(2)", "numLocalObsWindows = double(1)", 
                              "numPoints = double(0)", "numWindows = double(0)",
                              "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
       
        # dnormalizer
        registerDistributions(
            list(
                dnormalizer = list(
                    BUGSdist = "dnormalizer(logNormConstant)",
                    types = c("value = double(0)", "logNormConstant = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        # dpoisppAC
        registerDistributions(
            list(
                dpoisppAC = list(
                    BUGSdist = "dpoisppAC(lowerCoords, upperCoords, logIntensities, sumIntensity, habitatGrid, numGridRows, numGridCols, numPoints)",
                    types = c("value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)", "logIntensities = double(1)", 
                              "sumIntensity = double(0)", "habitatGrid = double(2)", "numGridRows = double(0)", "numGridCols = double(0)", "numPoints = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        
        # dpoisppDetection_normal
        registerDistributions(
            list(
                dpoisppDetection_normal = list(
                    BUGSdist = "dpoisppDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, 
                                                         windowIndices, numPoints, numWindows, indicator)",
                    types = c("value = double(2)", 
                              "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "windowIndices = double(1)",
                              "numPoints = double(0)", "numWindows = double(0)", "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        # dpoisppLocalDetection_normal
        registerDistributions(
            list(
                dpoisppLocalDetection_normal = list(
                    BUGSdist = "dpoisppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, windowIndices,
                                 habitatGridLocal, resizeFactor, localObsWindowIndices, numLocalObsWindows, numPoints, numWindows, indicator)",
                    types = c("value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)", "windowIndices = double(1)",  
                              "habitatGridLocal = double(2)", "resizeFactor = double(0)",
                              "localObsWindowIndices = double(2)", "numLocalObsWindows = double(1)", 
                              "numPoints = double(0)", "numWindows = double(0)",
                              "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = FALSE   
                )
            ),
            verbose = FALSE)
        
        
    })
}
