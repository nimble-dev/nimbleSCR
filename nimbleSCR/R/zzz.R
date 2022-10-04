.onAttach <- function(libname, pkgname) {
    
    
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
        
        
        # dmultiLocal_normal
        registerDistributions(
          list(
            dmultiLocal_normal = list(
              BUGSdist ='dmultiLocal_normal(detNums       , detIndices    , size, p0       , p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
              Rdist = c('dmultiLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                        
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                        'dmultiLocal_normal(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                        'dmultiLocal_normal(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)'
              ),
              types = c('value = double(1)', 'detIndices = double(1)', 'size = double(0)', 'p0Traps = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
              discrete = TRUE,
              mixedSizes = TRUE,
              pqAvail = FALSE
            )
          ),
          verbose = F)
        
        # dbinomLocal_normalPlateau
        registerDistributions(
            list(
                dbinomLocal_normalPlateau = list(
                    BUGSdist ='dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0       , p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                    Rdist = c('dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_normalPlateau(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_normalPlateau(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, w, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)'
                    ),
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0Traps = double(1)', 'w = double(0)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = F)
        
        # dbinomLocal_exp
        registerDistributions(
            list(
                dbinomLocal_exp = list(
                    BUGSdist ='dbinomLocal_exp(detNums       , detIndices    , size, p0       , p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                    Rdist = c('dbinomLocal_exp(detNums       , detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0 = -999, p0Traps    , rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums = -999, detIndices = s, size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined)',
                              'dbinomLocal_exp(detNums = -999, detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)',
                              'dbinomLocal_exp(detNums       , detIndices    , size, p0       , p0Traps = s, rate, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0)'
                    ),
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0Traps = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
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
                    mixedSizes = TRUE,
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
                mixedSizes = TRUE
            )
          ),
        verbose = FALSE)
        
        # dbernppACmovement_exp
        registerDistributions(
          list(
            dbernppACmovement_exp= list(
              BUGSdist = "dbernppACmovement_exp(lowerCoords, upperCoords, s, lambda, rate, baseIntensities, habitatGrid, numGridRows, numGridCols, numWindows)",
              Rdist = c("dbernppACmovement_exp(lowerCoords, upperCoords, s, lambda = -999, rate, baseIntensities, habitatGrid, numGridRows, numGridCols, numWindows)",
                        "dbernppACmovement_exp(lowerCoords, upperCoords, s, lambda, rate = -999, baseIntensities, habitatGrid, numGridRows, numGridCols, numWindows)"),
              types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                        "s = double(1)", "lambda = double(0)", "rate = double(0)", "baseIntensities = double(1)",
                        "habitatGrid = double(2)", "numGridRows = double(0)", "numGridCols = double(0)", "numWindows = double(0)"),
              discrete = TRUE,
              mixedSizes = TRUE,
              pqAvail = FALSE
            )
          ),
          verbose = FALSE)
        
        # dbernppDetection_normal
        registerDistributions(
            list(
                dbernppDetection_normal = list(
                    BUGSdist = "dbernppDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities,  numWindows, indicator)",
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "numWindows = double(0)", "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = TRUE)
            ),
            verbose = FALSE)
        
        # dbernppLocalACmovement_normal
        registerDistributions(
            list(
                dbernppLocalACmovement_normal = list(
                    BUGSdist = "dbernppLocalACmovement_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGrid, habitatGridLocal,
                                                               resizeFactor, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                    Rdist = c("dbernppLocalACmovement_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGrid, habitatGridLocal, resizeFactor    , localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                              "dbernppLocalACmovement_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGrid, habitatGridLocal, resizeFactor = 1, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)"),
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "habitatGrid = double(2)", "habitatGridLocal = double(2)",
                              "resizeFactor = double(0)","localHabWindowIndices  = double(2)", "numLocalHabWindows = double(1)",
                              "numGridRows = double(0)", "numGridCols = double(0)", "numWindows = double(0)" ),
                    pqAvail = FALSE,
                    mixedSizes = TRUE
                )
            ),
            verbose = FALSE)
        
        # dbernppLocalACmovement_exp
        registerDistributions(
          list(
            dbernppLocalACmovement_exp = list(
              BUGSdist = "dbernppLocalACmovement_exp(lowerCoords, upperCoords, s, lambda, rate, baseIntensities, habitatGrid, habitatGridLocal, resizeFactor, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
              Rdist = c("dbernppLocalACmovement_exp(lowerCoords, upperCoords, s, lambda = -999, rate,        baseIntensities, habitatGrid, habitatGridLocal, resizeFactor    , localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                        "dbernppLocalACmovement_exp(lowerCoords, upperCoords, s, lambda       , rate = -999, baseIntensities, habitatGrid, habitatGridLocal, resizeFactor    , localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                        "dbernppLocalACmovement_exp(lowerCoords, upperCoords, s, lambda = -999, rate,        baseIntensities, habitatGrid, habitatGridLocal, resizeFactor = 1, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)",
                        "dbernppLocalACmovement_exp(lowerCoords, upperCoords, s, lambda       , rate = -999, baseIntensities, habitatGrid, habitatGridLocal, resizeFactor = 1, localHabWindowIndices, numLocalHabWindows, numGridRows, numGridCols, numWindows)"),
              types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                        "s = double(1)", "lambda = double(0)", "rate = double(0)", "baseIntensities = double(1)",
                        "habitatGrid = double(2)", "habitatGridLocal = double(2)", "resizeFactor = double(0)",
                        "localHabWindowIndices = double(2)", "numLocalHabWindows = double(1)",
                        "numGridRows = double(0)", "numGridCols = double(0)", "numWindows = double(0)"),
              pqAvail = FALSE,
              mixedSizes = TRUE
            )
          ),
          verbose = FALSE)
        
        # dbernppLocalDetection_normal
        registerDistributions(
            list(
                dbernppLocalDetection_normal = list(
                    BUGSdist = "dbernppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, 
                                 habitatGridLocal, resizeFactor, localObsWindowIndices, numLocalObsWindows,  numWindows, indicator)",
                    Rdist = c("dbernppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGridLocal, resizeFactor    , localObsWindowIndices, numLocalObsWindows,  numWindows, indicator)",
                              "dbernppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGridLocal, resizeFactor = 1, localObsWindowIndices, numLocalObsWindows,  numWindows, indicator)"),
                    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)", 
                              "habitatGridLocal = double(2)", "resizeFactor = double(0)",
                              "localObsWindowIndices = double(2)", "numLocalObsWindows = double(1)", 
                              "numWindows = double(0)",
                              "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = TRUE)
            ),
            verbose = FALSE)
       
        # dnormalizer
        registerDistributions(
            list(
                dnormalizer = list(
                    BUGSdist = "dnormalizer(logNormConstant)",
                    types = c("value = double(0)", "logNormConstant = double(0)"),
                    pqAvail = FALSE
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
                    mixedSizes = TRUE
                )
            ),
            verbose = FALSE)
        
        # dpoisppDetection_normal
        registerDistributions(
            list(
                dpoisppDetection_normal = list(
                    BUGSdist = "dpoisppDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, 
                                numMaxPoints, numWindows, indicator)",
                    types = c("value = double(2)", 
                              "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",
                              "numMaxPoints = double(0)",
                              "numWindows = double(0)", "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = TRUE)
            ),
            verbose = FALSE)
        
        # dpoisppLocalDetection_normal
        registerDistributions(
            list(
                dpoisppLocalDetection_normal = list(
                    BUGSdist = "dpoisppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, 
                                 habitatGridLocal, resizeFactor, localObsWindowIndices, numLocalObsWindows, numMaxPoints, numWindows, indicator)",
                    Rdist = c("dpoisppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGridLocal, resizeFactor    , localObsWindowIndices, numLocalObsWindows, numMaxPoints, numWindows, indicator)",
                              "dpoisppLocalDetection_normal(lowerCoords, upperCoords, s, sd, baseIntensities, habitatGridLocal, resizeFactor = 1, localObsWindowIndices, numLocalObsWindows, numMaxPoints, numWindows, indicator)"),
                    types = c("value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)",
                              "s = double(1)", "sd = double(0)", "baseIntensities = double(1)",  
                              "habitatGridLocal = double(2)", "resizeFactor = double(0)",
                              "localObsWindowIndices = double(2)", "numLocalObsWindows = double(1)", 
                              "numMaxPoints = double(0)", "numWindows = double(0)",
                              "indicator = double(0)"),
                    pqAvail = FALSE,
                    mixedSizes = TRUE)
            ),
            verbose = FALSE)
        
        
    })
}
