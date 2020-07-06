.onAttach <- function(libname, pkgname) {
    suppressMessages({

        registerDistributions(
            list(
                dDispersal = list(
                    BUGSdist = 'dDispersal(S, lam)',
                    types = c('value = double(1)', 'S = double(1)', 'lam = double()'),
                    discrete = FALSE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                ),
                dhabitatMask = list(
                  BUGSdist = 'dhabitatMask(s, xmax, xmin, ymax, ymin, habitatMask)',
                  types = c('value = double()', 's = double(1)', 'xmax = double()', 'xmin = double()', 'ymax = double()', 'ymin = double()'. 'habitatMask = double(2'),
                  discrete = FALSE,
                  mixedSizes = FALSE,
                  pqAvail = FALSE
                ),
                dbinom_vector = list(
                  BUGSdist = 'dbinom_vector(size, prob)',
                  types = c('value = double(1)', 'size = double(1)', 'prob = double(1)'),
                  discrete = FALSE,
                  mixedSizes = FALSE,
                  pqAvail = FALSE
                ),
                dbinom_sparseLocalSCR = list(
                  BUGSdist = 'dbinom_sparseLocalSCR(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor, habitatGrid, indicator)',
                  types = c('x = double(1)', 'detNums = double()', 'detIndices = double(1)', 'size = double(1)', 'p0 = double()', 'sigma = double()', 's = double(1)', 'trapCoords = double(2)',
                            'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'resizeFactor = double()', 'habitatGrid = double(2)', 'indicator = double()'),
                  discrete = FALSE,
                  mixedSizes = TRUE,
                  pqAvail = FALSE
                ),
            ),
            verbose = FALSE)
    })
}
