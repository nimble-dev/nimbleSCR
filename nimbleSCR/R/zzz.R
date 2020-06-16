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
                )
            ),
            verbose = FALSE)
        
    })
}
