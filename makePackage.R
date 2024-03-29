rm(list=ls())
library(devtools)
library(roxygen2)
library(nimble, warn.conflicts = FALSE)

if(Sys.info()['user'] == 'dturek') {
    baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
    baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
    baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
    baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else if(Sys.info()['user'] == 'admin') {                    ## Soumen
    baseDir <- '~/GitHubSD/nimbleSCR/'         
} else stop('unknown user')



if(!('makePackage.R' %in% list.files(baseDir))) stop('')


if(!Sys.info()['user'] == 'admin') { 
    document(paste0(baseDir, 'nimbleSCR'))
} else if(Sys.info()['user'] == 'admin') {                    ## Soumen
    document('C:/Users/admin/OneDrive - Norwegian University of Life Sciences/Documents/GitHubSD/nimbleSCR/nimbleSCR')
}


if(.Platform$OS.type == 'windows') {
    system(paste0('R CMD build ', baseDir, 'nimbleSCR'))
} else {
    system(paste0('R CMD BUILD ', baseDir, 'nimbleSCR'))
}

if(!Sys.info()['user'] == 'admin') { 
    check(paste0(baseDir, 'nimbleSCR'))
} else if(Sys.info()['user'] == 'admin') {                    ## Soumen
    check('C:/Users/admin/OneDrive - Norwegian University of Life Sciences/Documents/GitHubSD/nimbleSCR/nimbleSCR')
}

suppressMessages(try(remove.packages('nimbleSCR'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
if(.Platform$OS.type == 'windows') {
    system(paste0('R CMD INSTALL --build ', lastTarFile))
} else {
    system(paste0('R CMD install ', lastTarFile))
}

## now quit R

## now restart R

library(nimbleSCR)


# # Run testthat tests
# if(Sys.info()['user'] == 'dturek') {
#     baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
# } else if(Sys.info()['user'] == 'pidu') {
#     baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
# } else if(Sys.info()['user'] == 'cymi') {
#     baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
# } else if(Sys.info()['user'] == 'arichbi') {
#     baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
# } else stop('unknown user')
# 
# # Utility function (fast!)
# source(file.path(baseDir,'nimbleSCR/test/testthat/testUtilityFunctions.R'))
# # distributions (slow)
# source(file.path(baseDir,'nimbleSCR/test/testthat/testDistributionFunctions.R'))
library(testthat)
##test_package('nimbleSCR')
devtools::test('nimbleSCR')

dDispersal_exp
?dDispersal_exp

rDispersal_exp
?rDispersal_exp

makeGrid
?makeGrid

dbinom_vector
?dbinom_vector

dbinomLocal_normal
?dbinomLocal_normal

getSparseY
?getSparseY

dbinomLocal_exp
?dbinomLocal_exp



## inspect package vignettes
(vignettes <- vignette(package = 'nimbleSCR'))
for(v in vignettes$results[, 'Item'])   print(vignette(v))



## time building of all vignettes:

library(nimbleSCR)
library(knitr)
library(rmarkdown)
 
setwd('~/github/nimble/nimbleSCR/nimbleSCR/vignettes')
f <- list.files()
f2 <- grep('\\.[Rr]md$', f, value = TRUE)
times <- numeric(length(f2))
names(times) <- f2
 
for(iii in 1:length(f2)) {
    message(f2[iii], ':')
    tm <- system.time(render(f2[iii]))
    times[iii] <- tm[3]
}
 
times
sum(times) / 60   ## total in minutes
