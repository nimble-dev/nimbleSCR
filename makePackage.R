
library(devtools)
library(roxygen2)
library(nimble, warn.conflicts = FALSE)

if(Sys.info()['user'] == 'dturek') {
    baseDir <- '~/github/nimble/nimbleSCR/'          ## Daniel
}
if(Sys.info()['user'] == 'pidu') {
    baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'   ## Pierre
}
if(Sys.info()['user'] == 'cymi') {
    baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'   ## Cyril
}
if(Sys.info()['user'] == 'arichbi') {
    baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'   ## Pierre
}

if(!('makePackage.R' %in% list.files(baseDir))) stop('')

document(paste0(baseDir, 'nimbleSCR'))

if(.Platform$OS.type == "windows") {
    system(paste0('R CMD build ', baseDir, 'nimbleSCR'))
} else {
    system(paste0('R CMD BUILD ', baseDir, 'nimbleSCR'))
}

check(paste0(baseDir, 'nimbleSCR'))

suppressMessages(try(remove.packages('nimbleSCR'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
if(.Platform$OS.type == "windows") {
    system(paste0('R CMD INSTALL --build ', lastTarFile))
} else {
    system(paste0('R CMD install ', lastTarFile))
}

## now quit R

## now restart R

library(nimbleSCR)

## inspect package vignettes
(vignettes <- vignette(package = 'nimbleSCR'))
for(v in vignettes$results[, 'Item'])   print(vignette(v))

dDispersal
?dDispersal

rDispersal
?rDispersal

makeGrid
?makeGrid

dbinom_vector
?dbinom_vector

dbinom_sparseLocalSCR
?dbinom_sparseLocalSCR

getSparseY
?getSparseY

getLocalTraps
?getLocalTraps

