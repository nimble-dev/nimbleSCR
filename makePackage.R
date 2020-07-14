
library(devtools)
library(roxygen2)
library(nimble, warn.conflicts = FALSE)
baseDir <- '~/github/nimble/nimbleSCR/'
# baseDir <- 'C:/Users/pidu/OneDrive - Norwegian University of Life Sciences/PROJECTS/WILDMAP/PACKAGES/nimbleSCR/'
if(!('makePackage.R' %in% list.files(baseDir))) stop('')

document(paste0(baseDir, 'nimbleSCR'))
system(paste0('R CMD BUILD ', baseDir, 'nimbleSCR'))

check(paste0(baseDir, 'nimbleSCR'))

suppressMessages(try(remove.packages('nimbleSCR'), silent = TRUE))
tarFiles <- grep('\\.tar\\.gz', list.files(baseDir, include.dirs = TRUE), value = TRUE)
lastTarFile <- tarFiles[length(tarFiles)]
message('installing package version ', gsub('\\.tar\\.gz$', '', lastTarFile))
system(paste0('R CMD install ', lastTarFile))

## now quit R

## now restart R

library(nimbleSCR)

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

