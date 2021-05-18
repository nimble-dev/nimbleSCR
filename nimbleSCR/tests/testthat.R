library(testthat)
library(nimbleSCR)

if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else stop('unknown user')

## test individual files 
testthat::test_file(file.path(baseDir,"nimbleSCR/tests/testthat/test-UtilityFunctions.R"))
##
testthat::test_file(file.path(baseDir,"nimbleSCR/tests/testthat/test-dpoisLocal_normal.R"))
testthat::test_file(file.path(baseDir,"nimbleSCR/tests/testthat/test-dbinomLocal_normal.R"))

#TEST ALL FILES 
#test("nimbleSCR")
