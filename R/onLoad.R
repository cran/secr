###############################################################################
## package 'secr'
## onLoad.R
## 2020-02-21
###############################################################################

.onLoad <- function (libname, pkgname) {
    ## also sets environment variable RCPP_PARALLEL_NUM_THREADS
    RcppParallel::setThreadOptions(min(2, RcppParallel::defaultNumThreads()))  
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package