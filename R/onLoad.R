###############################################################################
## package 'secr'
## onLoad.R
## 2020-02-21, 2021-03-13, 2021-04-13
###############################################################################

.onLoad <- function (libname, pkgname) {
    ## also sets environment variable RCPP_PARALLEL_NUM_THREADS
    RcppParallel::setThreadOptions(min(2, RcppParallel::defaultNumThreads()))
    
    ## following advice of Kevin Ushey 2020-03-18
    ## to avoid UBSAN errors from parallelFor
    ## Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
    
    ## possible improvement once RcppParallel updated 2021-03-14
    ## see https://github.com/RcppCore/RcppParallel/pull/133?_pjax=%23js-repo-pjax-container
    ## Sys.setenv(RCPP_PARALLEL_NUM_THREADS = 2)
    
}

## .onLoad is preferred if actions are required for single functions 
## that may be called without attaching package