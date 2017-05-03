## 2014-02-11
## 2016-06-04 annotation; memo -> message
## 2017-04-04 prefix argument

par.secr.fit <- function (arglist, ncores = 1, seed = 123, trace = TRUE,
                          logfile = "logfile.txt", prefix = "fit.") {
    ptm  <- proc.time()

    ## 'inherits' causes R to search in enclosing frames
    if (is.character(arglist))
        arglist <- mget(arglist, inherits = TRUE)
    
    ## force 'trace' to common value across all components of arglist
    arglist <- lapply(arglist, function (x) {x$trace <- trace; x})

    ## check for capthist, mask, dframe mentioned by name
    ## objects are exported to the worker processes as required
    getnames <- function(obj = 'capthist') {
        tmpnames <- sapply(arglist, function(x) if (is.character(x[[obj]])) x[[obj]] else '')
        unique(tmpnames)
    }
    data <- c(getnames('capthist'), getnames('mask'),getnames('dframe'),getnames('details'))
    data <- data[nchar(data)>0]

    ## default details savecall to FALSE across all components of arglist
    arglist <- lapply(arglist, function (x) {
        if (is.null(x$details))
            x$details <- list(savecall = FALSE)
        else if (!('savecall' %in% names(x$details))) {
            x$details[['savecall']] <- FALSE
        }
        x
    })
    
    ## individual fits must use ncores = 1
    if (ncores > 1) {
        ## force 'ncores' to 1 across all components of arglist
        arglist <- lapply(arglist, function (x) {x$ncores <- 1; x})
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big', outfile = logfile)
        clusterSetRNGStream(clust, seed)
        clusterExport(clust, c(data, 'secr.fit'), environment())
        output <- parLapply(clust, arglist, do.call, what = 'secr.fit')
        stopCluster(clust)
    }
    else {
        set.seed (seed)
        output <- lapply(arglist, do.call, what = 'secr.fit')
    }
    
    ## apply standard naming convention
    names(output) <- paste0(prefix, names(arglist))

    ## changed from memo() 2016-06-04
    message(paste('Completed in ', round((proc.time() - ptm)[3]/60,3), ' minutes at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"), sep=''))

    if (inherits(output[[1]], 'secr')) 
        secrlist(output)
    else 
        output
}

par.derived <- function (secrlist, ncores = 1, ...) {

    if (!inherits(secrlist, 'secrlist'))
        stop("requires secrlist input")
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        output <- parLapply(clust, secrlist, derived, ...)
        stopCluster(clust)
    }
    else {
        output <- lapply(secrlist, derived, ...)
    }
    names(output) <- names(secrlist)
    output
}

par.region.N <- function (secrlist, ncores = 1, ...) {

    if (!inherits(secrlist, 'secrlist'))
        stop("requires secrlist input")
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        output <- parLapply(clust, secrlist, region.N, ...)
        stopCluster(clust)
    }
    else {
        output <- lapply(secrlist, region.N, ...)
    }
    names(output) <- names(secrlist)
    output
}


