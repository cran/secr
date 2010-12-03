############################################################################################
## package 'secr'
## empirical.R
## design-based variance of derived density
## last changed 2010 06 20
############################################################################################

empirical.varD <- function ( object, esa = NULL, se.esa = NULL ) {
    if (inherits (object, 'secr')) {
        if (!is.list(object$capthist) )
            stop ("requires multi-session 'capthist' object")
        tempm <- object$model
        tempm$D <- NULL
        if (any (sapply(tempm, function (x) x != ~1)))
            stop ('assumes constant detection model')
        der <-  derived(object, se.esa = TRUE, se.D = FALSE)[[1]][,1:2]
        esa <- der['esa','estimate']
        se.esa <- der['esa','SE.estimate']
        nj <- sapply(object$capthist, nrow)    ## animals per session
    }
    else {
        nj <- object
        if (is.null(esa) | is.null(se.esa))
            stop("requires valid 'esa' and 'se.esa'")
    }
    n <- sum (nj)                          ## total animals
    K <- length(nj)                        ## number of sessions
    varn <- sum ((nj - n/K)^2) * K / (K-1)
    D <- n / esa / K                       ## overall density
    varD <- D^2 * (varn/n^2 + (se.esa/esa)^2)
    c(D = D, seD = varD^0.5, CVD = varD^0.5/D, CVn = varn^0.5/n, CVa = se.esa/esa )
}

## empirical.varD(trial)


