############################################################################################
## package 'secr'
## esa.R
## last changed 2009 06 29 2009 10 24
## 2010 02 25 dropped details$spherical
## 2010 02 26 fixed for nmix>1
## 2010 03 04 use scaled.detection from functions.R
## 2010 03 09 fixed bug : need to call scaled.detection when !is.null(real)
## 2011-04-04 added noccasions; debugged 2011-04-07
############################################################################################

esa <- function (object, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL)

# Return vector of 'a' for given g0, sigma, [z (if hazard fn) ] and session
# detectfn is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, mask, detectfn

## strictly doesn't need data, so better to avoid when object not available...
{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (ms(object$mask))
        mask <- object$mask[[sessnum]]
    else
        mask <- object$mask

    if (is.null(beta) & is.null(real))
        beta <- object$fit$par

    traps   <- attr(capthists, 'traps')  ## need session-specific traps
    if (!(detector(traps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for esa")
    dettype <- detectorcode(traps)
    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)
    constant <- !is.null(noccasions)    ## fix 2011-04-07
    if (is.null(noccasions)) {
        noccasions <- s
    }

    nmix    <- object$details$nmix
    if (is.null(nmix))
        nmix <- 1

    ##############################################
    ## adapt for marking occasions only 2009 10 24
    q <- attr(capthists, 'q')
    if (!is.null(q))
        if (q<s) s <- q
    ##############################################

    if (dettype %in% c(3,6)) {
        k <- c(table(polyID(traps)),0)
        K <- length(k)-1
    }
    else if (dettype %in% c(4,7)) {
        k <- c(table(transectID(traps)),0)
        K <- length(k)-1
    }
    else {
        k <- nrow(traps)
        K <- k
    }
    m      <- length(mask$x)            ## need session-specific mask...
    cell   <- attr(mask,'area')

    if (constant) {
        ## assume constant
        if (is.null(beta))
            real <- detectpar(object)
        else {
            real <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            real <- as.list(real)
            names(real) <- parnames(object$detectfn)
        }
        a <- cell * sum(pdot(X = mask, traps = traps, detectfn = object$detectfn,
                             detectpar = real, noccasions = noccasions))
        return(rep(a,n))
    }
    else {
        if (is.null(beta)) {
            if (is.null(real))
                stop ("requires real parameter values")
            PIA <- rep(1, n * s * K * nmix)    ## nmix added 2010 02 25

            ## new code 2010-11-26
            realparval0 <- matrix(rep(real, rep(n,length(real))), nrow = n)   ## UNTRANSFORMED
            ## replacing...
            ##        realparval0 <- matrix(real, nrow = 1)   ## UNTRANSFORMED
            ##        ncolPIA <- 1
            Xrealparval0 <- scaled.detection (realparval0, FALSE, object$details$scaleg0, NA)
        }
        else {
            ## allow for old design object
            if (length(dim(object$design0$PIA))==4)
                dim(object$design0$PIA) <- c(dim(object$design0$PIA),1)
            PIA <- object$design0$PIA[sessnum,,1:s,,,drop=F]
            ncolPIA <- dim(object$design0$PIA)[2]

            #############################################
            ## trick to allow for changed data 2009 11 20
            ## nmix>1 needs further testing 2010 02 26
            ## NOTE 2010-11-26 THIS LOOKS WEAK
            if (dim(PIA)[2] != n) {
                PIA <- array(rep(PIA[1,1,,,],n), dim=c(s,K,nmix,n))
                PIA <- aperm(PIA, c(4,1,2,3))   ## n,s,K,nmix
                ncolPIA <- n     ## 2010 02 26
            }
            #############################################

            realparval0 <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive

            Xrealparval0 <- scaled.detection (realparval0, FALSE,
                object$details$scaleg0, NA)
        }

        ## new code 2010-11-26
        used <- usage(traps)
        if (any(used==0))
        PIA <- PIA * rep(rep(t(used),rep(n,s*K)),nmix)
        ncolPIA <- n

        param <- object$details$param
        if (is.null(param))
            param <- 0    ## default Borchers & Efford (vs Gardner & Royle)
        gamma <- 1  ## DUMMY
        useD <- FALSE
        temp <- .C("integralprw1", PACKAGE = 'secr',
            as.integer(dettype),
            as.integer(param),
            as.double(Xrealparval0),
            as.integer(n),
            as.integer(s),
            as.integer(k),
            as.integer(m),
            as.integer(nmix),
            as.double(unlist(traps)),
            as.double(unlist(mask)),
            as.integer(nrow(Xrealparval0)), # rows in lookup
            as.integer(PIA),                # index of nc*,S,K to rows in realparval0
            as.integer(ncolPIA),            # ncol - if CL, ncolPIA = n, else ncolPIA = 1 or ngrp
            as.double(cell),
            as.double(gamma),
            as.integer(object$detectfn),
            as.integer(object$details$binomN),
            as.double(object$details$cutval),
            as.integer(useD),
            a=double(n),
            resultcode=integer(1)
       )
       if (temp$resultcode == 3)
           stop ("groups not implemented in external function 'integralprw1'")
       if (temp$resultcode != 0)
           stop ("error in external function 'integralprw1'")
       return(temp$a)
    }
}
############################################################################################
