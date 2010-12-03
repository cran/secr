############################################################################################
## package 'secr'
## esa.R
## last changed 2009 06 29 2009 10 24
## 2010 02 25 dropped details$spherical
## 2010 02 26 fixed for nmix>1
## 2010 03 04 use scaled.detection from functions.R
## 2010 03 09 fixed bug : need to call scaled.detection when !is.null(real)
############################################################################################

esa <- function (object, sessnum = 1, beta = NULL, real = NULL)

# Return vector of 'a' for given g0, sigma, [z (if hazard fn) ] and session
# detectfn is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, mask, detectfn

## strictly doesn't need data, so better to avoid when object not available...
{
    if (inherits(object$capthist, 'list')) capthists <- object$capthist[[sessnum]]
    else capthists <- object$capthist
    if (inherits(object$mask, 'list')) masks <- object$mask[[sessnum]]
    else masks <- object$mask

    if (is.null(beta) & is.null(real)) beta <- object$fit$par

    traps  <- attr(capthists, 'traps')  ## need session-specific traps
    dettype <- detectorcode(traps)
    n      <- nrow(capthists)
    s      <- ncol(capthists)
    nmix   <- object$details$nmix
    if (is.null(nmix)) nmix <- 1

    ##############################################
    ## adapt for marking occasions only 2009 10 24
    q <- attr(capthists, 'q')
    if (!is.null(q))
        if (q<s) s <- q
    ##############################################

    if (dettype == 6) {
        k <- c(table(polyID(traps)),0)
        K <- length(k)-1
    }
    else if (dettype == 7) {
        k <- c(table(transectID(traps)),0)
        K <- length(k)-1
    }
    else {
        k <- nrow(traps)
        K <- k
    }
    m      <- length(masks$x)            ## need session-specific mask...
    cell   <- attr(masks,'area')

    if (n==0)
        stop(paste("no data in 'capthist' for session", session))
    if (is.null(beta)) {
        if (is.null(real))
            stop("requires real parameter values")
        PIA <- rep(1, n * s * K * nmix)    ## nmix added 2010 02 25

## new code 2010-11-26
        realparval0 <- matrix(rep(real, rep(n,length(real))),nr=n)   ## UNTRANSFORMED
## replacing...
##        realparval0 <- matrix(real,nr=1)   ## UNTRANSFORMED
##        ncolPIA <- 1
        Xrealparval0 <- scaled.detection (realparval0, FALSE, object$details$scaleg0, NA)
        details <- list(binomN=0, cutval=0, spherical=FALSE)
    }
    else {
        ## allow for old design object
        if (length(dim(object$design0$PIA))==4)
            dim(object$design0$PIA) <- c(dim(object$design0$PIA),1)
        PIA <- object$design0$PIA[sessnum,,1:s,,,drop=F]   ## modified 1:s 2009 10 24; extra dim 2009 12 11
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

##        Xrealparval0 <- realparval0
##        details <- object$details
##        if (details$scaleg0) {
##            realnames <- dimnames(realparval0)[2]
##            sigmaindex <- match('sigma', realnames)
##            g0index <- match('g0', realnames)
##            if (is.na(g0index) | is.na(sigmaindex)) stop('scaleg0 requires both g0 and sigma in model')
##            Xrealparval0[,g0index] <- Xrealparval0[,g0index] / Xrealparval0[,sigmaindex]^2
##        }
## simplified 2010 03 04:
        Xrealparval0 <- scaled.detection (realparval0, FALSE, object$details$scaleg0, NA)

    }

    ## new code 2010-11-26
    used <- usage(traps)
    if (any(used==0))
    PIA <- PIA * rep(rep(t(used),rep(n,s*K)),nmix)
    ncolPIA <- n

    ## space <- attr(traps(capthists), 'spacing')
    ## searchcell <- ifelse (is.null(space), 1, space^2 / 10000)  ## ha

    temp <- .C("integralprw1", PACKAGE = 'secr',
      as.integer(dettype),
      as.double(Xrealparval0),
      as.integer(n),
      as.integer(s),
      as.integer(k),
      as.integer(m),
      as.integer(nmix),
      as.double(unlist(traps)),
      as.double(unlist(masks)),
      as.integer(nrow(Xrealparval0)), # rows in lookup
      as.integer(PIA),                # index of nc*,S,K to rows in realparval0
      as.integer(ncolPIA),            # ncol - if CL, ncolPIA = n, else ncolPIA = 1 or ngrp
      as.double(cell),
      as.integer(object$detectfn),
      as.integer(object$details$binomN),
      as.double(object$details$cutval),
      a=double(n),
      resultcode=integer(1)
   )

   if (temp$resultcode == 3)
       stop ("groups not implemented in external function 'integralprw1'")
   if (temp$resultcode != 0)
       stop ("error in external function 'integralprw1'")
   temp$a
}
############################################################################################

