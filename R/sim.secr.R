############################################################################################
## package 'secr'
## sim.secr.R
## Copyright (C) 2009, 2010 Murray Efford
## last changed 2009 09 07; 2009 09 19 verify = FALSE; 2009 09 26 remove 'attach'
## 2009 11 03 nullval
## 2010 02 05 removed 'spherical' argument (now detectfn 11)
## 2010 03 10 debugged simsecr in secr.c
## 2010 03 10 debugged dummyCH
## 2010 06 30 memory allocation error in sim.detect
############################################################################################

simulate.secr <- function (object, nsim = 1, seed = NULL, chat = 1, ...)
## if CL, condition on n? what about distribution of covariates over n?
## locate according to IHP with lambda(X) controlled by f(X|covar), assuming homog Poisson
## i.e. use f(X|covar)/max(f(X|covar)) to reject until meet quota n?
## or f(X, covar | detected)?
## TOO HARD - cf MARK

{

##  check input
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    if (object$CL)
        stop ("not implemented for conditional likelihood")
    if (!all(sapply(object$fixed, is.null)))
        stop ("not implemented for fixed parameters")
    if (is.null(object$D))
        stop("old 'secr' object does not have 'D' component")

## setup
    # dim(object$D)[1] is number of mask points
    ngrp <- dim(object$D)[2]
    nsession <- dim(object$D)[3]
    if (!is.null(object$groups)) {
        ## individual covariates for foundation of g
        di <- disinteraction (object$capthist, object$groups)
    }
    ## we will grow this list - inefficient for very large nsim
    sesscapt <- list()

    ##################
    ## set random seed
    ## copied from simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ##################

    ## loop over replicates
    for (i in 1:nsim) {
        sesspopn <- list()
        for (sessnum in 1:nsession) {
            if (nsession==1) mask <- object$mask
            else mask <- object$mask[[sessnum]]
            popn <- list()
            for (g in 1:ngrp) {
                if (object$model$D == ~1) {
                    density <- object$D[1,g,sessnum]   ## homogeneous
                    mod2D <- 'poisson'
                }
                else {
                    density <- object$D[,g,sessnum]    ## inhomogeneous
                    mod2D <- 'IHP'
                }
                if (chat > 1)
                    density <- density / chat
                popn[[g]] <- sim.popn (D = density, core = mask, model2D = mod2D)

                ## following line replaces any previous individual covariates
                ## ---groups only---
                if (!is.null(object$groups)) {
                    covariates(popn[[g]]) <- di[rep(g, nrow(popn[[g]])),]
                }
            }
            sesspopn[[sessnum]] <- rbind.popn(popn)   ## combine groups in one popn object
        }
        sesscapt[[i]] <- sim.detect(object, object$fit$par, sesspopn)

        ## experimental
        if (chat>1)
            sesscapt[[i]] <- replicate (sesscapt[[i]], chat)
    }
    attr(sesscapt,'seed') <- RNGstate   ## save random seed
    class(sesscapt) <- c('list', 'secrdata')
    sesscapt
}
############################################################################################

sim.secr <- function (object, nsim = 1,
    extractfn = function(x) c(deviance=deviance(x), df=df.residual(x)),
    seed = NULL, data = NULL, tracelevel = 1, hessian = 'none',
    start = object$fit$par)  {

## parametric bootstrap simulations based on a fitted secr object
## extractfn is a required function to extract values from an secr fit
## it should return a vector of named values that does not vary in length
## 'hessian' overrides value in object$details
## set hessian='auto' if extractfn requires variance-covariance matrix

    cl   <- match.call(expand.dots = TRUE)
    cl <- paste(names(cl)[-1],cl[-1], sep=' = ', collapse=', ' )
    cl <- paste('sim.secr(', cl, ')')

    if (is.null(extractfn)) extractfn <- trim
    test <- extractfn(object)

   if (is.numeric(test)) {
        n.extract <- length(test)
        if (n.extract<=0)
            stop ("invalid 'extractfn'")
    }
    detectnames <- names(object$design0[[1]])   ## names of real detection parameters
    details <- replace(object$details, 'hessian', hessian)
    tracelevel <- as.integer(tracelevel)
    details$trace <- tracelevel > 1
    min.detections <- 1
    i <- 0

    if (detector(traps(object$capthist)) == 'single') {
        memo ('multi-catch likelihood used for single-catch traps', tracelevel>0)
    }

    if (is.null(data)) {
        memo ('sim.secr simulating detections...', tracelevel>0)
        data <- simulate(object, nsim = nsim, seed = seed)
    }
    else {
        if (any(class(data) != c('list','secrdata')))
            stop("invalid data")
    }

    fitmodel <- function (sc) {
        i <<- i+1
        memo (paste('sim.secr fitting replicate',i,'...'), tracelevel>0)
        nc <-  sum(counts(sc)$'M(t+1)'[,'Total'])
        if (nc >= min.detections) {
            tempfit <- suppressWarnings( secr.fit(sc, model = object$model, mask = object$mask,
                CL = object$CL, detectfn = object$detectfn, start = start, link = object$link,
                fixed = object$fixed, timecov = object$timecov, sessioncov = object$sessioncov,
                groups = object$groups, dframe = object$dframe, details = details,
                method = object$fit$method, verify = FALSE) )
            extractfn(tempfit)
        }
        else if (is.list(test)) list() else rep(NA, n.extract)
    }

    if (is.numeric(test)) {
        output <- data.frame(t(sapply (data, fitmodel)))
    }
    else {
        output <- lapply (data, fitmodel)
        class(output) <- c('list','secrlist')
    }

    attr(output,'seed') <- attr(data,'seed')
    attr(output,'call') <- cl
    attr(output,'extractfn') <- extractfn
    output
}
############################################################################################

print.secrdata <- function(x,...) {
## suggestion of Rolf Turner 19 Jan 2009 for printing without attributes
    attributes(x) <- NULL
    print(x)
}
############################################################################################

print.secrlist <- function(x,...) {
## suggestion of Rolf Turner 19 Jan 2009 for printing without attributes
    attributes(x) <- NULL
    print(x)
}
############################################################################################

sim.detect <- function (object, beta, popnlist, renumber = TRUE)
## popnlist is always a list of popn objects
{

#    stop('sim.detect not yet working in version 1.4.1, sorry')

    Markov <- 'B' %in% object$vars
    dummycapthist<- function (capthist, pop, fillvalue=1) {
        if (inherits(capthist, 'list')) {
            output <- list()
            for (i in 1:nsession)
                output[[i]] <- dummycapthist (capthist[[i]],
                    pop=pop[i], fillvalue = fillvalue)
            class(output) <- c('list','capthist')
            session(output) <- session(capthist)   ## 2010 03 10
            output
        }
        else {
            newdim <- dim(capthist)
            newdim[1] <- nrow(pop[[1]])
            output <- array(fillvalue, dim = newdim)
            class(output) <- 'capthist'
            traps(output) <- traps(capthist)
            session(output) <- session(capthist)
            covariates(output) <- covariates(pop[[1]])
            output
        }
    }
    ## setup
    MS <- inherits(object$capthist,'list')
    N <- sapply(popnlist, nrow)
    sessionlevels <- session(object$capthist)   ## was names(capthist) 2009 08 15
    nsession <- length(sessionlevels)

    ## design matrices etc.
    dummyCH <- dummycapthist(object$capthist, popnlist, fillvalue = 0)
    design0 <- secr.design.MS (dummyCH, object$model, object$timecov, object$sessioncov, object$groups, object$dframe)
    realparval0 <- makerealparameters (design0, beta, object$parindx, object$link, object$fixed)  # naive

    if ('b' %in% object$vars) {
        dummyCH <- dummycapthist(object$capthist, popnlist, fillvalue = 1)
        design1 <- secr.design.MS (dummyCH, object$model, object$timecov, object$sessioncov, object$groups, object$dframe)
        realparval1 <- makerealparameters (design1, beta, object$parindx, object$link, object$fixed)  # caught before
    }
    else {   ## faster
        design1 <- design0
        realparval1 <- realparval0
    }

    output <- list()
    for (sessnum in 1:nsession) {
        ## in multi-session case must get session-specific data from lists
        if (MS) {
            s <- ncol(object$capthist[[sessnum]])
            session.traps    <- traps(object$capthist)[[sessnum]]
        }
        else {
            s <- ncol(object$capthist)
            session.traps    <- traps(object$capthist)
        }

        dettype <- detectorcode(session.traps, MLonly = FALSE)
        if (dettype < -1)
            stop ("detector type ",
                  detector(session.traps),
                  " not implemented")

        if (dettype %in% c(2)) {                      # count detectors
            ## 2010-12-01
            ## suspend use of traps attribute until resolve vector/scalar
            ## binomN <- attr(session.traps, 'binomN')   # <0 negBin, 0 Poisson, >0 binomial
            ## if (is.null(binomN))...
            binomN <- object$details$binomN
            if (is.null(binomN))
                binomN <- 0                           # default Poisson
        }
        else binomN <- 0   # not used, just place holder

        if (dettype %in% c(6,7)) {
            k <- c(table(polyID(session.traps)),0)
            K <- length(k)-1
        }
        else {
            k <- nrow(session.traps)
            K <- k
        }
        trps  <- unlist(session.traps, use.names=F)
        sessg <- min (sessnum, design1$R)
        session.animals <- unlist(popnlist[[sessnum]])

        #------------------------------------------
        # allow for scaling of detection parameters

        Xrealparval1  <- realparval1
        Xrealparval0 <- realparval0
        ## D assumed constant over mask, groups
        sigmaindex <- 2
        g0index <- 1
        if (object$details$scalesigma) {   ## assuming previous check that scalesigma OK...
            Xrealparval1[,sigmaindex] <- Xrealparval1[,sigmaindex] / D[1,1,sessnum]^0.5
            Xrealparval0[,sigmaindex] <- Xrealparval0[,sigmaindex] / D[1,1,sessnum]^0.5
        }
        if (object$details$scaleg0)    {   ## assuming previous check that scaleg0 OK...
            Xrealparval1[,g0index] <- Xrealparval1[,g0index] / Xrealparval1[,sigmaindex]^2
            Xrealparval0[,g0index] <- Xrealparval0[,g0index] / Xrealparval0[,sigmaindex]^2
        }
        #------------------------------------------
        ## simulate this session...
        used <- unlist(usage(session.traps))
        if (is.null(used)) used <- rep(1,s*K)
        NR <- N[sessnum]

## 2010-06-30 following call tries to allocate too much memory
## 2010-12-02 fixed by reducing safety margin for detectedXY, signal
##            from 100 to 10

        temp <- .C('simsecr', PACKAGE = 'secr',
            as.integer(dettype),
            as.double(Xrealparval0),
            as.double(Xrealparval1),
            as.integer(nrow(Xrealparval0)),                # number of rows in lookup table, naive
            as.integer(nrow(Xrealparval1)),                # ditto, caught before
            as.integer(design0$PIA[sessg,1:NR,1:s,1:K,]),  # index of N,S,K to rows in Xrealparval0
            as.integer(design1$PIA[sessg,1:NR,1:s,1:K,]),  # index of N,S,K to rows in Xrealparval1
            as.integer(NR),
            as.integer(s),
            as.integer(k),
            as.integer(object$details$nmix),
            as.double(session.animals),
            as.double(trps),
            as.integer(used),
            as.integer(Markov),
            as.integer(binomN),                # used only for count detector checked 2010-12-01
            as.double(object$details$cutval),  # detection threshold on transformed scale
            as.integer(object$detectfn),
            n = integer(1),
            caught = integer(NR),
            detectedXY = double (NR*s*K*20),  # safety margin 10 detections per animal per detector per occasion
            signal = double (NR*s*K*10),      # safety margin 10 detections per animal per detector per occasion
            value = integer(NR*s*K),
            resultcode = integer(1)
        )

        if (temp$resultcode != 0) {
          if ((temp$resultcode == 2) && (dettype %in% (6:7)))
              stop (">100 detections per animal per polygon per occasion")
          else
              stop ("simulated detection failed, code ",
                    as.character(temp$resultcode))
        }
        if (dettype %in% c(-1,0)) {
            w <- array(dim=c(s,temp$n), dimnames = list(1:s, NULL))
            if (temp$n>0)  w[,] <- temp$value[1:(temp$n*s)]
            w <- t(w)
        }
        else {
            w <- array(dim=c(s, K, temp$n), dimnames = list(1:s,NULL,NULL))
            if (temp$n>0)  w[,,] <- temp$value[1:(temp$n*s*K)]
            w <- aperm(w, c(3,1,2))
        }
        class(w)   <- 'capthist'    # NOT data.frame
        traps(w)   <- session.traps
        session(w) <- sessionlevels[sessnum]

        if (!is.null(covariates(popnlist))) {
            covariates(w) <- subset(covariates(popnlist[[sessnum]]), subset =
                as.logical(temp$caught))
        }
        if ((dettype %in% c(5)) && (temp$n>0)) {
            nd <- sum(abs(w))
            attr(w, 'signal') <- temp$signal[1:nd]
            attr(w, 'cutval') <- object$details$cutval
        }
        else {
            attr(w, 'signal') <- NULL
            attr(w, 'cutval') <- NULL
        }
        if ((dettype %in% c(6:7)) && (temp$n>0)) {
            nd <- sum(abs(w))
            detectedXY <- data.frame(matrix(nc=2, temp$detectedXY[1:(2*nd)]))
            names(detectedXY) <- c('x','y')
            attr(w, 'detectedXY') <- detectedXY
        }
        else
            attr(w, 'detectedXY') <- NULL

        if (renumber & (temp$n>0)) rownames(w) <- 1:temp$n
        else rownames(w)          <- (1:N[sessnum])[as.logical(temp$caught)]
        output[[sessnum]] <- w
    }

    if (nsession==1) output <- output[[1]]
    else {
        names(output) <- sessionlevels
        class(output) <- c('list','capthist')
    }
    output
}
############################################################################################

