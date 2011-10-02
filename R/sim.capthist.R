###############################################################################
## package 'secr'
## sim.capthist.R
## simulate capture histories
## 2009 10 08 sim.resight
## 2010 07 01 allow alphanumeric detection functions
## 2010 10 09 annular normal and cumulative lognormal detection functions
## 2011 03 19 allow zero detections
## 2011 03 27 multiple sessions
## 2011 09 09 p.available; session numbers
## 2011 09 30 presence, unmarked treated as special case of proximity
###############################################################################

expand <- function (x, n, q = 0, default = 1) {
    if (is.null(x)) rep(default, n)
    else {
        y <- numeric(n)
        if ((length(x)==2) && (q>0))
            y[] <- rep(x, c(q,n-q))
        else
            y[] <- x
        y
    }
}

sim.capthist <- function (
    traps,
    popn = list(D = 5, buffer = 100, Ndist = 'poisson'),
    detectfn = 0,
    detectpar = list(),
    noccasions = 5,
    nsessions = 1,
    binomN = NULL,
    p.available = 1,
    renumber = TRUE,
    seed = NULL,
    maxperpoly = 100
    )

#
# Simulate detection of known population
#

## A note on sort order of records  2009 12 01

## Simulation routines return the primary detection data in 'ski' order,
## meaning that occasion (s) changes fastest and individual (i) slowest -
## this allows efficient sequential output as new animals are detected.
## In R the order is changed to 'isk' as this is pictorially more natural.

## Secondary data (xy locations, signal strength) are generated in the order
## 'kis' (detector (k) changing fastest), because the simulation routines use -
## for (s=0; s<*ss; s++)
##   for (i=0; i<*N; i++)
##     for (k=0; k<*kk; k++) {
##     ...
##     }
## or similar loops.
## For consistency, all data are sorted to 'isk' order before returning to R.
## Secondary simulated data are re-sorted into 'ksi' order by creating the
## index 'start' within C as required for prwipolygon & prwitransect.
## The functions 'animalID', 'occasion', and 'trap' are the safest way to
## retrieve values for detections in isk order.

{
    poplist <- inherits(popn,'popn') & inherits(popn,'list')
    if (poplist | (nsessions > 1)) {

        if (poplist & (nsessions>1) & (length(popn) != nsessions))
            stop ("incompatible use of popn list and nsessions>1")

        ## supplied with spatiotemporal population
        R <- ifelse (poplist, length(popn), nsessions)
        output <- vector(R, mode='list')
        nocc <- numeric(R)
        nocc[] <- noccasions
        if (!inherits(popn,'popn')) # generate if not provided
        {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            popn <- sim.popn (popn$D, core=traps, buffer=popn$buffer,
                covariates=NULL, Ndist = popn$Ndist)
        }
        if (poplist) {
            if (any(p.available) != 1)
                warning ("incomplete availability not implemented ",
                         "for population lists")
        }
        else {
            if (!(length(p.available) %in% 1:2))
                stop ("p.available must be vector of 1 or 2 probabilities")
            availability <- 'random'
            if (length(p.available) == 1)
                ## random temporary emigration
                available <- runif(nrow(popn)) < p.available
            else {
                ## Markovian temporary emigration
                availability <- 'Markov'
                equilibrium.p <- (p.available[2] / (1-p.available[1]+p.available[2]))
                available <- runif(nrow(popn)) < equilibrium.p
            }

        }
        for (t in 1:R) {
            if (poplist)
                temppop <- popn[[t]]
            else {
                temppop <- subset(popn, available)
                ## update availability in preparation for next session
                if (availability == 'random') {
                    available <- runif(nrow(popn)) < p.available
                }
                else {
                    p.vect <- ifelse(available, p.available[1], p.available[2])
                    available <- runif(nrow(popn)) < p.vect
                }
            }
            output[[t]] <- sim.capthist(traps, temppop, detectfn, detectpar,
                nocc[t], 1, binomN, 1, renumber, seed, maxperpoly)
            session( output[[t]] ) <- t   ## added 2011-09-09

        }
        class(output) <- c('list','capthist')
        names(output) <- 1:R
        output
    }

    else {

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

        if (is.null(detector(traps)))
            stop ("'traps' lacks detector type")

        usage <- usage(traps)
        if (is.null(usage))
            usage <- matrix (1, nrow = ndetector(traps), ncol = noccasions)
        else {
            if (nrow(usage) != ndetector(traps))
                stop ("invalid usage matrix; number of rows ",
                      "must match number of detectors")
            if (ncol(usage) != noccasions) {
                noccasions <- ncol(usage)
                warning ("'noccasions' does not match usage ",
                         "attribute of 'traps'; ignored")
            }
        }


        validatepar <- function (x, xrange) {
            xname <- deparse(substitute(x))
            if (is.null(x))
                stop ("no value for ", xname)
            if (any(is.na(x)))
                stop ("NA is not a valid value for ", xname)
            if (any(x < xrange[1]))
                warning ("value for ", xname, " is less than minimum ",
                    xrange[1])
            if (any(x > xrange[2]))
                warning ("value for ", xname, " is greater than maximum ",
                    xrange[2])
        }

        ## added 2010-07-01
        if (is.character(detectfn))
            detectfn <- detectionfunctionnumber(detectfn)
        if (detector(traps) %in% c('cue','signal')) {
            if (detectfn != 10)
                warning ("forcing detection function = 10 for signal detectors")
            detectfn <- 10
        }

        ## Detection function parameters

        ##    0  halfnormal
        ##    1  hazard rate
        ##    2  exponential
        ##    3  compound halfnormal
        ##    4  uniform
        ##    5  w-exponential
        ##    6  annular normal
        ##    7  cumulative lognormal
        ##    9  binary signal strength (b0 = (beta0-c)/sdS, b1 = beta1/sdS)
        ##   10  signal strength (signal detectors only)
        ##   11  signal strength with spherical spreading (signal detectors only)

        ## extended for uniform (detectfn=4) 2010-06-13
        if (detectfn %in% c(0:4))  defaults <- list(g0 = 0.2, sigma = 25, z = 1)
        if (detectfn %in% c(5))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(6))    defaults <- list(g0 = 0.2, sigma = 25, w = 10)
        if (detectfn %in% c(7))    defaults <- list(g0 = 0.2, sigma = 25, z = 5)
        if (detectfn %in% c(9))  defaults <- list(b0 = 1, b1=-0.1, cutval = 60,
            tx = 'identity')
        if (detectfn %in% c(10,11))  defaults <- list(beta0 = 90, beta1=-0.2,
            sdS = 2, cutval = 60, sdM = 0, tx = 'identity')
        else defaults <- c(defaults, list(truncate = 1e+10))

        ## changed 2011-01-28
        defaults$binomN <- 0
        if (detector(traps) == 'proximity') defaults$binomN <- 1
        if (detector(traps) == 'signal') defaults$binomN <- 1
        if (detector(traps) == 'cue') defaults$cuerate <- 3
        if (!is.null(binomN)) detectpar$binomN <- binomN
        detectpar <- replacedefaults(defaults, detectpar)

        if (detectfn %in% c(0,1,2,3,4,5,6,7)) {
            g0    <- expand(detectpar$g0, noccasions)
            sigma <- expand(detectpar$sigma, noccasions)
            z     <- expand(ifelse(detectfn %in% c(5,6),
                detectpar$w, detectpar$z), noccasions)
            if ((detector(traps) %in% .localstuff$countdetectors) &
               (detectpar$binomN != 1))
                validatepar(g0, c(0,Inf))
            else validatepar(g0, c(0,1))
            validatepar(sigma, c(1e-10,Inf))
            validatepar(z, c(0,Inf))
        }

        # Acoustic detection function parameters
        if (detectfn %in% c(10,11)) {
            tx <- detectpar$tx
            cutval <- detectpar$cutval
            sdM <- detectpar$sdM
            beta0 <- expand(detectpar$beta0, noccasions)
            beta1 <- expand(detectpar$beta1, noccasions)
            sdS   <- expand(detectpar$sdS, noccasions)
            validatepar(beta0, c(-Inf,Inf))
            validatepar(beta1, c(-Inf,Inf))
            validatepar(sdS, c(0,Inf))
        }
        else if (detectfn %in% c(9)) {
            g0 <- expand(detectpar$b0, noccasions)
            sigma <- expand(detectpar$b1, noccasions)
            z <- 0
            cutval <- detectpar$cutval
        }
        else {
            cutval <- NULL
            truncate <- ifelse(is.null(detectpar$truncate),
                1e+10, detectpar$truncate)
            validatepar(truncate, c(1e-10, Inf)) ## must be positive
        }
        if (!inherits(popn,'popn')) # generate if not provided
        {
            popn <- replacedefaults(list(D = 5, buffer = 100,
                Ndist = 'poisson'), popn)
            popn <- sim.popn (popn$D, core=traps, buffer=popn$buffer,
                covariates=NULL, Ndist = popn$Ndist)
        }

        ################################
        ################################
        if (detector(traps) == 'cue') {
            rowi <- 1:nrow(popn)
            ncue <- rpois(nrow(popn), detectpar$cuerate)
            OK <- ncue>0
            group <- rowi[OK]   ## drop zeros
            ncue <- ncue[OK]
            ind <- rep(group, ncue)
            N <- length(ind)
            animals <- unlist(popn[ind,])
        }
        else {
            N <- nrow(popn)
            animals <- unlist(popn)
        }
        k <- nrow(traps)
        ################################

        if (detector(traps) %in% c('single','multi')) {
            simfunctionname <- paste('trapping', detector(traps), sep='')
            temp <- .C(simfunctionname, PACKAGE = 'secr',
                as.double(g0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                as.double(truncate^2),
                n = integer(1),
                caught = integer(N),
                value=integer(N*noccasions),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to '", simfunctionname, "' failed")
            w <- matrix(ncol = temp$n, nrow = noccasions, dimnames =
                 list(1:noccasions, NULL))
            if (temp$n > 0) w[,] <- temp$value[1:(temp$n*noccasions)]
            w <- t(w)
        }

        else
        if (detector(traps) %in% c('polygonX','transectX')) {
            simfunctionname <- paste('trapping', detector(traps), sep='')
            if (simfunctionname == 'trappingpolygonX') {
                nk <- length(levels(polyID(traps)))
                k <- table(polyID(traps))
            }
            else {
                nk <- length(levels(transectID(traps)))
                k <- table(transectID(traps))
            }
            temp <- .C(simfunctionname, PACKAGE = 'secr',
                as.double(g0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(nk),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                as.double(truncate^2),
                n = integer(1),
                caught = integer(N),
                detectedXY = double (N*noccasions),
                value = integer(N*noccasions),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to '", simfunctionname, "' failed")

            w <- matrix(ncol = temp$n, nrow = noccasions, dimnames =
                list(1:noccasions, NULL))
            if (temp$n > 0) {
                w[,] <- temp$value[1:(temp$n*noccasions)]
            }
            w <- t(w)

            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w) > 0)
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                attr(w, 'detectedXY') <- detectedXY
            }
            else
                attr(w, 'detectedXY') <- NULL

        }

        else
        ## includes presence 2011-09-26
        if (detector(traps) %in% c('proximity', 'count', 'presence','unmarked')) {
            binomN <- switch(detector(traps), proximity=1,
                             count=detectpar$binomN, presence=1, unmarked = 1)
            temp <- .C("trappingcount", PACKAGE = 'secr',
                as.double(g0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                as.double(truncate^2),
                as.integer(binomN),
                n = integer(1),
                caught = integer(N),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingcount' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                list(1:noccasions,NULL, NULL))
            if (temp$n>0) {
                w[,,] <- temp$value[1:(temp$n*noccasions*k)]
            }
            w <- aperm(w, c(3,1,2))
        }
        else

        if (detector(traps) %in% c('signal','cue')) {
            if (detectpar$binomN != 1)
                stop ("binomN != 1 not yet implemented for signal detectors")
            temp <- .C("trappingsignal", PACKAGE = 'secr',
                as.double(beta0),
                as.double(beta1),
                as.double(sdS),
                as.double(cutval),
                as.double(sdM),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                n = integer(1),
                caught = integer(N),
                signal = double(N*noccasions*k),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingsignal' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                 list(1:noccasions,NULL,NULL))
            if (temp$n>0)  {
                w[,,] <- temp$value[1:(temp$n * noccasions * k)]
            }
            w <- aperm(w, c(3,1,2))
            if (temp$n>0)  {
                attr(w, 'signal') <- temp$signal[1:sum(w)]
            }
        }
        else
        if (detector(traps) == 'times') {
            temp <- .C("trappingtimes", PACKAGE = 'secr',
                as.double(g0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                as.double(truncate^2),
                n = integer(1),
                caught = integer(N),
                times = double(N*noccasions*k),
                value = integer(N*noccasions*k),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0)
                stop ("call to 'trappingtimes' failed")
            w <- array(dim=c(noccasions, k, temp$n), dimnames =
                list(1:noccasions,NULL,NULL))
            if (temp$n>0)  {
                w[,,] <- temp$value[1:(temp$n * noccasions * k)]
            }
            w <- aperm(w, c(3,1,2))
            if (temp$n>0)  {
                attr(w, 'times') <- temp$times[1:sum(w)]
            }
        }
        else if (detector(traps) %in% c('polygon','transect')) {
            simfunctionname <- paste('trapping', detector(traps), sep='')
            if (simfunctionname == 'trappingpolygon') {
                nk <- length(levels(polyID(traps)))
                k <- table(polyID(traps))
            }
            else {
                nk <- length(levels(transectID(traps)))
                k <- table(transectID(traps))
            }

            temp <- .C(simfunctionname, PACKAGE = 'secr',
                as.double(g0),
                as.double(sigma),
                as.double(z),
                as.integer(noccasions),
                as.integer(nk),
                as.integer(k),
                as.integer(N),
                as.double(animals),
                as.double(unlist(traps)),
                as.integer(usage),
                as.integer(detectfn),
                as.double(truncate^2),
                as.integer(detectpar$binomN),
                as.integer(maxperpoly),
                n = integer(1),
                caught = integer(N),
                detectedXY = double (N*noccasions*nk*maxperpoly*2),
                value = integer(N*noccasions*nk),
                resultcode = integer(1)
            )
            if (temp$resultcode != 0) {
                if (temp$resultcode == 2)
                    stop ("more than ", maxperpoly, "  detections per animal",
                          " per polygon per occasion")
                else
                    stop ("call to ",simfunctionname, " failed")
            }
            w <- array(dim=c(noccasions, nk, temp$n),
                dimnames = list(1:noccasions, levels(polyID(traps)), NULL))

            if (temp$n > 0) {
                w[,,] <- temp$value[1:prod(dim(w))]
            }
            w <- aperm(w, c(3,1,2))

            if (temp$n > 0) {
                ## put XY coordinates in attribute
                nd <- sum(abs(w))
                detectedXY <- data.frame(matrix(ncol = 2,
                    temp$detectedXY[1:(2*nd)]))
                names(detectedXY) <- c('x','y')
                attr(w, 'detectedXY') <- detectedXY
            }
            else
                attr(w, 'detectedXY') <- NULL
        }
        else stop ('Unrecognised detector type')

        if (!is.null(covariates(popn))) {
            covariates(w) <- covariates(popn)[as.logical(temp$caught),, drop=F]
        }
        if (detector(traps)=='cue') {
            group <- ind[as.logical(temp$caught)]
            covariates(w) <- data.frame(animal = factor(group))   ## overrides
        }

        class(w)             <- 'capthist'    ## NOT data.frame
        traps(w)             <- traps
        attr(w, 'cutval')    <- cutval
        attr(w, 'seed')      <- RNGstate      ## save random seed
        attr(w, 'detectpar') <- detectpar
        session(w)           <- '1'           ## dummy session values for now

        if (renumber && (temp$n>0)) rownames(w) <- 1:temp$n
        ##   else rownames(w) <- (1:N)[as.logical(temp$caught)]
        ## 2011-04-02 BUG FIX
        else {
#            rown <- (1:N)[temp$caught>0]
# test 2011-09-11
rown <- rownames(popn)[temp$caught > 0]
            caught <- temp$caught[temp$caught>0]
            rownames(w) <- rown[order(caught)]
        }

        w
    }
}
############################################################################################

sim.resight <- function (traps, ..., q = 1, pID = 1, unmarked = TRUE,
    nonID = TRUE) {

    defaultpar <- list(noccasion = 5)
    dots <- list(...)
    dots <- replace (defaultpar, names(dots), dots)
    dots$detectpar$g0 <- expand (dots$detectpar$g0, dots$noccasion, q)
    dots$detectpar$sigma <- expand (dots$detectpar$sigma, dots$noccasion, q)
    dots$detectpar$z <- expand (dots$detectpar$z, dots$noccasion, q)

    capthist <- do.call('sim.capthist', c(list(traps = traps), dots))

    ## transform simulated capthist object into resighting data
    S <- ncol(capthist)
    K <- nrow(traps(capthist))

    if (S <= q)
        stop ("no sighting intervals")
    if (!(detector(traps(capthist)) %in% c('proximity')))
        stop ("only for proximity detectors")

    ## sighting only unmarked animals
    marked <- apply(capthist[,1:q,,drop = FALSE], 1, sum) > 0
    ## always get warning about no detections on first occasion
    suppressWarnings(R <- subset(capthist, subset=!marked, dropnull=F))
    tempM <- subset(capthist, subset=marked, dropnull=F)
    nM <- nrow(tempM)
    ID <- abind(tempM[,1:q, , drop=F],
        array(runif(nrow(tempM)*(S-q)*K) < pID, dim=c(nM,S-q,K)), along=2)
    ## marking and sighting, marked animals
    M <- ifelse(ID, tempM,0)    ## ID
    U <- ifelse(!ID, tempM,0)   ## notID
    dimnames(M)[[2]] <- 1:S
    countfn <- function(x) {
        x <- x * col(x)   ## x 0/1
        tabulate(x[x>0], nbins = K)
    }

    if (unmarked) {
        Tu <- apply(R, 2, countfn)  ## not marked
        row.names(Tu) <- row.names(traps(capthist))
    }
    else Tu <- NULL

    if (nonID) {
        Tm <- apply(U, 2, countfn)  ## notID, marked
        dimnames(Tm) <- dimnames(Tu)
    }
    else Tm <- NULL

    class(M) <- 'capthist'
    session(M) <- session(capthist)
    traps(M) <- traps(capthist)

    attr(M, 'Tu') <- Tu
    attr(M, 'Tm') <- Tm
    attr(M, 'q') <- q

    M
}
############################################################################################

