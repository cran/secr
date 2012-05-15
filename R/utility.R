###############################################################################
# Global variables in namespace
#
## define a local environment (in namespace?) for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html
## 2012-04-13 pointsInPolygon extended to allow mask 'polygon'
###############################################################################

.localstuff <- new.env()
.localstuff$validdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                'cue', 'unmarked','presence')
.localstuff$simpledetectors <- c('single','multi','proximity','count')
.localstuff$individualdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                     'cue')
.localstuff$pointdetectors <- c('single','multi','proximity','count',
    'signal', 'signalnoise', 'cue', 'unmarked','presence')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
## 'signal' is not a count detector 2011-02-01
.localstuff$countdetectors <- c('count','polygon','transect','unmarked')
.localstuff$detectors3D <- c('proximity','count','signal','signalnoise','polygon',
                             'transect','times','cue','unmarked','presence')
.localstuff$iter <- 0
.localstuff$detectionfunctions <-
        c('halfnormal',
      'hazard rate',
      'exponential',
      'compound halfnormal',
      'uniform',
      'w exponential',
      'annular normal',
      'cumulative lognormal',
      'cumulative gamma',
      'binary signal strength',
      'signal strength',
      'signal strength spherical',
      'signal-noise',
      'signal-noise spherical')

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}
detectionfunctionnumber <- function (detname) {
    dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn))
        stop ("unrecognised detection function ", detname)
    dfn-1
}
parnames <- function (detectfn) {
    parnames <- switch (detectfn+1,
        c('g0','sigma'),
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','z'),
        c('g0','sigma'),
        c('g0','sigma','w'),
        c('g0','sigma','w'),
        c('g0','sigma','z'),
        c('g0','sigma','z'),
        c('b0','b1'),
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS'),    ## include cutval?
        c('beta0','beta1', 'sdS','muN','sdN'),
        c('beta0','beta1', 'sdS','muN','sdN'),
        ,,,,,,
        c('g0','sigma')    ## 20
    )
}

valid.detectfn <- function (detectfn, valid = c(0:3,5:13)) {
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ("invalid detection function")
    detectfn
}

valid.detectpar <- function (detectpar, detectfn) {
    if (is.null(detectpar) | is.null(detectfn))
        stop ("requires valid 'detectpar' and 'detectfn'")
    if (!all(parnames(detectfn) %in% names(detectpar)))
        stop ("requires 'detectpar' ", paste(parnames(detectfn), collapse=','),
            " for ", detectionfunctionname(detectfn), " detectfn")
    detectpar[parnames(detectfn)]
}

detectorcode <- function (object, MLonly = TRUE) {
    ## numeric detector code from traps object
    detcode <- switch (detector(object),
        single = -1,
        multi = 0,
        proximity = 1,
        count = 2,
        polygonX = 3,
        transectX = 4,
        signal = 5,
        polygon = 6,
        transect = 7,
        times = 8,
        cue = 9,
        unmarked = 10,
        presence = 11,
        signalnoise = 12,
        -2)
    if (MLonly) {
        detcode <- ifelse (detcode==-1, 0, detcode)
        if (detcode<0)
            stop ("Unrecognised detector type")
    }
    detcode
}

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)),
        replace = T, size = n)
    else rep(tN,n)
}

ndetector <- function (traps) {
    if (detector(traps) %in% c('polygon', 'transect','polygonX', 'transectX'))
        length(levels(polyID(traps)))
    else
        nrow(traps)
}

memo <- function (text, trace) {
    ## could use message(text), but does not immediately flush console
    if (trace) { cat (text, '\n')
    flush.console() }
}

insertdim <- function (x, dimx, dims) {
  ## make vector of values
  ## using x repeated so as to fill array
  ## with dim = dims and the x values occupying dimension(s) dimx
  olddim <- 1:length(dims)
  olddim <- c(olddim[dimx], olddim[-dimx])
  temp <- array (dim=c(dims[dimx], dims[-dimx]))
  tempval <- array(dim=dims[dimx])
  if (length(x) > length(tempval))
      tempval[] <- x[1:length(tempval)]
  else
      tempval[] <- x     ## repeat as needed
  temp[] <- tempval  ## repeat as needed
  if (is.factor(x))
    factor(levels(x), levels=levels(x))[aperm(temp, order(olddim))]   ## 2010 02 25
  else
    as.vector(aperm(temp, order(olddim)))
}

pad1 <- function (x, n) {
## pad x to length n with dummy (first value)
    if (is.factor(x)) {
        xc <- as.character(x)
        xNA <- c(xc, rep(xc[1], n-length(xc)))
        out <- factor(xNA, levels=levels(x))
    }
    else out <- c(x, rep(x[1], n-length(x)))
    out
}

padarray <- function (x, dims) {
    temp <- array(dim=dims)
    dimx <- dim(x)
    if (length(dimx)<2 | length(dimx)>3)
        stop ("invalid array")
    if (length(dimx)>2) temp[1:dimx[1], 1:dimx[2], 1:dimx[3]] <- x
    else temp[1:dimx[1], 1:dimx[2]] <- x
    temp
}

## regularize a list of formulae
## added 2009 08 05
stdform <- function (flist) {
    LHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==2) '' else trms[2]
    }
    RHS <- function (form) {
        trms <- as.character (form)
        if (length(trms)==3) as.formula(paste(trms[c(1,3)])) else form
    }
    lhs <- sapply(flist, LHS)
    temp <- lapply(flist, RHS)
    if (is.null(names(flist))) names(temp) <- lhs
    else names(temp) <- ifelse(names(flist) == '', lhs, names(flist))
    temp
}

## Start of miscellaneous functions

invlogit <- function (y) 1/(1+exp(-y))   # plogis(y)
logit    <- function (x) log(x/(1-x))    # qlogis(x), except for invalid argument
sine     <- function (x) asin (x*2-1)
invsine  <- function (y) (sin(y)+1) / 2
odds     <- function (x) x / (1-x)
invodds  <- function (y) y / (1+y)

lnbinomial <- function (x,size,prob) {
  lgamma (size+1) - lgamma (size-x+1) - lgamma (x+1) +
      x * log(prob) + (size-x) * log (1-prob)
}

var.in.model <- function(v,m) v %in% unlist(lapply(m, all.vars))

get.nmix <- function (model) {
    model$D <- NULL  ## ignore density model
    nmix <- 1
    if (any(var.in.model('h2', model))) {
        nmix <- 2
        if (any(var.in.model('h3', model)))
            stop ("do not combine h2 and h3")
    }
    if (any(var.in.model('h3', model))) {
        nmix <- 3
        warning ("implementation of 3-part mixtures is not reliable")
    }
    nmix
}

add.cl <- function (df, alpha, loginterval, lowerbound = 0) {

## add lognormal or standard Wald intervals to dataframe with columns
## 'estimate' and 'SE.estimate'
## lowerbound added 2011-07-15
    z <- abs(qnorm(1-alpha/2))
    if (loginterval) {
        delta <- df$estimate - lowerbound
        df$lcl <- delta / exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
        df$ucl <- delta * exp(z * sqrt(log(1 + (df$SE.estimate /
                        delta)^2))) + lowerbound
    }
    else {
        df$lcl <- pmax(lowerbound, df$estimate - z * df$SE.estimate)
        df$ucl <- df$estimate + z * df$SE.estimate
    }
    df
}

###############################################################################

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 +
            beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
#        (0.5 - detpar$b0) / detpar$b1
        - 1 / detpar$b1   ## 2010-11-01
    }
    else stop ("unrecognised detectfn")
}

###############################################################################

pointsInPolygon <- function (xy, poly, logical = TRUE) {
    xy <- matrix(unlist(xy), ncol = 2)  ## in case dataframe
    if (inherits(poly, "SpatialPolygonsDataFrame")) {
        if (!require (sp))
            stop ("package 'sp' required for pointsInPolygon()")
        xy <- SpatialPoints(xy)
        OK <- overlay (xy, poly)
        !is.na(OK)
    }
    else if (inherits(poly, 'mask')) {   2012-04-13
        if (ms(poly))
            stop ("multi-session masks not supported")
        sp <- spacing(poly)
        minx <- min(poly$x)
        miny <- min(poly$y)
        mask <- sweep(poly, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        mask <- round(mask/sp) + 1
        xy <- sweep(xy, MARGIN = 2, FUN = '+', STATS = c(-minx, -miny))
        xy <- round(xy/sp) + 1
        xy[xy<0] <- NA
        xy[,1][xy[,1]>max(mask$x)] <- NA
        xy[,2][xy[,2]>max(mask$y)] <- NA
        maskmatrix <- matrix(0, ncol = max(mask$y), nrow = max(mask$x))
        maskmatrix[as.matrix(mask)] <- 1:nrow(mask)
        inside <- maskmatrix[as.matrix(xy)]
        inside[is.na(inside)] <- 0
        if (logical)
            inside <- inside > 0
        inside
    }
    else {
        checkone <- function (xy1) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (xy1),
                as.integer (0),
                as.integer (np-1),
                as.integer (np),
                as.double (poly),
                result = integer(1))
            as.logical(temp$result)
        }
        poly <- matrix(unlist(poly), ncol = 2)  ## in case dataframe
        np <- nrow(poly)
        apply(xy, 1, checkone)
    }
}
###############################################################################
## logical for whether object specifies userDfn

userD <- function (object) {
    if (!inherits(object, 'secr'))
        stop ("requires secr fitted model")
    !is.null(object$details$userDfn)
}

###############################################################################

## mean and SD if x numeric
## was statfn 2011-11-08
getMeanSD <- function(xy) {
    MeanSD <- function (x) {
        if (is.numeric(x))
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
        else
            c(NA,NA)
    }
    as.data.frame (apply(xy, 2, MeanSD))
}
###############################################################################
