############################################################################################
# Global variables in namespace
#
## define a local environment (in namespace?) for temporary variables e.g. iter
## see e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html

.localstuff <- new.env()
.localstuff$validdetectors <- c('single','multi','proximity','count', 'polygonX', 'transectX',
                                'signal', 'polygon', 'transect', 'times','cue')
.localstuff$pointdetectors <- c('single','multi','proximity','count','signal','cue')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
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
      'signal strength spherical')

## 'signal' is not a count detector 2011-02-01
.localstuff$countdetectors <- c('count','polygon','transect')

detectionfunctionname <- function (fn) {
    .localstuff$detectionfunctions[fn+1]
}
detectionfunctionnumber <- function (detname) {
    dfn <- match (tolower(detname), .localstuff$detectionfunctions)
    if (is.na(dfn)) stop(paste("unrecognised detection function", detname))
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
        c('beta0','beta1', 'sdS'),     ## include cutval?
        ,,,,,,,,
        c('g0','sigma')    ## 20
    )
}

valid.detectfn <- function (detectfn, valid = c(0:3,5:11)) {
    if (is.null(detectfn))
        stop ("requires 'detectfn'")
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if (!(detectfn %in% valid))
        stop ('invalid detection function')
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
        polygon=6,
        transect=7,
        times=8,
        cue=9,
        -2)
    if (MLonly) {
        detcode <- ifelse (detcode==-1, 0, detcode)
        if (detcode<0) stop('Unrecognised detector type')
    }
    detcode
}

replacedefaults <- function (default, user) replace(default, names(user), user)

discreteN <- function (n, N) {
    tN <- trunc(N)
    if (N != tN) tN + sample (x = c(1,0), prob = c(N-tN, 1-(N-tN)), replace = T, size = n)
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
