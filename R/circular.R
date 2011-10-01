################################################################################## package 'secr'
## circular.R
## "circular error probable"
## last changed 2011 06 12
################################################################################
# Functions from plot.secr.R
# HN <- function (r, pars, cutval)
# HZ <- function (r, pars, cutval)
# EX <- function (r, pars, cutval)
# UN <- function (r, pars, cutval)
# CHN <- function (r, pars, cutval)
# WEX <- function (r, pars, cutval)
# ANN <- function (r, pars, cutval)
# CLN <- function (r, pars, cutval)
# CG <- function (r, pars, cutval)
# CN <- function (r, pars, cutval)
# BSS <- function (r, pars, cutval)
# SS <- function (r, pars, cutval)
# SSS <- function (r, pars, cutval)

## source ('d:\\density secr 2.1\\secr\\r\\CEP.R')
# source ('d:\\density secr 2.1\\secr\\r\\plot.secr.R')  # dfn's
# source ('d:\\density secr 2.1\\secr\\r\\utility.R')    # parnames
# source ('d:\\density secr 2.1\\secr\\r\\pdot.R')    # spatial scale

circular.r <- function (p = 0.95, detectfn = 0, sigma = 1, detectpar = NULL) {

    ## convert character detectfn to numeric code
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (detectfn %in% c(0,2,3)) {
        ## if input is a named list
        if (!is.null(detectpar)) {
            sigma <- detectpar$sigma
        }
        detectpar <- list(sigma = sigma)
    }
    else
        if (is.null(detectpar))
            stop ("require detectpar in list format except for ",
                  "halfnormal, exponential and uniform")
    detectpar$g0 <- 1  ## always
    detectpar <- detectpar[parnames(detectfn)]  ## correct order
    pars <- unlist(detectpar)
    cutval <- ifelse (detectfn %in% c(9,10,11), detectpar$cutval, NA)
    scale <- spatialscale (detectpar, detectfn) ## see pdot.R; assumes cutval in detectpar

    ## use formula for halfnormal
    if (detectfn == 0) {
        (-2*log(1-p))^0.5 * sigma
    }
    ## uniform is dead easy
    else if (detectfn == 3) {
        p^0.5 * sigma
    }
    ## otherwise integrate
    else {
        dfn <- switch (detectfn+1, HN, HZ, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS)
        rdfn <- function (r, pars, cutval)
            r * dfn(r, pars, cutval)
        I1 <- integrate (rdfn, 0, Inf, pars, cutval)$value

        fnr <- function (r, this.p) {
            I2 <- integrate (rdfn, 0, r, pars, cutval)$value
            I2 / I1 - this.p
        }
        getroot <- function (p) uniroot(fnr, c(0,100*scale), this.p = p)$root
        sapply(p, getroot)
    }
}

circular.p <- function (r = 1, detectfn = 0, sigma = 1, detectpar = NULL) {

    ## convert character detectfn to numeric code
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (detectfn %in% c(0,2,3)) {
        ## if input is a named list
        if (!is.null(detectpar)) {
            sigma <- detectpar$sigma
        }
        detectpar <- list(sigma = sigma)
    }
    else
        if (is.null(detectpar))
            stop ("require detectpar in list format except for ",
                  "halfnormal, exponential and uniform")
    detectpar$g0 <- 1  ## always
    detectpar <- detectpar[parnames(detectfn)]  ## correct order
    pars <- unlist(detectpar)
    cutval <- ifelse (detectfn %in% c(9,10,11), detectpar$cutval, NA)

    scale <- spatialscale (detectpar, detectfn) ## see pdot.R; assumes cutval in detectpar


    ## use formula for halfnormal
    if (detectfn == 0) {
        1 - exp(-(r/sigma)^2 / 2)
    }
    ## uniform is dead easy
    else if (detectfn == 3) {
        (r/sigma)^2
    }
    ## otherwise integrate
    else {
        dfn <- switch (detectfn+1, HN, HZ, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS)
        rdfn <- function (r, pars, cutval)
            r * dfn(r, pars, cutval)
        I1 <- integrate (rdfn, 0, Inf, pars, cutval)$value

        fnr <- function (r) {
            I2 <- integrate (rdfn, 0, r, pars, cutval)$value
            I2 / I1
        }
        sapply(r, fnr)
    }
}

# plot(seq(0,5,0.1), circular.p(seq(0,5,0.1), detectfn=0), type='l', xlab='radius', ylab='p')
# lines(seq(0,5,0.1),circular.p(seq(0,5,0.1), detectfn=1, detectpar=list(sigma=1, z=4)), col='blue')
# lines(seq(0,5,0.1),circular.p(seq(0,5,0.1), detectfn=2), col='red')



