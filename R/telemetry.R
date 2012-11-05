###############################################################################
## package 'secr'
## telemetry.R
## towards fitting detection function to telemetry data
## 2012-10-19,20,21,22
###############################################################################

distances <- function (X, Y) {
    ## X and Y are 2-column matrices of coordinates
    onerow <- function (xy) {
        d <- function(xy2) {
            sqrt(sum((xy2 - xy)^2))
        }
        apply(Y, 1, d)
    }
    t(apply(X, 1, onerow))
}

telemetryloglik <- function(CH, detectfn, sigma, z) {
    ## assume one occasion
    dfn <- switch (detectfn+1, HN, HZ, EX, CHN, UN, WEX, ANN, CLN, CG, BSS, SS, SSS)
    rdfn <- function (r, pars) r * dfn(r, pars, 0)
    # for HN cf Johnson & Wichern 4.4, noting det(vcv) for 2-D involves sigma^2
    g0 <- switch (detectfn+1, 1 / (sigma^2 * 2 * pi),
                  1/integrate (rdfn, 0, sigma * 20, pars=c(1, sigma, z))$value)
    parm <- c(g0, sigma, z)
    xylist <- attr(CH,'xylist')
    if (is.null(xylist))
        xylist <- split(xy(CH), animalID(CH))
    oneID <- function (xy) {
        centrexy <- matrix(apply(xy, 2, mean), nrow = 1)
        d <- distances(centrexy, xy)
        dfn(d, parm, 0)
    }
    L1 <- unlist(lapply(xylist, oneID))
    list(value = sum(log(L1)), resultcode = 0)
}

addTelemetry <- function (detectionCH, telemetryCH) {
    ## combine capture histories from telemetry and hair snags etc.
    if (ms(detectionCH) | ms(telemetryCH))
        stop("addTelemetry is not ready for multi-session inputs")
    if (!detector(traps(detectionCH)) %in% c('proximity','count'))
        stop ("detectionCH should be of detector type 'proximity' or 'count'")
    if (!detector(traps(telemetryCH))=='telemetry')
        stop ("telemetryCH should be of detector type 'telemetry'")
    telemID <- animalID(telemetryCH)
    telemxy <- xy(telemetryCH)
    xylist <- split(telemxy, telemID)
    proxID <- row.names(detectionCH)
    OK <- names(xylist) %in% proxID
    # construct empty histories for telemetry animals not caught
    zerohist <- array(0, dim = c(sum(!OK), dim(detectionCH)[2:3]))
    dimnames(zerohist)[[1]] <- names(xylist)[!OK]
    class(zerohist) <- 'capthist'
    traps(zerohist) <- traps(detectionCH)
    covariates(zerohist) <- covariates(telemetryCH)[!OK,,drop = FALSE]
    # combine with true histories
    CH <- rbind.capthist(detectionCH, zerohist, renumber = FALSE, verify=FALSE)
    attr(CH, 'xylist') <- xylist # for both non-empty and empty CH
    CH
}

outsidemask <- function(capthist, mask, threshold = spacing(mask) / sqrt(2)) {
    xylist <- attr(capthist, 'xylist')
    dfun <- function(xy) {
        centre <- matrix(apply(xy, 2, mean), ncol = 2)
        distancetotrap(centre, mask)
    }
    sapply(xylist, dfun) > threshold
}

telemloglik <- function(capthist, traps, mask, detectfn = 0, detectpar) {
    ## experimental
    n <- dim(capthist)[1]
    J <- dim(capthist)[2]
    K <- dim(capthist)[3]
    ## detection functions are in plot.secr.r
    g <- switch (detectfn+1, HN, HZ, EX, CHN, UN, WEX, ANN, CLN, CG,
                 BSS, SS, SSS, SN, SNS)
    xylist <- attr(capthist, 'xylist')
    centrexy <- t(sapply(xylist, apply, 2, mean))

    ## animal x mask
    dfn <- function (xy) {
        centre <- matrix(apply(xy, 2, mean), ncol = 2)
        vcv <- var(xy)/nrow(xy)
        detS <- det(vcv)^0.5  ## sqrt(generalised variance)
        tempmask <- sweep (mask, MARGIN = 2, STATS = centre, FUN = '-')
        close <- apply(tempmask,1, function(x) sum(x^2)) < (detS*30)
        if (sum(close)<1) {
            close <- nearesttrap(matrix(0, ncol = 2, nrow = 1), tempmask)
        }
        tempmask <- tempmask[close,, drop = FALSE]
        # require(mvtnorm)
        # dbvn <- apply(tempmask, 1, dmvnorm, mean=c(0,0), sigma=vcv, log=FALSE)
        # mymvn is faster and produces identical result
        invS <- solve(vcv)
        mymvn <- function(XY) exp(-(XY %*% invS %*% XY)/2) / 2/pi/det(vcv)
        dbvn <- apply(tempmask,1,mymvn)
        dbvn <- dbvn / sum(dbvn)
        list(dbvn=dbvn, mask = mask[close,, drop=FALSE])
    }
    pmask <-  lapply(xylist, dfn)
    loglik <- 0
    for (id in names(xylist)) {
        m <- nrow(pmask[[id]]$mask)
        dtrap <- distances (traps, pmask[[id]]$mask)
        if (m==1) dtrap <- t(dtrap)  ## kludge to deal with vector
        gkm <- g(dtrap, unlist(detectpar), 0)   # 3rd arg is dummy for cutval
        prw <- matrix(1, nrow = K, ncol = m)
         for (j in 1:J) {
             wij <- capthist[id,j,]  ## detection sites on occasion j, assume 0/1
             prw <- prw *
                 (sweep(gkm, STATS=wij, MARGIN=1, FUN = '*') +     # detected
                 sweep(1-gkm, STATS=1-wij, MARGIN=1, FUN = '*'))   # not detected
         }
        prwm <- apply(prw,2,prod)            ## product over detectors
        prwm <- prwm * pmask[[id]]$dbvn       ## m-length vectors
        loglik <- loglik + log(sum(prwm))    ## sum prob over mask_i points
    }
    loglik
}

# telemloglik (temp, traps(temp), make.mask(traps(temp), buffer=200), detectpar=c(0.2,30))
