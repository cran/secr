###############################################################################
# Global variables in namespace
#
## define a local environment (in namespace?) for temporary variables e.g. iter
## e.g. Roger Peng https://stat.ethz.ch/pipermail/r-devel/2009-March/052883.html
## 2012-04-13 pointsInPolygon extended to allow mask 'polygon'
## 2012-07-25 nclusters() added
## 2012-09-04 leadingzero moved here
## 2012-10-21 HN, HZ etc moved here
## 2012-10-28 model.string moved here
## 2012-11-03 gradient moved here
## 2012-11-03 several other functions moved here from functions.r
###############################################################################

.localstuff <- new.env()
.localstuff$validdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                'cue', 'unmarked','presence','telemetry')
.localstuff$simpledetectors <- c('single','multi','proximity','count')
.localstuff$individualdetectors <- c('single','multi','proximity','count',
    'polygonX', 'transectX', 'signal', 'signalnoise', 'polygon', 'transect', 'times',
                                     'cue','telemetry')
.localstuff$pointdetectors <- c('single','multi','proximity','count',
    'signal', 'signalnoise', 'cue', 'unmarked','presence')
.localstuff$polydetectors <- c('polygon','transect','polygonX','transectX','telemetry')
.localstuff$exclusivedetectors <- c('single','multi','polygonX','transectX')
## 'signal' is not a count detector 2011-02-01
.localstuff$countdetectors <- c('count','polygon','transect','unmarked','telemetry')
.localstuff$detectors3D <- c('proximity','count','signal','signalnoise','polygon',
                             'transect','times','cue','unmarked','presence','telemetry')
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
        telemetry = 13,
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
    if (detector(traps) %in% c('polygon', 'transect','polygonX', 'transectX','telemetry'))
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
############################################################################################
## moved from methods.r 2012-10-28

model.string <- function (model, userDfn) {
    if (!is.null(userDfn)) {
        if (!is.null(model$D))
            model$D <- paste('~userD', userDfn('name'), sep='.')
    }
    temp <- paste (names(model), as.character(model), collapse=' ', sep='')
    temp
}
fixed.string <- function (fixed) {
    if (is.null(fixed) | length(fixed)==0) 'none'
    else paste (names(fixed), as.character(fixed), collapse=', ', sep=' = ')
}
############################################################################################

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

nclusters <- function (capthist) {
    if (ms(capthist)) {
	lapply(capthist, nclusters)
    }
    else 	{
        nmash <- attr(capthist, 'n.mash')
        ifelse (is.null(nmash), 1, length(nmash))
    }
}
###############################################################################
## moved here from make.grid 2012-09-04

# leadingzero <- function (x) {
#     formatC(x, width=max(nchar(x)), flag="0")  ## returns character value
# }

## clunky but effective re-write 2012-09-04
leadingzero <- function (x) {
    x <- as.character(x)
    w <- max(nchar(x))
    n0 <- function(n) paste(rep('0',n), collapse='')
    paste(sapply(w-nchar(x), n0), x, sep='')
}

###############################################################################

    HN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]
        g0 * exp (-r^2 / 2 / sigma^2)
    }

    HZ <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * (1 - exp (-(r / sigma)^-z))
    }

    EX <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * exp (-r / sigma)
    }
    UN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]
        ifelse (r<=sigma, g0, 0)
    }
    CHN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * ( 1 - (1 - exp (-r^2 / 2 / sigma^2)) ^ z )
    }
    WEX <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
        ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
    }
    ANN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
        g0 * exp (-(r-w)^2 / 2 / sigma^2)
    }
    CLN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        CV2 <- (z/sigma)^2
        sdlog <- log(1 + CV2)^0.5
        meanlog <- log(sigma) - sdlog^2/2
        g0 * plnorm(r, meanlog, sdlog, lower.tail = FALSE)
    }
    CG <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * pgamma(r, shape=z, scale=sigma/z, lower.tail = FALSE)
    }
    CN <- function (r, pars, cutval) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        x <- z * (r - sigma)
        g0 * (1 + (1 - exp(x)) / (1 + exp(x)))/2
    }
    BSS <- function (r, pars, cutval) {
        b0 <- pars[1]; b1 <- pars[2]
        gam <- -(b0 + b1 * r);
        pnorm (gam, mean=0, sd=1, lower.tail=FALSE)
    }
    SS <- function (r, pars, cutval) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(cutval))
            stop ("require 'details$cutval' for signal strength plot")
        mu <- beta0 + beta1 * r
        1 - pnorm (q=cutval, mean=mu, sd=sdS)
    }
    SSS <- function (r, pars, cutval) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(cutval))
            stop ("require 'details$cutval' for signal strength plot")
        ## spherical so assume distance r measured from 1 m
        mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
        mu[r<1] <- beta0
        1 - pnorm (q=cutval, mean=mu, sd=sdS)
    }
    SN <- function (r, pars, cutval) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
        muN <- pars[4]; sdN <- pars[5]
        muS <- beta0 + beta1 * r
        1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
    }

    SNS <- function (r, pars, cutval) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3];
        muN <- pars[4]; sdN <- pars[5]
        ## spherical so assume distance r measured from 1 m
        muS <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
        muS[r<1] <- beta0
        1 - pnorm (q=cutval, mean=muS-muN, sd=sqrt(sdS^2+sdN^2))
    }
############################################################################################

gradient <- function (pars, fun, eps=0.001, ...)
## quick & dirty 2009 09 14
## used by plot.secr for delta method limits
{
  est <- pars
  g   <- pars
  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- fun (est, ...)
      est[i]  <- temp + delta
      fplus   <- fun (est, ...)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
}
############################################################################################

distancetotrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coor in col 2
    X <- matrix(unlist(X), ncol = 2)
    nxy <- nrow(X)
    ## 2011-10-14
    detecttype <- detector(traps)
    detecttype <- ifelse (is.null(detecttype), "", detecttype)
    if (detecttype %in% .localstuff$polydetectors) {
        ## approximate only
        traps <- split(traps, polyID(traps))
        trpi <- function (i, n=100) {
            intrp <- function (j) {
                tmp <- data.frame(traps[[i]][j:(j+1),])
                if (tmp$x[1] == tmp$x[2])
                    data.frame(x=rep(tmp$x[1], n),
                               y=seq(tmp$y[1], tmp$y[2], length=n))
                else
                    data.frame(approx(tmp, n=n))
            }
            tmp <- lapply(1:(nrow(traps[[i]])-1),intrp)
            do.call(rbind, tmp)
        }
        trps <- do.call(rbind, lapply(1:length(traps), trpi))
        trps <- matrix(unlist(trps), ncol = 2)
    }
    else
        trps <- traps
    temp <- .C("nearest",  PACKAGE = 'secr',
        as.integer(nxy),
        as.double(X),
        as.integer(nrow(trps)),
        as.double(unlist(trps)),
        index = integer(nxy),
        distance = double(nxy)
    )
    if (detecttype %in% c('polygon', 'polygonX')) {
        inside <- lapply(traps, pointsInPolygon, xy=X)
        inside <- do.call(rbind, inside)
        temp$distance [apply(inside,2,any)] <- 0
    }
    temp$distance
}

nearesttrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2
    X <- matrix(unlist(X), ncol = 2)
    nxy <- nrow(X)
    temp <- .C("nearest",  PACKAGE = 'secr',
        as.integer(nxy),
        as.double(X),
        as.integer(nrow(traps)),
        as.double(unlist(traps)),
        index = integer(nxy),
        distance = double(nxy)
    )
    temp$index
}

transform <- function (x, link) {
  switch (link,
          identity = x,
          log = log(x),
          neglog = log(-x),
          logit = logit(x),
          odds = odds(x),
          sin = sine(x)
  )
}

untransform <- function (beta, link) {
  switch (link,
          identity = beta,
          log = exp(beta),
          neglog = -exp(beta),
          logit = invlogit(beta),
          odds = invodds(beta),
          sin = invsine(beta))
}

se.untransform <- function (beta, sebeta, link) {
  switch (link,
          identity = sebeta,
          log = exp(beta) * sqrt(exp(sebeta^2)-1),
          neglog = exp(beta) * sqrt(exp(sebeta^2)-1),
          logit = invlogit(beta) * (1-invlogit(beta)) * sebeta,
          sin = NA)         ####!!!!
}


mlogit.untransform <- function (beta, mix) {
    nmix <- max(mix)    ## assume zero-based
    ##  b <- beta[match(2:nmix, mix)] -- 2010 02 26
    b <- beta[2:nmix]    ## 2010 02 26
    pmix <- numeric(nmix)
    pmix[2:nmix] <- exp(b) / (1+sum(exp(b)))
    pmix[1] <- 1 - sum(pmix[2:nmix])
    pmix[mix]   ## same length as input
}

# vector version of transform()
Xtransform <- function (real, linkfn, varnames) {
  out <- real
  for (i in 1:length(real)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = real[i],
                  log = log(real[i]),
                  neglog = log(-real[i]),
                  logit = logit(real[i]),
                  odds = odds(real[i]),
                  sin = sine(real[i]))
  }
  out
}
se.Xtransform <- function (real, sereal, linkfn, varnames) {
  out <- real
  for (i in 1:length(real)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = sereal[i],
                  log = log((sereal[i]/real[i])^2 + 1)^0.5,
                  neglog = log((sereal[i]/-real[i])^2 + 1)^0.5,
                  logit = sereal[i] / real[i] / (1 - real[i]),
                  sin = NA)
  }
  out
}

# vector version of untransform()
Xuntransform <- function (beta, linkfn, varnames) {
  out <- beta
  for (i in 1:length(beta)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = beta[i],
                  log = exp(beta[i]),
                  neglog = -exp(beta[i]),
                  logit = invlogit(beta[i]),
                  odds = invodds(beta[i]),
                  sin = invsine(beta[i]))
  }
  out
}

se.Xuntransform <- function (beta, sebeta, linkfn, varnames)
# Approximate translation of SE to untransformed scale
# Delta method cf Lebreton et al 1992 p 77
{
  out <- beta
  if (length(beta)!=length(sebeta))
      stop ("'beta' and 'sebeta' do not match")
  if (!all(varnames %in% names(linkfn)))
      stop ("'linkfn' component missing for at least one real variable")
  for (i in 1:length(beta)) {
      vn <- varnames[i]
      out[i] <- switch (linkfn[[vn]],
                  identity = sebeta[i],
                  log = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                  neglog = exp(beta[i]) * sqrt(exp(sebeta[i]^2)-1),
                  logit = invlogit(beta[i]) * (1-invlogit(beta[i])) * sebeta[i],
                  sin = NA)         ####!!!!
  }
  out
}

## End of miscellaneous functions
############################################################################################

group.levels <- function (capthist, groups, sep='.') {
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.levels, groups, sep)   ## sep added 2010 02 24
        sort(unique(unlist(temp)))  ## vector of global levels
    }
    else {
        if (is.null(groups)) 0
        else {
            temp <- as.data.frame(covariates(capthist)[,groups])
            if (ncol(temp) != length(groups))
                stop ("one or more grouping variables is missing ",
                      "from covariates(capthist)")
            sort(levels(interaction(temp, drop=T, sep=sep)))  # omit null combinations, sort as with default of factor levels
        }
    }
}
############################################################################################

n.occasion <- function (capthist) {
## return the number of sampling occasions for each session in capthist
    if (inherits(capthist, 'list')) {
        sapply(capthist, n.occasion)
    }
    else {
        ncol(capthist)
    }
}

############################################################################################

group.factor <- function (capthist, groups, sep='.')
## convert a set of grouping factors to a single factor (g)
## levels common to all sessions
{
    if (inherits(capthist, 'list')) {
        temp <- lapply(capthist, group.factor, groups)  ## recursive call
        grouplevels <- group.levels(capthist, groups)
        if (length(grouplevels)<2)
            temp
        else
            # list; force shared factor levels on each component
            lapply (temp, factor, levels=grouplevels)
    }
    else {
        if (is.null(groups) | (length(groups)==0) )
            return (factor(rep(1, nrow(capthist))))
        temp <- as.data.frame(covariates(capthist)[,groups])
        if (ncol(temp) != length(groups))
            stop ("one or more grouping variables is missing from ",
                  "covariates(capthist)")
        temp <- interaction(temp, drop=T, sep=sep)  # omit null combinations
        temp
    }
}
############################################################################################
