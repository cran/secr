fxi.secr <- function (object, i = 1, sessnum = 1, X, normal = TRUE) {

# Return scaled Pr(wi|X) for one nominated detection history,
# where X holds coordinates of points

    MS <- ms(object)
    if (MS) sessionlevels <- session(object$capthist)   ## was names(capthist) 2009 08 15
    else sessionlevels <- 1

    X <- matrix(unlist(X), ncol = 2)
    beta <- coef(object)$beta
    details <- object$details
    if (is.null(details$param)) details$param <- 0
    realparval  <- makerealparameters (object$design, beta, object$parindx,
                                       object$link, object$fixed)

    #--------------------------------------------------------------------

    D.modelled <- !object$CL & is.null(object$fixed$D)
    if (D.modelled)
        if (object$model$D != ~1)
            stop ("fxi.secr requires uniform density")
    #--------------------------------------------------------------------

    ####################################################

    ## in multi-session case must get session-specific data from lists
    if (MS) {
        session.capthist <- object$capthist[[sessnum]]
        session.traps    <- traps(object$capthist)[[sessnum]]
        session.mask     <- object$mask[[sessnum]]
        session.xy       <- 0
    }
    else {
        session.capthist <- object$capthist
        session.traps    <- traps(object$capthist)
        session.mask     <- object$mask
    }

    #--------------------------------------------------------------------

    used <- usage(session.traps)
    if (!is.null(used)) {
        if (sum (used) != (ndetector(session.traps) * ncol(used)))
            stop ("fxi.secr is not implemented for incomplete usage")
    }

    if (nrow(session.capthist)==0)
        stop ("no data for session ", sessnum)

    if (is.character(i))
        i <- match(i, row.names(session.capthist))
    if (is.na(i) | (i<1) | (i>nrow(session.capthist)))
        stop ("invalid i in fxi.secr")

    nc   <- nrow(session.capthist)
    s    <- ncol(session.capthist)
    m    <- nrow(session.mask)
    cell <- attr(session.mask,'area')
    session.mask <- as.matrix(session.mask[,1:2])

    dettype <- detectorcode(session.traps)

    if (dettype == 9) {
    # Groups defined by 'animal' covariate
        grp  <- group.factor (session.capthist, 'animal')
        ngrp <- max(1,length(group.levels(session.capthist, 'animal')))
    }
    else {
        grp <- rep(1,nrow(session.capthist))
        ngrp <- 1
    }

    if ((dettype == 5) | (dettype == 9)) {    # signal strength
        session.signal <- signal(session.capthist)
        session.signal <- switch( details$tx,
            log = log(session.signal),
            logit = logit(session.signal),
            identity = session.signal
        )
        if (object$detectfn == 11)
            stop ("fxi.secr does not work with spherical spreading in 2.0")
    }
    else
        session.signal <- 0

    if (dettype %in% c(3,4,6,7)) {
        k <- table(polyID(session.traps))
        K <- length(k)
        k <- c(k,0)
        session.xy <- xy(session.capthist)
    }
    else {
        k <- nrow(session.traps)
        K <- k
        session.xy <- 0
    }

    if (dettype == 8) {    # times -- phony use of 'signal'
        session.signal <- times(session.capthist)
    }

    trps  <- unlist(session.traps, use.names=F)

    #------------------------------------------
    # allow for scaling of detection parameters

    D.modelled <- FALSE   ## uniform
    Dtemp <- ifelse (D.modelled, object$D[1,1,sessnum], NA)
    Xrealparval <- scaled.detection (realparval, details$scalesigma, details$scaleg0, Dtemp)

    if (!all(is.finite(Xrealparval))) {
        cat ('beta vector :', beta, '\n')
        stop ("invalid parameters in 'pwuniform'")
    }

    sessg <- 1   ## max one group
    indices <- object$design$PIA[sessg,1:nc,1:s,1:K,]

    if (detector(session.traps)=='cue') {
        cuerate <- exp(coef(object)['cuerate','beta'])
    }
    else {
        cuerate <- 1
    }

# print(i)
# print(nrow(X))
# print(as.double(X))
# print(object$CL)
# print(dettype)
# print(details$param)
# print(head(session.capthist))
# print(unlist(session.xy))
# print(object$detectfn)
# cat('nc',nc,'\n')
# cat('s',s,'\n')
# cat('k',k,'\n')
# cat('m',m,'\n')
# cat('ngrp',ngrp,'\n')
# cat('details$binomN',details$binomN,'\n')
# cat('details$cutval',details$cutval,'\n')
# cat('details$minprob',details$minprob,'\n')
# cat('details$nmix',details$nmix,'\n')
# cat('normal',normal,'\n')
# cat('cuerate',cuerate,'\n')
# print(table(indices))
# print(dim(indices))
# print(Xrealparval)
# print(head(trps))
# print(cell)
# print(nrow(session.mask))

    temp <- .C('pwuniform', PACKAGE = 'secr',
        as.integer(i),                 # number of detection history within capthist
        as.integer(nrow(X)),
        as.double(X),
        as.integer(object$CL),         # 0 = full, 1 = CL
        as.integer(dettype),           # 0 = multicatch, 1 = proximity, etc
        as.integer(details$param),     # 0 = Borchers & Efford, 1 Gardner & Royle
        as.integer(session.capthist),
        as.double(unlist(session.xy)),
        as.double(session.signal),
        as.integer(grp),
        as.integer(nc),
        as.integer(s),
        as.integer(k),
        as.integer(m),
        as.integer(ngrp),
        as.integer(details$nmix),
        as.double(trps),
        as.double(session.mask),
        as.double(Xrealparval),
        as.integer(nrow(Xrealparval)), # number of rows in lookup table
        as.integer(indices),           # index of nc,S,K,mix to rows in Xrealparval
        as.double(cell),
        as.double(cuerate),
        as.integer(normal),
        as.integer(object$detectfn),
        as.integer(details$binomN),
        as.double(details$cutval),
        as.double(details$minprob),
        value=double(nrow(X)),
        resultcode=integer(1))

    if (temp$resultcode != 0)
        stop ("error in pwuniform")
     ####################################################

    temp$value
}
###############################################################################

fxi.contour <- function (object, i = 1, sessnum = 1, border = 100, nx = 64,
    levels = NULL, p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE, fitmode =
    FALSE, plotmode = FALSE, normal = TRUE, ...) {

    if (ms(object)) {
        session.traps <- traps(object$capthist)[[sessnum]]
    }
    else {
        session.traps <- traps(object$capthist)
    }


    tempmask <- make.mask (session.traps, border, nx = nx, type = 'traprect')
    xlevels <- unique(tempmask$x)
    ylevels <- unique(tempmask$y)

    fxi <- function (ni) {
        z <- fxi.secr(object, i = ni, X = tempmask, normal = normal)
        if (is.null(levels)) {
            temp <- sort(z, decreasing = T)
            levels <- approx (x = cumsum(temp)/sum(temp), y = temp, xout= p)$y
            labels <- p
        }
        else
            labels <- levels
        templines <- contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
        ## extra effort to apply correct labels
        getlevels <- function(clines.i) sapply(clines.i, function(q) q$level)
        label.levels <- function (island) {
            which.levels <- match (getlevels(island), levels)
            names(island) <- labels[which.levels]
            island
        }
        templines <- label.levels(templines)

        wh <- which.max(unlist(lapply(templines, function(y) y$level)))

        if (length(templines) > 0) {   ## 2011-04-14
            cc <- templines[[wh]]
            cc <- data.frame(cc[c('x','y')])
            templines$mode <- data.frame(x=mean(cc$x), y=mean(cc$y))
            if (fitmode)
                templines$mode <- fxi.mode(object, start=templines$mode, i = ni)
            if (plt) {
                labels <- labels[!is.na(levels)]
                levels <- levels[!is.na(levels)]
                contour (xlevels, ylevels, matrix(z, nrow = nx), add = add,
                         levels = levels, labels = labels, ...)
                if (plotmode) {
                    points(templines$mode, col = 'red', pch = 16)
                }

            }
        }
        templines
    }
    temp <- lapply(i, fxi)
    if (plt)
        invisible(temp)
    else
        temp
}

###############################################################################

fxi.mode <- function (object, i = 1, sessnum = 1, start = NULL, ...) {
    if (ms(object))
        session.capthist <- object$capthist[[sessnum]]
    else
        session.capthist <- object$capthist
    start <- unlist(start)
    if (is.null(start)) {
        session.traps <- traps(session.capthist)
        animal <- animalID(session.capthist, names=F) == i
        trp <- trap(session.capthist)[animal][1]
        start <- unlist(traps(session.capthist)[trp,])
    }
    if (is.character(i))
        i <- match(i, row.names(session.capthist))
    if (is.na(i) | (i<1) | (i>nrow(session.capthist)))
        stop ("invalid i in fxi.secr")

    fn <- function(xy,i) -fxi.secr(object, i = i, sessnum = sessnum, X = xy, normal = FALSE)
    temp <- nlm(f = fn, p = start, i = i, typsize = start, ...)$estimate
    data.frame(x=temp[1], y=temp[2])
}

###############################################################################