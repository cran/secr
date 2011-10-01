###############################################################################
## package 'secr'
## functions.R
## miscellaneous functions
## last changed
## 2009 12 10 mixtures
## 2010 02 05 update signal-strength detection functions
## 2010 02 25 insertdim factor handling improved to retain ordering of levels
## 2010-10-09 allow detectfn 6,7
## 2010-12-02
## 2011-01-04 param=1 GR parameterisation in secrloglik
## 2011-01-23 distancetotrap updated for polygons
## 2011-02-06 distancetotrap updated for polygonX
## 2011-03-19 shift logmultinom to secrloglik to allow depends on session detector
## 2011-09-26 new detector type checks
## 2011-09-26 experimental presence detector
###############################################################################

###############################################################################

## distancetotrap <- function (X, traps) {
##     ## X should be 2-column dataframe, mask, matrix or similar
##     ## with x coord in col 1 and y coor in col 2
##     X <- matrix(unlist(X), ncol = 2)
##     disttotrap <- function (xy) {
##         temp <- .C("nearest",  PACKAGE = 'secr',
##         as.double(xy),
##         as.integer(nrow(traps)),
##         as.double(unlist(traps)),
##         index = integer(1),
##         distance = double(1)
##         )
##         temp$distance
##     }
##     apply(X, 1, disttotrap)
## }
##
## nearesttrap <- function (X, traps) {
##     ## X should be 2-column dataframe, mask, matrix or similar
##     ## with x coord in col 1 and y coor in col 2
##     X <- matrix(unlist(X), ncol = 2)
##     nearest <- function (xy) {
##         temp <- .C("nearest",  PACKAGE = 'secr',
##         as.double(xy),
##         as.integer(nrow(traps)),
##         as.double(unlist(traps)),
##         index = integer(1),
##         distance = double(1)
##         )
##         temp$index
##     }
##     apply(X, 1, nearest)
## }

distancetotrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coor in col 2
    X <- matrix(unlist(X), ncol = 2)
    nxy <- nrow(X)
    if (detector(traps) %in% .localstuff$polydetectors) {
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
    if (detector(traps) %in% c('polygon', 'polygonX')) {
        inside <- lapply(traps, insidepoly, xy=X)
        inside <- do.call(rbind, inside)
        temp$distance [apply(inside,2,any)] <- 0
    }
    temp$distance
}

nearesttrap <- function (X, traps) {
    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coor in col 2
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

CLdensity <- function (beta, object, individuals, sessnum)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
# Return the density for given g0, sigma, z in beta
# Only 1 session
{
    sum(1 / esa (object, sessnum, beta)[individuals])
}
############################################################################################

CLgradient <- function (object, individuals, sessnum, eps=0.001)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
{
  beta <- object$fit$par
  if (object$detectfn %in% c(0,2,9)) {    ## halfnormal, exponential and binary SS have just 2 parameters
      est <- beta[1:2]
      g   <- double(2)
  } else {                       ## other detectfn have 3 parameters
      est <- beta[1:3]
      g   <- double(3)
  }

  est <- beta
  g   <- beta

  ## consider replacing this with packaged, optimized gradient function (nlme fnHess?)

  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- CLdensity (est, object, individuals, sessnum)
      est[i]  <- temp + delta
      fplus   <- CLdensity (est, object, individuals, sessnum)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
}

############################################################################################

CLmeanesa <- function (beta, object, individuals, sessnum, noccasions = NULL)
## object is a fitted secr object (CL=T)
## individuals is vector indexing the subset of a to be used
# Return the weighted mean esa for given g0, sigma, z in beta
# Only 1 session
{
## mean (esa (object, sessnum, beta)[individuals])
## modified 2010-11-30 after suggestion of DLB

##  noccasions = NULL added 2011-04-04

    a <- esa (object, sessnum, beta, noccasions=noccasions)[individuals]
    length(a) / sum (1/a)
}
############################################################################################

esagradient <- function (object, individuals, sessnum, noccasions = NULL, eps=0.001)
##  noccasion = NULL added 2011-04-04
{
  beta <- object$fit$par
  if (object$detectfn %in% c(0,2,9)) {    ## halfnormal, exponential and binary SS
                                          ## have just 2 parameters
      est <- beta[1:2]
      g   <- double(2)
  } else {                       ## other detectfn have 3 parameters
      est <- beta[1:3]
      g   <- double(3)
  }

  est <- beta
  g   <- beta

  ## consider replacing this with fdHess from package nlme

  for (i in 1:length(est))
  {
      temp     <- est[i]
      if (temp != 0.0) delta <- eps * abs(temp)
      else             delta <- eps
      est[i]  <- temp - delta
      fminus  <- CLmeanesa (est, object, individuals, sessnum, noccasions)
      est[i]  <- temp + delta
      fplus   <- CLmeanesa (est, object, individuals, sessnum, noccasions)
      g[i]    <- (fplus - fminus) / (2.0 * delta)
      est[i]  <- temp;
  }
  g
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
        lev <- group.levels(capthist, groups)
        if (length(lev)<2) temp
        else lapply (temp, factor, levels=lev)  # list; force shared factor levels on each component
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

disinteraction <- function (capthist, groups, sep='.') {
    ngv <- length(groups)
    f <- group.levels(capthist, groups, sep=sep)
    if (ngv>1)
        temp <- matrix(unlist(strsplit(as.character(f), sep, fixed=TRUE)),
                       byrow = T, ncol = ngv)
    else temp <- f
    temp <- data.frame(temp)
    names(temp) <- groups
    temp
}

############################################################################################

secr.lpredictor <- function (model, newdata, indx, beta, field, beta.vcv=NULL) {
    vars <- all.vars(model)
    if (any(!(vars %in% names(newdata))))
        stop ("one or more model covariates not found in 'newdata'")
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata),dimnames=list(NULL,c('estimate','se')))
    mat <- model.matrix(model, data=newdata)
    ## drop pmix beta0 column from design matrix (always zero)
    if (field=='pmix') {
        mat <- mat[,-1,drop=FALSE]
    }
    lpred[,1] <- mat %*% beta[indx]
    if (is.null(beta.vcv)) return ( cbind(newdata,lpred) )
    else {

        vcv <- beta.vcv[indx,indx]    ## OR maybe all betas?
        nrw <- 1:nrow(mat)
        vcv <- apply(expand.grid(nrw, nrw), 1, function(ij)
            mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F]))  # link scale
        vcv <- matrix (vcv, nrow = nrw)
        lpred[,2] <- diag(vcv)^0.5
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

############################################################################################

## work in progress
# secr.deriv.inv.link <- function (model, x, beta.vcv) {
#     vars <- all.vars(model)
#     if (any(!(vars %in% names(newdata))))
#         stop ('some model parameters not found in newdata')
#     newdata <- as.data.frame(newdata)
#     mat <- model.matrix(model, data=newdata)
#     apply(mat,1,function(x) sum(outer(x,x) * beta.vcv[indx,indx]))  # link scale
#
#
# "deriv.inverse.link" <-
# function(real,x,link)
# {
# real=as.vector(real)
# switch(link,
#
# logit=x*real*(1-real),
# log=x*real,
# identity=x,
# sin=x*cos(asin(2*real-1))/2,
# )
# }
#
# }

############################################################################################

makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        if (names(link)[i] == 'pmix') {
            mlogit.untransform(Yp, design$parameterTable[,'pmix'])
        }
        else {
            Yp <- untransform(Yp, link[[i]])
            Yp[design$parameterTable[,i]]   ## replicate as required
        }
    }

    ## construct matrix of detection parameters
    nrealpar  <- length(design$designMatrices)
    parindx$D <- NULL ## detection parameters only
    link$D    <- NULL ## detection parameters only
    detectionparameters <- names(link)
    fixed.dp <- fixed[names(fixed) %in% detectionparameters]
    if (length(fixed.dp)>0)
        for (a in names(fixed.dp))     ## bug fixed by adding this line 2011-09-28
            link[[a]] <- NULL
    if (length(link) != nrealpar)
        stop ("number of links does not match design matrices")

    temp <- sapply (1:nrealpar, modelfn)
    if (nrow(design$parameterTable)==1) temp <- t(temp)
    nrw <- nrow(temp)
    ## make new matrix and insert columns in right place
    temp2 <- as.data.frame(matrix(nrow = nrw, ncol = length(detectionparameters)))
    names(temp2) <- detectionparameters
    temp2[ , names(design$designMatrices)] <- temp          ## modelled
    if (!is.null(fixed.dp) & length(fixed.dp)>0)
    ##  temp2[ , names(fixed.dp)] <- t(sapply(fixed.dp, rep, nrw))    ## fixed
    temp2[ , names(fixed.dp)] <- sapply(fixed.dp, rep, nrw)    ## fixed
    as.matrix(temp2)
}
############################################################################################

## scale detection parameters as requested

scaled.detection <- function (realparval, scalesigma, scaleg0, D) {
    realnames <- dimnames(realparval)[[2]]
    sigmaindex <- match('sigma', realnames)
    g0index <- match('g0', realnames)
    if (scalesigma) {   ## assuming previous check that scalesigma OK...
        if (is.na(D)) sigmaindex <- NA
        if (is.na(sigmaindex))
            stop ("'scalesigma' requires both 'sigma' and 'D' in model")
        realparval[,sigmaindex]  <- realparval[,sigmaindex] / D^0.5
    }
    if (scaleg0)    {   ## assuming previous check that scaleg0 OK...
        if (is.na(g0index) | is.na(sigmaindex))
            stop ("'scaleg0' requires both 'g0' and 'sigma' in model")
        realparval[,g0index]  <- realparval[,g0index] / realparval[,sigmaindex]^2
    }
    realparval
}
###############################################################################

getD <- function (designD, beta, mask, parindx, link, fixed, MS, ngrp,
    nsession) {
    dimD <- attr(designD,'dimD')
    if (!is.null(fixed$D)) {
        if (MS) nmask <- max(sapply(mask, nrow))
        else nmask <- nrow(mask)
        dimD <- c(nmask, ngrp, nsession)
        D <- rep(fixed$D, prod(dimD))
        D <- array(D, dim=dimD )
    }
    else {
        if (is.function(designD))
            D <- designD(mask, beta[parindx$D], )
        else
            D <- designD %*% beta[parindx$D]   # linear predictor
        D <- array(untransform (D, link$D), dim = dimD )
        D[D<0] <- 0                            # silently truncate at zero
    }
    D
}
###############################################################################

secr.loglikfn <- function (beta, parindx, link, fixed, designD, design,
    design0 = NULL, capthist, mask, detectfn = 0, CL = T, groups = NULL,
    details, logmult, dig = 3, betaw = 10)

# Return the negative log likelihood for inhomogeneous Poisson spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function
#    0 = halfnormal, 1 = hazard, 2 = exponential etc.
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
#
# details$trace=T sends a one-line report to the screen

{
    MS <- ms(capthist)    ## updated 2011-01-07
    if (MS)
        sessionlevels <- session(capthist)
    else
        sessionlevels <- 1
    nsession <- length(sessionlevels)

    #--------------------------------------------------------------------
    # Groups
    grp  <- group.factor (capthist, groups)
    ngrp <- max(1,length(group.levels(capthist, groups)))

    #--------------------------------------------------------------------
    # Fixed beta
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        fb[is.na(fb)] <- beta
        beta <- fb    ## complete
    }

    #--------------------------------------------------------------------
    # Detection parameters
    realparval  <- makerealparameters (design, beta, parindx, link, fixed)
    realparval0 <- makerealparameters (design0, beta, parindx, link, fixed)

    #--------------------------------------------------------------------

    minprob <- details$minprob
    if (is.null(minprob)) minprob <- 1e-50

    #--------------------------------------------------------------------
    # Density
    D.modelled <- !CL & is.null(fixed$D)

    if (!CL) {
         D <- getD (designD, beta, mask, parindx, link, fixed, MS,
                    ngrp, nsession)
        if (sum(D) <= 0)
            stop ("invalid density")
    }
    #--------------------------------------------------------------------

    loglik <- 0

    ####################################################
    ## start loop over sessions

    for (sessnum in 1:nsession) {

        ## in multi-session case must get session-specific data from lists
        if (MS) {
            session.capthist <- capthist[[sessnum]]
            session.traps    <- traps(capthist)[[sessnum]]
            session.mask     <- mask[[sessnum]]
            session.grp      <- grp[[sessnum]]
            session.xy <- 0
        }
        else {
            session.capthist <- capthist
            session.traps    <- traps(capthist)
            session.mask     <- mask
            session.grp      <- grp
        }

        nc   <- nrow(session.capthist)
        s    <- ncol(session.capthist)
        m    <- nrow(session.mask)
        cell <- attr(session.mask,'area')
        session.mask <- as.matrix(session.mask[,1:2])

        dettype <- detectorcode(session.traps)

        if (dettype %in% c(5,9)) {    # cue or signal strength
            session.signal <- signal(session.capthist)
            session.signal <- switch( details$tx,
                log = log(session.signal),
                logit = logit(session.signal),
                identity = session.signal
            )
        }
        else
            session.signal <- 0

        if (dettype == 9)   ## cue
            cuerate <- exp(beta[max(unlist(parindx))+1])   ## fudge: last
        else
            cuerate <- 1 ## dummy

        if (dettype %in% c(3,6)) {
            k <- table(polyID(session.traps))
            K <- length(k)
            k <- c(k,0)   ## zero terminate
            session.xy <- xy(session.capthist)
        }
        else
        if (dettype %in% c(4,7)) {
            k <- table(transectID(session.traps))
            K <- length(k)
            k <- c(k,0) ## zero terminate
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
        ## differentiate so density & g do not both need to use sessions
        if (CL)
            density <- 0
        else
            density <- D[1:m,,min(dim(D)[3],sessnum)]

        sessg <- min (sessnum, design$R)

        #------------------------------------------
        # allow for scaling of detection parameters

        Dtemp <- ifelse (D.modelled, D[1,1,sessnum], NA)
        Xrealparval <- scaled.detection (realparval, details$scalesigma, details$scaleg0, Dtemp)
        Xrealparval0 <- scaled.detection (realparval0, details$scalesigma, details$scaleg0, Dtemp)

        #-----------------------------------------
        # check valid parameter values
        if (!all(is.finite(Xrealparval))) {
            cat ('beta vector :', beta, '\n')
            warning ("extreme 'beta' in 'secr.loglikfn' ",
                "(try smaller stepmax in nlm Newton-Raphson?)")
            return (1e10)
        }
        if (!all(is.finite(density))) {
            cat ('densities :', head(density), '\n')
            warning ("bad densities in 'secr.loglikfn' ",
                "(try different optimisation method, link, or model?")
            return (1e10)
        }

        #------------------------------------------
        ##################################
        # debugging
        # cat('dettype = ', dettype, '\n')
        # print(Xrealparval)
        # print(summary(session.capthist))
        # cat('Dim\n')
        # cat ('nc = ', nc, '\n')
        # cat ('k = ', k, '\n')
        # cat ('m = ', m, '\n')
        # cat ('s = ', s, '\n')
        # cat ('ngrp = ', ngrp, '\n')
        # cat ('session.grp = ', session.grp, '\n')
        # cat('Mask\n')
        # print(summary(session.mask))
        # PIA = parameter index array
        # cat ('design$PIA','\n')
        # print(table(design$PIA[sessg,1:nc,1:s,1:k]))
        # cat ('design0$PIA','\n')
        # print(design0$PIA)
        # print(table(design0$PIA[sessg,max(1,ngrp),1:s,1:k]))
        # cat ('sessD = ', sessD, '\n')
        # cat ('table(D[1:m,,sessD])\n')
        # print (table(D[1:m,,sessD]))
        # cat ('sessg = ', sessg, '\n')
        # cat ('cell  = ', cell, '\n')
        # cat ('cuerate  = ', cuerate, '\n')
        # cat ('CL = ', CL, '\n')
        # cat ('density ', density, '\n')
        # stop('Debug stop')
        ##################################
        ## For conditional likelihood, supply a value for each
        ##  animal, not just groups

        K <- ifelse (detector(session.traps) %in% .localstuff$polydetectors, length(k)-1, k)
        if (CL) tempPIA0 <- design0$PIA[sessg,1:nc,1:s,1:K, ]
        else tempPIA0 <- design0$PIA[sessg,1:ngrp,1:s,1:K, ]     ## drop=FALSE unnecessary?
        if (is.numeric(details$distribution)) {
            if (details$distribution < nc)
                stop ("superpopulation (details$distribution) ",
                      "is less than number observed")
            distrib <- details$distribution
        }
        else
            distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)

        #-------------------------------------------
        # experimental 'unmarked' detector type
        #-------------------------------------------
        if (dettype == 10) {
            ## does not allow within-session models
            index <- ifelse(nrow(Xrealparval)==1, 1, sessnum)
            g0 <- Xrealparval[index,1]
            sigma <- Xrealparval[index,2]
            z <- ifelse (detectfn %in% c(1,3,5,6,7,8), Xrealparval[index,3], 1)
            temp <- .C('unmarkedloglik', PACKAGE = 'secr',
                as.integer(session.capthist), ## n,s,k array
                as.integer(nrow(session.capthist)),
                as.integer(s),          ## number of sampling occasions
                as.integer(K),          ## number of detectors
                as.double(density[1]),  ## Parameter value
                as.double(g0),          ## Parameter value
                as.double(sigma),       ## Parameter value
                as.double(z),           ## Parameter value
                as.integer(detectfn),
                as.integer(0),
                value = double(1),
                resultcode = integer(1)
            )
        }
        #-------------------------------------------
        # experimental 'presence' detector type
        #-------------------------------------------
        else if (dettype == 11) {
            ## details$presence sets type, which may be
            ## simple            Royle-Nichols
            ## integrated
            ## pairwise (not this version)

            if (is.null(details$presence))
                details$presence <- 'integrated'
            if (!(details$presence %in% c('simple', 'pairwise', 'integrated')))
                stop("details$presence should be one of ",
                     "'simple', 'pairwise', 'integrated'")

            ## does not allow within-session models
            index <- ifelse(nrow(Xrealparval)==1, 1, sessnum)
            g0 <- Xrealparval[index,1]
            sigma <- Xrealparval[index,2]
            z <- ifelse (detectfn %in% c(1,3,5,6,7,8), Xrealparval[index,3], 1)

            if (details$presence == 'pairwise') {
                stop("'pairwise' temporarily blocked")
                # pairwise likelihood allowing for overlap
                # thanks to Martin Ridout
                # require (occsim)
                CH <- apply(session.capthist, 2:3, sum)>0   ## collapse individuals
                yy <- apply(CH,2,sum)                       ## vector of number at each detector yk
                param <- c(density[1], g0, sigma/100)
                dx <- session.traps[,1]/100
                dy <- session.traps[,2]/100
                temp <- .C('loglik2', PACKAGE = 'occsim',
                    as.double(param),  ## parameter vector
                    as.integer(s),     ## number of sampling occasions
                    as.double(dx),     ## x-coordinates of detectors
                    as.double(dy),     ## y-coordinates of detectors
                    as.integer(yy),    ## per detector frequencies (vector length k)
                    as.integer(K),     ## number of detectors
                    value = double(1)
                )
                temp$resultcode <- 0
            }
            else {
                if ((.localstuff$iter < 1) & (details$presence == 'simple') &
                    !(detectfn %in% c(1,4,7,8)))
                    warning ("simple presence requires detectfn for which",
                             " sigma = radius (1,4,7,8)", .call = FALSE)
                temp <- .C('presenceloglik', PACKAGE = 'secr',
                    as.integer(session.capthist), ## n,s,k array
                    as.integer(nrow(session.capthist)),
                    as.integer(s),          ## number of sampling occasions
                    as.integer(K),          ## number of detectors
                    as.double(density[1]),  ## Parameter value
                    as.double(g0),          ## Parameter value
                    as.double(sigma),       ## Parameter value
                    as.double(z),           ## Parameter value
                    as.integer(detectfn),
                    as.integer(details$presence == 'integrated'),
                    value = double(1),
                    resultcode = integer(1))
            }
        }
        #-------------------------------------------
        # typical call (not 'presence' or 'unmarked'
        #-------------------------------------------
        else {
            temp <- .C('secrloglik', PACKAGE = 'secr',
                as.integer(CL),       # 0 = full, 1 = CL
                as.integer(dettype),  # 0 = multicatch, 1 = proximity, etc
                as.integer(details$param), # 0 = Borchers & Efford, 1 Gardner & Royle
                as.integer(distrib),  # Poisson = 0 vs binomial = 1 (distribution of n)
                as.integer(session.capthist),
                as.double(unlist(session.xy)),         # polygon or transect detection locations
                as.double(session.signal),
                as.integer(session.grp),
                as.integer(nc),
                as.integer(s),
                as.integer(k), # may be zero-terminated vector for parts of polygon or transect detector
                as.integer(m),
                as.integer(ngrp),
                as.integer(details$nmix),
                as.double(trps),
                as.double(session.mask),
                as.double(density),                            # density at each mask point x ngrp cols
                as.double(Xrealparval),
                as.double(Xrealparval0),
                as.integer(nrow(Xrealparval)),                 # number of rows in lookup table
                as.integer(nrow(Xrealparval0)),                # ditto, naive
                as.integer(design$PIA[sessg,1:nc,1:s,1:K,]),   # index of nc,S,K,mix to rows in Xrealparval
                as.integer(tempPIA0),                      # index of ngrp,S,K,mix to rows in Xrealparval0
                as.double(cell),                               # mask cell area
                as.double(cuerate),                            # miscellaneous parameter
                as.integer(detectfn),
                as.integer(details$binomN),
                as.double(details$cutval),
                as.double(minprob),
                a=double(nc),
                value=double(1),
                resultcode=integer(1))
        }

        if (temp$resultcode != 0)
            loglik <- NA
        else
            loglik <- loglik + temp$value
    } ## end loop over sessions

    ####################################################

##    if (!(dettype==9))
##    loglik <- loglik + logmult   ## term calc before and saved as global variable
## replaced 2011-03-19 with...
## unclear whether this is correct wrt groups
    if (logmult & detector(session.traps) %in% .localstuff$simpledetectors)
        loglik <- loglik + logmultinom(session.capthist,
                                       group.factor(session.capthist, groups))

    .localstuff$iter <- .localstuff$iter + 1   ## moved outside loop 2011-09-28
    if (details$trace) {

        ## allow for fixed beta parameters 2009 10 19
        if (!is.null(details$fixedbeta))
            beta <- beta[is.na(details$fixedbeta)]

        cat(format(.localstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=10),
            formatC(beta, format='f', digits=dig+1, width=betaw),
        '\n')

        flush.console()
    }

   if (is.finite(loglik)) -loglik   # return the negative loglikelihood
   else 1e10

}
############################################################################################

## NOT WORKING FOR A LONG TIME
MRsecr.loglikfn <- function (beta, parindx, link, fixed, designD, design, design0 = NULL,
                           capthist, mask, detectfn = 0, CL = F, groups = NULL,
                           details, logmult, dig = 3, betaw = 10)

# Return the negative log likelihood for inhomogeneous Poisson spatial capture-recapture model
# Mark-resight version

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
#
# trace=T sends a one-line report to the screen

{

    MS <- inherits(capthist,'list')
    if (MS) sessionlevels <- session(capthist)
    else sessionlevels <- 1
    nsession <- length(sessionlevels)

    if (!detectfn %in% c(0:3,5,6,7,8))
        stop ("'detectfn' can only take values ",
              "0 (halfnormal), ",
              "1 (hazard-rate), ",
              "2 (exponential), ",
              "3 (compound halfnormal, ",
              "5 (w-exponential), ",
              "6 (annular normal), ",
              "7 (cumulative lognormal), ",
              "8 (cumulative gamma) for now")

    #--------------------------------------------------------------------
    # Groups
    grp  <- group.factor (capthist, groups)

    ngrp <- max(1,length(group.levels(capthist, groups)))
    #--------------------------------------------------------------------
    # Detection parameters

    realparval  <- makerealparameters (design, beta, parindx, link, fixed)
    realparval0 <- makerealparameters (design0, beta, parindx, link, fixed)  # naive

    ## place holders for now 2009 09 27, 2009 10 25
    binomN <- 0

    #--------------------------------------------------------------------
    # Density

    D.modelled <- is.null(fixed$D)
    D <- getD (designD, beta, mask, parindx, link, fixed, MS, ngrp, nsession)
    if (sum(D)<=0)
        stop ("invalid density")

    #--------------------------------------------------------------------

    loglik <- 0

    ####################################################
    ## start loop over sessions

    for (sessnum in 1:nsession) {

        ## in multi-session case must get session-specific data from lists
        if (MS) {
            session.capthist <- capthist[[sessnum]]
            session.traps    <- traps(capthist)[[sessnum]]
            session.mask     <- mask[[sessnum]]
            session.grp      <- grp[[sessnum]]
        }
        else {
            session.capthist <- capthist
            session.traps    <- traps(capthist)
            session.mask     <- mask
            session.grp      <- grp
        }

        nc   <- nrow(session.capthist)
        s    <- ncol(session.capthist)
        m    <- nrow(session.mask)
        q    <- attr(session.capthist, 'q')
        if (is.null(q)) q <- 1
        cell <- attr(session.mask,'area')
        session.mask <- as.matrix(session.mask[,1:2])
## OUT-DATED
        session.detector <- detector(session.traps)
        if (session.detector %in% c('multi', 'single')) dettype <- 0
        else if (session.detector == 'proximity') dettype <- 1
        else if (session.detector == 'signal') dettype <- 2
        else if (session.detector == 'count') dettype <- 3
        else if (session.detector == 'areasearch') dettype <- 4
        else stop ("unrecognised detector type")

        k     <- nrow(session.traps)
        trps  <- unlist(session.traps, use.names=F)

        ###############################################################################
        ## Specific to mark-resight

        session.Tu  <- attr(session.capthist, 'Tu')[,-(1:q), drop=FALSE]
        if (is.null(session.Tu))
            session.Tu <- -1  ## marked nonID not recorded
        session.Tm  <- attr(session.capthist, 'Tm')[,-(1:q), drop=FALSE]
        if (is.null(session.Tm)) {
            session.Tm <- -1  ## marked nonID not recorded
            pID <- 1          ## for the main likelihood component
        }
        else pID <- realparval[1,'pID']  ## assumed constant for now


        ###############################################################################

        #--------------------------------
        density <- D[1:m,,min(dim(D)[3],sessnum)]
        sessg <- min (sessnum, design$R)

        #------------------------------------------
        # allow for scaling of detection parameters

        Xrealparval  <- realparval
        Xrealparval0 <- realparval0
        sigmaindex <- 2
        g0index <- 1
        if (details$scalesigma) {   ## assuming previous check that scalesigma OK...
            Xrealparval[,sigmaindex]  <- Xrealparval[,sigmaindex] / D[1,1,sessnum]^0.5
            Xrealparval0[,sigmaindex] <- Xrealparval0[,sigmaindex] / D[1,1,sessnum]^0.5
        }
        if (details$scaleg0)    {   ## assuming previous check that scaleg0 OK...
            Xrealparval[,g0index]  <- Xrealparval[,g0index] / Xrealparval[,sigmaindex]^2
            Xrealparval0[,g0index] <- Xrealparval0[,g0index] / Xrealparval0[,sigmaindex]^2
        }
        #------------------------------------------

        ##################################
        # debugging
        # cat('dettype = ', dettype, '\n')
        # print(realparval)
        # print(summary(session.capthist))
        # cat('Dim\n')
        # cat ('nc = ', nc, '\n')
        # cat ('k = ', k, '\n')
        # cat ('m = ', m, '\n')
        # cat ('s = ', s, '\n')
        # cat ('q = ', q, '\n')
        # cat ('ngrp = ', ngrp, '\n')
        # cat ('session.grp = ', session.grp, '\n')
        # cat('Mask\n')
        # print(summary(session.mask))
        ## PIA = parameter index array
        # cat ('design$PIA','\n')
        # print(table(design$PIA[sessg,1:nc,1:s,1:k]))
        ## cat ('design0$PIA','\n')
        ## print(design0$PIA)
        # print(table(design0$PIA[sessg,max(1,ngrp),1:s,1:k]))
        # cat ('sessD = ', sessD, '\n')
        # cat ('table(D[1:m,,sessD])\n')
        # print (table(D[1:m,,sessD]))
        # cat ('sessg = ', sessg, '\n')
        # cat ('cell  = ', cell, '\n')
        # cat ('density ', density, '\n')
        # stop('Debug stop')
        ##################################

        ## For conditional likelihood, supply a value for each
        ##  animal, not just groups
        tempPIA0 <- design0$PIA[sessg,1:ngrp,1:s,1:k]

        ## 2009 08 21
        distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)
        temp <- .C('MRsecrloglik', ## PACKAGE = 'secr',
            as.integer(dettype),  # 0 = multicatch, 1 = proximity, 2 = signal, 3 = count etc.
            as.integer(distrib),  # Poisson = 0 vs binomial = 1 (distribution of n)
            as.integer(session.capthist),
            as.integer(session.Tu),
            as.integer(session.Tm),
            as.integer(session.grp),
            as.integer(nc),
            as.integer(s),
            as.integer(q),
            as.integer(k),
            as.integer(m),
            as.integer(ngrp),
            as.double(trps),
            as.double(session.mask),
            as.double(density),                                   # density at each mask point x ngrp cols
            as.double(Xrealparval),
            as.double(Xrealparval0),
            as.integer(nrow(Xrealparval)),                        # number of rows in lookup table
            as.integer(nrow(Xrealparval0)),                       # ditto, naive
            as.integer(design$PIA[sessg,1:nc,1:s,1:k]),           # index of nc,S,K to rows in Xrealparval
            as.integer(tempPIA0),                                 # index of ngrp,S,K to rows in Xrealparval0
            as.double(pID),                                       # AVOID CHANGING C CODE FOR NOW 2009 10 07
            as.double(cell),
            as.integer(detectfn),
            as.integer(binomN),
            as.double(details$minprob),
            value=double(1),
            resultcode=integer(1))

        if (temp$resultcode != 0)
            stop ("error in external function 'secrloglik'")

        loglik <- loglik + temp$value
    }
    ## end loop over sessions
    ####################################################

    loglik <- loglik + logmult   ## term calc before and saved as global variable

    if (details$trace) {
        .localstuff$iter <- .localstuff$iter + 1
        cat(format(.localstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=10),
            formatC(beta, format='f', digits=dig+1, width=betaw),
        '\n')
        flush.console()
    }

   if (is.finite(loglik)) -loglik   # return the negative loglikelihood
   else 1e10

}
############################################################################################


make.lookup <- function (tempmat) {

    ## should add something to protect make.lookup from bad data...
    nrw <- nrow(tempmat)
    ncl <- ncol(tempmat)

    temp <- .C("makelookup",  PACKAGE='secr',
        as.double(tempmat),
        as.integer(nrw),
        as.integer(ncl),
        unique = integer(1),
        y      = double(nrw * ncl),
        index  = integer(nrw),
        result = integer(1))

    if (temp$result != 0)
        stop ("error in external function 'makelookup'; ",
              "perhaps problem is too large")
    lookup <- matrix(temp$y[1:(ncl*temp$unique)], nrow = temp$unique, byrow = T)

    colnames(lookup) <- colnames(tempmat)
    list (lookup=lookup, index=temp$index)
}
###############################################################################
