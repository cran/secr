###############################################################################
## package 'secr'
## secrloglik.R
## likelihood evaluation functions
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
## 2012-01-30 allow missing values in signal -- experimental
## 2012-02-07 add noise to signal if detector = signalnoise
## 2012-02-13 tweaked signalnoise
## 2012-10-19,22 telemetry
## 2012-11-12 renamed from functions.R
## 2012-11-12 neglik argument for secr.loglikfn
## 2012-11-12 'fixed' argument renamed to fixedpar to avoid conflict with maxLik
## 2012-12-17 pass usage to secrloglik
## 2013-01-08 details$ignoreusage added
## 2013-01-23 minprob passed as component of details
## 2013-07-02 param=2 esa parameterisation in secrloglik
## 2013-07-19 param=3,5 a0 parameterisation in secrloglik
## 2013-11-18 telemetry substantially improved
## 2014-03-25 param=4,5,6 sigmak parameterisation in secrloglik
## 2014-08-18 dropped 'c' column in reparameterize.sigmak
## 2014-08-19 moved secr.lpredictor to utility.R
## 2014-10-13 noneuc handled as mask-pixel-specific parameter
## 2014-12-04 improved handling of nc = 0
## 2015-05-04 reparameterize handles linear mask
## 2015-10-04 markocc
## 2015-10-16 makerealparameters bug with all fixed dp
## 2015-12-06 finalizing mark-resight additions, including chat
## 2016-10-16 makerealparameters shifted to utility.R
## 2016-10-28 userdist may be session-specific
## 2017-01-06 telemetrytype attribute of traps
## 2017-03-07 improved handling of resultcode

###############################################################################

disinteraction <- function (capthist, groups, sep='.') {
    ngv <- length(groups)
    grouplevels <- group.levels(capthist, groups, sep=sep)
    if (ngv>1)
        temp <- matrix(unlist(strsplit(as.character(grouplevels), sep, fixed=TRUE)),
                       byrow = T, ncol = ngv)
    else temp <- grouplevels
    temp <- data.frame(temp)
    names(temp) <- groups
    temp
}

############################################################################################

reparameterize.sigmak <- function (realparval, D, linear) {
    ## D,sigmak parameterisation 2014-03-12
    ## vector D must match rows in realparval
    realnames <- dimnames(realparval)[[2]]
    sigmakindex <- match('sigmak', realnames)
    cindex <- match('c', realnames)
    if (is.na(cindex))
        cval <- 0
    else
        cval <- realparval[,cindex]
    if (!('sigmak' %in% realnames))
        stop ("'param = 4:6 ' requires 'sigmak' in model")
    if (linear)
        realparval[,sigmakindex] <- realparval[,sigmakindex] / D + cval
    else
        realparval[,sigmakindex] <- realparval[,sigmakindex] / D^0.5 + cval
    dimnames(realparval)[[2]][sigmakindex] <- 'sigma'
    realparval <- realparval[, -cindex, drop = FALSE]   ## added 2014-08-18
    realparval
}
###############################################################################

reparameterize.esa <- function (realparval, mask, traps, detectfn, nocc) {

    ## esa, sigma parameterisation 2013-07-01
    ## could whole fn be coded in C?
    g0fromesa <- function (a, sigma, z, lower = 0, upper = 1) {
        fx <- function(g0) {
            ## pdot accounts for 'usage'
            ## pdot selects appropriate g0/lambda0 according to detectfn
            (sum(pdot(mask, traps, detectfn = detectfn,
                          detectpar = list(g0 = g0, lambda0 = g0, sigma = sigma, z = z),
                      noccasions = nocc)) * cell) - a
        }
        tmp <- try(uniroot(fx, lower=lower, upper=upper), silent = TRUE)
        ## debug if (inherits(tmp, 'try-error')) print(c(fx(lower),fx(upper)))
        if (inherits(tmp, 'try-error')) 0.0001 else tmp$root
    }
    cell <- getcellsize(mask)
    if (inherits(mask, 'linearmask'))
        stop ('esa parameterization not available for linear masks')
    realnames <- dimnames(realparval)[[2]]
    sigmaindex <- match('sigma', realnames)
    esaindex <- match('esa', realnames)
    z <- ifelse (ndetectpar(detectfn) == 3, realparval[, match('z', realnames)], 1)
    if (is.na(esaindex) | is.na(sigmaindex))
        stop ("'param = 2' requires both 'esa' and 'sigma' in model")
    realparval[,esaindex]  <- unlist(mapply(g0fromesa,
                                            realparval[,esaindex],  ## a
                                            realparval[,sigmaindex],
                                            z
                                            ))
    realparval
}
###############################################################################

reparameterize.a0 <- function (realparval, detectfn, linear) {
    ## a0, sigma parameterisation 2013-07-19
    realnames <- dimnames(realparval)[[2]]
    sigmaindex <- match('sigma', realnames)
    a0index <- match('a0', realnames)
    if (! all (c('a0','sigma') %in% realnames))
        stop ("'param = 3 or 5 ' requires both 'a0' and 'sigma' in model")
    if (!(detectfn %in% c(0:8, 14:19)))
        stop ('invalid combination of param = 3 or 5 and detectfn')

    if (linear)
        lambda0 <- realparval[,a0index] / realparval[,sigmaindex] * 1000
    else
        lambda0 <- realparval[,a0index] / 2 / pi / realparval[,sigmaindex]^2 * 10000
    realparval[,a0index] <- if (detectfn %in% 0:8) 1-exp(-lambda0) else lambda0
    dimnames(realparval)[[2]][a0index] <- if (detectfn<9) 'g0' else 'lambda0'
    realparval
}
###############################################################################

reparameterize <- function (realparval, detectfn, details, mask, traps, D, s) {

    ##----------------------------------------------
    ## allow for scaling of detection in one session

    ## D is scalar density or NA
    ## s is number of occasions

    linear <- inherits(mask, 'linearmask')
    if (details$param %in% 4:6) {
        ## does not allow varying density surface
        ## cf scaled.detection()
        realparval <- reparameterize.sigmak (realparval, D, linear)
    }

    if (details$param %in% c(2,6)) {
        realparval <- reparameterize.esa (realparval, mask, traps, detectfn, s)
    }
    else if (details$param %in% c(3,5)) {
        realparval <- reparameterize.a0 (realparval, detectfn, linear)
    }

    if (all(detector(traps) == 'telemetry'))
        realparval[,'lambda0'] <- 1.0
    
    realparval
}
###############################################################################

getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels, parameter = 'D') {
  ## adapted 2014-10-13 to apply to either 'D' or 'noneuc'
  if (!is.function(designD)) {
      if ((is.null(designD) | nrow(designD)==0) & (is.null(fixed[[parameter]]))) return(NULL)
  }
  if (ms(mask))
      nmask <- max(sapply(mask, nrow))
  else
      nmask <- nrow(mask)
  ngroup <- length(grouplevels)
  nsession <- length(sessionlevels)
  D <- array(dim = c(nmask, ngroup, nsession))
  dimnames(D) <- list(1:nrow(D), grouplevels, sessionlevels)
  if (!is.null(fixed[[parameter]])) {
      D[,,] <- fixed[[parameter]]
  }
  else {
      if (is.function(designD))
          D[,,] <- designD(beta[parindx[[parameter]]], mask, ngroup, nsession)
      else {
          D[,,] <- designD %*% beta[parindx[[parameter]]]   # linear predictor
          D[,,] <- untransform (D, link[[parameter]])
      }
      # silently truncate D at zero
      # allow non-positive noneuc
      if (parameter == 'D')
          D[D<0] <- 0
  }
  D
}
###############################################################################

secr.loglikfn <- function (beta, parindx, link, fixedpar, designD, designNE, design,
    design0, capthist, mask, detectfn, CL, hcov, groups, details, logmult, ncores,
    clust, dig = 3, betaw = 10, neglik = TRUE)

# Return the negative log likelihood for inhomogeneous Poisson spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function
#    0 = halfnormal, 1 = hazard, 2 = exponential etc.
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
# details$trace=T sends a one-line report to the screen

{
    ## for models fitted with old versions we need to fill in these values
    if (is.null(details$minprob)) details$minprob <- 1e-200
    if (is.null(details$telemetryscale)) details$telemetryscale <- 1
    if (is.null(details$debug)) details$debug <- 0   ## 2016-11-25
    if (is.null(details$ignoreusage)) details$ignoreusage <- FALSE  ## 2013-01-23
    if (is.null(details$unmash)) details$unmash <- FALSE ## 2013-01-23
    if (is.null(details$normalize)) details$normalize <- FALSE ## 2013-11-10
    if (is.null(details$nsim)) details$nsim <- 0 ## 2015-11-23
    if (is.null(details$chat)) details$chat <- c(1,1,1) ## 2017-04-11
    if (is.null(details$knownmarks)) details$knownmarks <- TRUE ## 2015-11-23
    if (is.null(details$splitmarked)) details$splitmarked <- FALSE ## 2016-11-23

    if (ms(capthist))
        sessionlevels <- session(capthist)
    else
        sessionlevels <- 1
    nsession <- length(sessionlevels)
    if ((ncores>1) & missing(clust))
        stop("not ready for multicore here")
    if ((ncores>1) & (nsession==1))  ## 2016-02-20
        warning("ncores>1 has no benefit when nsession = 1")

    nmix <- details$nmix

    #--------------------------------------------------------------------
    # Groups
    grp  <- group.factor (capthist, groups)
    ngroup <- max(1,length(group.levels(capthist, groups)))

    #--------------------------------------------------------------------
    # Fixed beta
    beta <- fullbeta(beta, details$fixedbeta)
    #--------------------------------------------------------------------
    # Detection parameters
    detparindx <- parindx[!(names(parindx) %in% c('D', 'noneuc'))]
    detlink <- link[!(names(link) %in% c('D', 'noneuc'))]

    realparval  <- makerealparameters (design, beta, detparindx, detlink, fixedpar)
    realparval0 <- makerealparameters (design0, beta, detparindx, detlink, fixedpar)
    
    #--------------------------------------------------------------------
    # Density
 
    D.modelled <- !CL & is.null(fixedpar$D)
    if (!CL ) {
        grplevels <- if (length(grp)>0) levels(grp[[1]]) else 1
      D <- getD (designD, beta, mask, parindx, link, fixedpar,
                 grplevels, sessionlevels, parameter = 'D')

      if (!is.na(sumD <- sum(D)))
          if (sumD <= 0)
              warning ("invalid density <= 0")
    }
    #--------------------------------------------------------------------
    # Non-Euclidean distance parameter
    NE <- getD (designNE, beta, mask, parindx, link, fixedpar,
                levels(grp[[1]]), sessionlevels, parameter = 'noneuc')
    userdistnames <- getuserdistnames(details$userdist)

    #--------------------------------------------------------------------

    ###############################################################################################
    ###############################################################################################
    sessionLL <- function (sessnum, like = 0) {

        ## function returns +logLik

        ## in multi-session case must get session-specific data from lists
        if (ms(capthist)) {
            session.capthist <- capthist[[sessnum]]
            session.traps    <- traps(capthist)[[sessnum]]
            session.mask     <- mask[[sessnum]] ## subset(mask, session = sessnum)?
            session.grp      <- grp[[sessnum]]
            session.xy       <- 0
        }
        else {
            session.capthist <- capthist
            session.traps    <- traps(session.capthist)
            session.mask     <- mask
            session.grp      <- grp
        }
        nc   <- nrow(session.capthist)
        s    <- ncol(session.capthist)
        m    <- nrow(session.mask)
        sessg <- min (sessnum, design$R)

        ## from secr 3.0 this is a vector, length s
        dettype <- detectorcode(session.traps, MLonly = TRUE, noccasions = s)
        if (all(dettype==13)) like <- 1  ## force conditional likelihood for this session only
        
        ## miscparm is used to package beta parameters that are not modelled
        ## and hence do not have a beta index specified by parindx.
        ## This includes the signal threshold and the mean and sd of noise.

        ## miscparm is passed to the C likelihood code, and also as mask attribute to userdistfn

        nmiscparm <- length(details$miscparm)
        miscparm <- numeric(max(4, nmiscparm))  ## 2015-02-21
        if (detectfn %in% c(12,13))            ## experimental signal-noise
            miscparm[1:3] <- c(details$cutval, beta[max(unlist(parindx))+(1:2)])   ## fudge: last 2
        else if (detectfn %in% c(10,11))        ## Dawson & Efford 2009 models
            miscparm[1] <- details$cutval
        else if (nmiscparm > 0)
            miscparm[1:nmiscparm] <- beta[max(unlist(parindx)) + (1:nmiscparm)]

        ###################################################################
        ## mark-resight data
        MRdata <- markresight(session.capthist, session.mask, CL, fixedpar, 
                              details$chat, sessnum, details$markresight)
        ###################################################################

        ## like == 0 is default; adjust if using conditional likelihood or allsighting
        if (CL & (like == 0)) like <- 1
        if (MRdata$allsighting) {
            if (CL) stop ("CL cannot be used with sighting-only data", call. = FALSE)
            like <- if (details$knownmarks) 5 else 6
        }
        else if (MRdata$anysighting & (like == 0) & (details$splitmarked)) {
            like <- 2
        }
        
        ##############################################
        
        ## 2013-11-10
        if ((detectfn %in% 14:19) & details$normalize) {
            if (!is.null(details$userdist))
                stop("normalization incompatible with userdist")
            miscparm[1] <- 1   ## normalize
            miscparm[3] <- 1   ## scale to mean denom 1.0 across mask
            if (!is.null(details$usecov)) {
                miscparm[2] <- 1   ## use covariate
                alpha2 <- beta[parindx[['lambda0']][2]]
                alpha2 <- ifelse (is.na(alpha2),0, alpha2)
                z <- covariates(session.mask)[,details$usecov]
                if (is.null(z))
                    stop("usecov mask covariate not valid")
                session.mask <- cbind(session.mask, exp(z * alpha2))
                session.mask <- as.matrix(session.mask[,c(1:3,3)]) ## double last col
            }
            else {
                session.mask <- as.matrix(session.mask[,c(1:2,2)])
            }
        }
        
        ## knownclass is computed from current session.capthist, so doesn't need ID
        knownclass <- getknownclass(session.capthist, nmix, hcov)
        LL <- 0
        
        # 2017-01-04 
        
        telemcode <- switch (telemetrytype(session.traps), 
            none = 0, independent = 1, dependent = 2, concurrent = 3, 0)
        
        if (all(dettype %in% c(5,12))) {    # signal strength, signalnoise
            session.signal <- signal(session.capthist)
            if (all(dettype==12))
                session.signal <- c (session.signal, noise(session.capthist))
            session.signal <- switch( details$tx,
                                      log = log(session.signal),
                                      logit = logit(session.signal),
                                      identity = session.signal
            )
            ## 2012-01-30 code missing values as negative for C code
            session.signal[is.na(session.signal)] <- -1
        }
        else {
            session.signal <- 0
        }
        
        if (all(dettype %in% c(3,6))) {    # polygonX, polygon
            k <- table(polyID(session.traps))
            K <- length(k)
            k <- c(k,0)   ## zero terminate
            session.xy <- xy(session.capthist)
        }
        else {
            if (all(dettype %in% c(4,7))) {
                k <- table(transectID(session.traps))
                K <- length(k)
                k <- c(k,0) ## zero terminate
                session.xy <- xy(session.capthist)
            }
            else {
                k <- nrow(session.traps)
                K <- k
                ## 2016-12-29
                if (any(dettype == 13)) {
                    ## ensure order matches
                    ## should have null histories in capthist
                    
                    telem <- telemetryxy(session.capthist)
                    ord <- match(names(telem), rownames(session.capthist), nomatch = 0)
                    newtelem <- vector('list',nrow(session.capthist))
                    newtelem[ord] <- telem
                    session.xy <- do.call(rbind, newtelem)
                                    }
                else {
                    session.xy <- 0
                }
            }
        }
        
        trps <- unlist(session.traps, use.names=F)
        
        usge <- usage(session.traps)
        if (is.null(usge) | details$ignoreusage) {
            usge <- matrix(1, nrow = K, ncol = s)
        }
        #---------------------------------------------------
        
        ## differentiate so density & g do not both need to use sessions
        if (like %in% c(1,3,4))
            density <- 0   ## do not model density
        else {
            density <- D[1:m,,min(dim(D)[3],sessnum)]
            ## optional scaling by session-specific number of clusters
            ## 2012-07-24
            if (details$unmash) {
                nmash <- attr(session.capthist, 'n.mash')
                if (!is.null(nmash))
                    density <- density * length(nmash)
            }
        }
        
        ##------------------------------------------
        ## allow for scaling of detection
        
        ## DOES NOT ALLOW FOR GROUP VARIATION IN DENSITY
        ## Dtemp <- if (D.modelled) D[1,1,sessnum] else NA
        ## 2014-09-30
        Dtemp <- if (D.modelled) mean(D[,1,sessnum]) else NA
        
        ## more thoughts 2015-05-05
        ## could generalize by
        ## -- making Dtemp a vector of length equal rows in realparval
        ## -- matching either
        ##      first group (as before)
        ##      sum of all groups
        ##      own group [PROBLEM: locating group of each realparval row]
        ## in all cases density is the mean over mask points
        
        ## CHECK use of Dtemp in
        ##  regionN.R
        ##  sim.secr.R
        
        ## PERHAPS for consistency make a function to construct Dtemp vector
        ## given mask, model, group matching rule (first, sum, own)
        
        Xrealparval <- reparameterize (realparval, detectfn, details,
                                       session.mask, session.traps, Dtemp, s)
        Xrealparval0 <- reparameterize (realparval0, detectfn, details,
                                        session.mask, session.traps, Dtemp, s)
        
        ##-----------------------------------------
        ## check valid parameter values
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
        
        ##--------------------------------------------------------------
        ## 2014-08-27, 2014-09-01
        ## Apply user-provided distance function or matrix
        ## do not use if detector is one of
        ## polygonX, polygon, transectX, transect, telemetry
        ## The matrix distmat is passed to the C function secrloglik, etc.
        if (is.null(details$userdist))
            distmat <- -1
        else {
            
            if (is.null(covariates(session.mask)))
                covariates(session.mask) <- data.frame(row.names = 1:m)
            if ('noneuc' %in% userdistnames)
                covariates(session.mask)$noneuc <- NE[1:m,,min(dim(NE)[3],sessnum)]
            if ('D' %in% userdistnames)
                covariates(session.mask)$D <- density
            ## pass miscellaneous unmodelled parameter(s)
            if (nmiscparm > 0)
                attr(session.mask, 'miscparm') <- miscparm
            
            distmat <- valid.userdist (details$userdist,
                                       detector(session.traps),
                                       xy1 = session.traps,
                                       xy2 = session.mask,
                                       mask = session.mask,
                                       sessnum = sessnum)
            baddist <- (!is.finite(distmat)) | (distmat<0) | is.na(distmat)
            if (any(baddist)) {
                warning ("replacing infinite, negative and NA userdist values with 1e10")
                distmat[baddist] <- 1e10
            }
        }
        
        ##--------------------------------------------------------------
        ## For conditional likelihood, supply a value for each animal,
        ## not just groups. This changes the dimensions of PIA0, so
        ## need to be wary of 'ncol' in C code.
        
        ## k-1 because we have zero-terminated these vectors
        K <- ifelse (all(detector(session.traps) %in% .localstuff$polydetectors), length(k)-1, k)
        
        PIA <- design$PIA[sessg,1:nc,1:s,1:K,]
        
        if (CL)
            PIA0 <- design0$PIA[sessg,1:nc,1:s,1:K, ]   ## CL
        else
            PIA0 <- design0$PIA[sessg,1:ngroup,1:s,1:K, ]     ## drop=FALSE unnecessary?
        if (is.numeric(details$distribution)) {
            if (details$distribution < nc)
                stop ("superpopulation (details$distribution) ",
                      "is less than number observed")
            distrib <- details$distribution
        }
        else
            distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)
        
        ##-------------------------------------------
        ## experimental 'unmarked' detector type
        ##-------------------------------------------
        if (all(dettype == 10)) {
            ## does not allow within-session models
            index <- ifelse(nrow(Xrealparval)==1, 1, sessnum)
            g0 <- Xrealparval[index,1]
            sigma <- Xrealparval[index,2]
            z <- ifelse (detectfn %in% c(1,3,5,6,7,8), Xrealparval[index,3], 1)
            temp <- .C('unmarkedloglik', # PACKAGE = 'secr',
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
        ##-------------------------------------------
        ## experimental 'presence' detector type
        ##-------------------------------------------
        else if (all(dettype == 11) & is.null(MRdata$Tn)) {
            ## details$presence sets type, which may be
            ## simple      (Royle-Nichols)
            ## integrated
            ## pairwise
            ## summedhazard (2017)
            
            if (is.null(details$presence)) details$presence <- 'summedhazard'
            details$presence <- match.arg(details$presence, c('summedhazard','simple', 'pairwise', 'integrated'))
            
            ## does not allow within-session models
            index <- ifelse(nrow(Xrealparval)==1, 1, sessnum)
            g0 <- Xrealparval[index,1]
            lambda0 <- Xrealparval[index,1]
            sigma <- Xrealparval[index,2]
            z <- ifelse (detectfn %in% c(1,3,5,6,7,8), Xrealparval[index,3], 1)
            if ((.localstuff$iter < 1) & (details$presence %in% c('simple','pairwise')) &
                !(detectfn %in% c(1,4,7,8)))
                warning ("simple presence requires detectfn for which",
                         " sigma = radius (4 or possibly 1,7,8)")
            type <- switch (details$presence, simple = 0, integrated = 1, pairwise = 2, summedhazard = 3, 4)
            if (type == 3)   ## summedhazard
                temp <- .C('presenceloglik2017', # PACKAGE = 'secr',
                           as.integer(session.capthist), ## n,s,k array
                           as.integer(nrow(session.capthist)),
                           as.integer(s),          ## number of sampling occasions
                           as.integer(K),          ## number of detectors
                           as.integer(m),
                           as.double(trps),        ## detector locations
                           as.double(unlist(session.mask)),
                           as.double(getcellsize(session.mask)),
                           as.double(density[1]),  ## Parameter value
                           as.double(lambda0),          ## Parameter value
                           as.double(sigma),       ## Parameter value
                           as.double(z),           ## Parameter value
                           as.integer(detectfn),
                           value = double(1),
                           resultcode = integer(1)
                )
            else 
                temp <- .C('presenceloglik', # PACKAGE = 'secr',
                           as.integer(session.capthist), ## n,s,k array
                           as.integer(nrow(session.capthist)),
                           as.integer(s),          ## number of sampling occasions
                           as.integer(K),          ## number of detectors
                           as.double(trps),        ## detector locations
                           as.double(density[1]),  ## Parameter value
                           as.double(g0),          ## Parameter value
                           as.double(sigma),       ## Parameter value
                           as.double(z),           ## Parameter value
                           as.integer(detectfn),
                           as.integer(type),
                           value = double(1),
                           resultcode = integer(1)
                )
        }
        ##-------------------------------------------
        ## 'telemetry' detector type i.e. sigma and related parameters only
        ## likelihood component L_T
        ## revised 2013-11-18 to allow covariates etc.
        ##-------------------------------------------
        # else if (all(dettype == 13)) {
        #     PIA <- design$PIA[sessnum,ID,,,,drop = FALSE]
        #     dim(PIA) <- dim(PIA)[-1]  ## drop session dimension
        #     temp <- telemetry.LT(session.capthist, detectfn, Xrealparval,
        #         PIA, nmix, knownclass, uppersigma = 20)
        # 
        # }
        ##--------------------------------------------
        ## typical call (not 'presence' or 'unmarked')
        ##--------------------------------------------
        else {
            if (usge[1]==0 & nmix>1)
                stop ("mixture models fail when the first detector is not ","
                        used on the first day")
            
            #######################################################################
            ## option to estimate sighting overdispersion by simulation and exit */
            if (!is.null(details$nsim)) {
                if (details$nsim > 0) {
                    
                    if (like == 1)
                        stop("simulation for overdispersion requires full likelihood (not CL)")
                    
                    temp <- .C('chat', # PACKAGE = 'secr',
                               as.integer(like),
                               as.integer(dettype),
                               as.integer(distrib),
                               as.integer(session.capthist),
                               as.integer(session.grp),
                               as.integer(nc),
                               as.integer(s),
                               as.integer(k),
                               as.integer(m),
                               as.integer(ngroup),
                               as.integer(nmix),
                               as.integer(knownclass),
                               as.double(trps),
                               as.double(distmat),
                               as.double(usge),
                               as.integer(MRdata$markocc),
                               as.double(MRdata$pi.mask),
                               as.double(unlist(session.mask)),
                               as.double(density),
                               as.double(Xrealparval0),
                               as.integer(nrow(Xrealparval0)),
                               as.integer(PIA0),
                               as.double(getcellsize(session.mask)),
                               as.double(miscparm),
                               as.integer(detectfn),
                               as.integer(expandbinomN(details$binomN, dettype)),
                               as.integer (details$nsim),
                               chat = double(3),
                               resultcode = integer(1))
                    if (temp$resultcode == 0) {
                        names(chat) <- c('Tu','Tm','Tn')
                        return(temp$chat)
                    }
                    else  {
                        warning ("chat calculation failed, resultcode = ", temp$resultcode)
                        return (1)
                    }
                }
            }
            #####################################################################
            if (details$debug>=3) browser()
          
            temp <- .C('secrloglik', # PACKAGE = 'secr',
                       as.integer(like),          # 0 = full, 1 = CL, 2 = full, splitmarked,
                       # 3 = concurrent, 4 = concurrent  ????
                       # CL, 5 = sighting only known n0, 6 = sighting only, unknown n0
                       as.integer(dettype),       # 0 = multicatch, 1 = proximity, etc
                       as.integer(distrib),       # Poisson = 0 vs binomial = 1 (distribution of n)
                       as.integer(session.capthist),
                       as.double(unlist(session.xy)),  # polygon or transect or telemetry detection locations
                       as.double(session.signal),
                       as.integer(session.grp),
                       as.integer(nc),
                       as.integer(s),
                       as.integer(k), # may be zero-terminated vector for parts of polygon or transect detector
                       as.integer(m),
                       as.integer(ngroup),
                       as.integer(nmix),
                       as.integer(knownclass),     
                       as.double(trps),
                       as.double(distmat),         
                       as.double(usge),
                       as.integer(telemcode),   # 0 none, 1 IT, 2 DT, 3 CT 
                       as.integer(MRdata$markocc), 
                       as.integer(MRdata$Tu),
                       as.integer(MRdata$Tm),
                       as.integer(MRdata$Tn),
                       as.double(MRdata$chat),
                       as.double(unlist(session.mask)),
                       as.double(density),                    # density at each mask point x ngroup cols
                       ##                           as.double(pi.mask),                    # individual probability density if [1]>=0
                       as.double(MRdata$pi.mask),               # individual probability density if [1]>=0
                       as.double(Xrealparval),
                       as.double(Xrealparval0),
                       as.integer(nrow(Xrealparval)),            # number of rows in lookup table
                       as.integer(nrow(Xrealparval0)),           # ditto, naive
                       as.integer(PIA),                          # index nc,S,K,mix to rows Xrealparval
                       as.integer(PIA0),                         # index ngroup,S,K,mix to rows Xrealparval0
                       as.double(getcellsize(session.mask)),      # mask cell area or length
                       as.double(miscparm),                       # miscellaneous parameter
                       as.integer(detectfn),
                       as.integer(expandbinomN(details$binomN, dettype)),
                       as.double(details$minprob),
                       as.double(details$telemetryscale),
                       as.integer(details$debug),
                       components = double(6),                    # likelihood components 
                       value      = double(1),
                       resultcode = integer(1))
        }
        
        ####################################################
        if (details$debug>=1) {
            errorno <- match(temp$resultcode, c(0,1,5,6,8,9,11,51,53,54))
            cat("secrloglik resultcode ", temp$resultcode, " ", 
                switch(errorno, 
                       "successful completion",                      # 0
                       "generic failure",                            # 1     
                       "negative pimask[0] in filla0",               # 5 
                       "bad a0 in filla0",                           # 6 
                       "n exceeds N in binomial",                    # 8 
                       "non-finite value in secrloglik",             # 9 
                       "chat must be positive",                      # 11 
                       "negative D-Dmarked",                         # 51 
                       "zero sighting probability when T number >0", # 53 
                       "very negative Tlik in Tsightinglik"          # 54 
                ),"\n")   
            cat ('Likelihood components', round(temp$comp,4), '\n')
        }
        ####################################################
        
        LL <- if ((temp$resultcode != 0) | (temp$value < -1e9)) NA else LL + temp$value
        if (details$debug>=3 & temp$resultcode != 0) browser()
        if (temp$value == 0 & any(dettype==13)) {
            warning("zero likelihood with telemetry data suggests numerical problem",
                    " - try larger telemetryscale")
        }
        
        ####################################################
        ## unclear whether this is correct wrt groups
        if (logmult & all(detector(session.traps) %in% .localstuff$simpledetectors)) {
            LL <- LL + logmultinom(session.capthist,
                                   group.factor(session.capthist, groups))
        }
        LL
    } ## end sessionLL
    
   ###############################################################################################
    if (details$nsim > 0) {    ## overdispersion of sightings simulations only
        chat <- t(sapply (1:nsession, sessionLL))
        if (dim(chat)[2] == 3)  # otherwise C code has returned 1e10
            dimnames(chat) <- list(sessionlevels, c('Tu','Tm','Tn'))
        chat
    }
    else {
        if (ncores > 1) {
             clusterExport(clust, c("realparval", "details", "grp", "beta",
                                    "parindx"), envir = environment())
            loglik <- sum(parSapply(clust, 1:nsession, sessionLL))

        }
        else {
            loglik <- sum(sapply (1:nsession, sessionLL))
        }

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
        loglik <- ifelse(is.finite(loglik), loglik, -1e10)
        ifelse (neglik, -loglik, loglik)
    }
}
############################################################################################

