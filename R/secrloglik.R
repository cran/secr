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

secr.lpredictor <- function (model, newdata, indx, beta, field, beta.vcv=NULL) {
    vars <- all.vars(model)
    if (any(!(vars %in% names(newdata))))
        stop ("one or more model covariates not found in 'newdata'")
    newdata <- as.data.frame(newdata)
    lpred <- matrix(ncol = 2, nrow = nrow(newdata),dimnames=list(NULL,c('estimate','se')))
    mat <- model.matrix(model, data = newdata)

    ## drop pmix beta0 column from design matrix (always zero)
    if (field=='pmix') {
        mat <- mat[,-1,drop=FALSE]
    }
    lpred[,1] <- mat %*% beta[indx]

    ## replaced for robustness 2013-06-06
    ## if (is.null(beta.vcv)) return ( cbind(newdata,lpred) )
    if (is.null(beta.vcv) | (any(is.na(beta[indx])))) return ( cbind(newdata,lpred) )
    else {
        vcv <- beta.vcv[indx,indx, drop = FALSE]    ## OR maybe all betas?
## bug fix 2013-04-14
## nrw <- 1:nrow(mat)
## vcv <- apply(expand.grid(nrw, nrw), 1, function(ij)
        nrw <- nrow(mat)
        vcv <- apply(expand.grid(1:nrw, 1:nrw), 1, function(ij)
            mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F]))  # link scale
        vcv <- matrix (vcv, nrow = nrw)
        lpred[,2] <- diag(vcv)^0.5
        temp <- cbind(newdata,lpred)
        attr(temp, 'vcv') <- vcv
        return(temp)
    }
}

############################################################################################

makerealparameters <- function (design, beta, parindx, link, fixed) {
    modelfn <- function(i) {
        ## linear predictor for real parameter i
        Yp <- design$designMatrices[[i]] %*% beta[parindx[[i]]]
        if (names(link)[i] == 'pmix') {
            ## 2013-04-14 index of class groups (pmix sum to 1.0 within latentmodel)
            h2 <- grep('.h2',dimnames(design$designMatrices[[i]])[[2]], fixed=T)
            h3 <- grep('.h3',dimnames(design$designMatrices[[i]])[[2]], fixed=T)
            tmp <- design$designMatrices[[i]][,-c(h2,h3), drop = FALSE]
            tmph <- design$designMatrices[[i]][,c(h2,h3), drop = FALSE]
            latentmodel <- as.numeric(factor(apply(tmp,1,paste, collapse='')))
            refclass <- apply(tmph,1,sum) == 0
            Yp[refclass] <- NA
            Yp <- mlogit.untransform(Yp, latentmodel)
            Yp[design$parameterTable[,i]]
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
    OK <- names(fixed) %in% detectionparameters
    fixed.dp <- fixed[OK]
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

getD <- function (designD, beta, mask, parindx, link, fixed,
                  grouplevels, sessionlevels) {
    if (ms(mask))
        nmask <- max(sapply(mask, nrow))
    else
        nmask <- nrow(mask)
    ngroup <- length(grouplevels)
    nsession <- length(sessionlevels)
    D <- array(dim = c(nmask, ngroup, nsession))
    dimnames(D) <- list(1:nrow(D), grouplevels, sessionlevels)
    if (!is.null(fixed$D)) {
        D[,,] <- fixed$D
    }
    else {
        if (is.function(designD))
            D[,,] <- designD(beta[parindx$D], mask, ngroup, nsession)
        else {
            D[,,] <- designD %*% beta[parindx$D]   # linear predictor
            D[,,] <- untransform (D, link$D)
        }
        D[D<0] <- 0                            # silently truncate at zero
    }
    D
}
###############################################################################

secr.loglikfn <- function (beta, parindx, link, fixedpar, designD, design,
    design0, capthist, mask, detectfn, CL, hcov, groups, details, logmult, ncores,
    clust, dig = 3, betaw = 10, neglik = TRUE)

# Return the negative log likelihood for inhomogeneous Poisson spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function
#    0 = halfnormal, 1 = hazard, 2 = exponential etc.
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
# details$trace=T sends a one-line report to the screen

{
    ## for models fitted with old versions we need to fillin these values
    if (is.null(details$minprob)) details$minprob <- 1e-50
    if (is.null(details$debug)) details$debug <- FALSE   ## 2012-10-28
    if (is.null(details$ignoreusage)) details$ignoreusage <- FALSE  ## 2013-01-23
    if (is.null(details$unmash)) details$unmash <- FALSE ## 2013-01-23

    if (ms(capthist))
        sessionlevels <- session(capthist)
    else
        sessionlevels <- 1
    nsession <- length(sessionlevels)
    if ((ncores>1) & missing(clust))
        stop("not ready for multicore here")

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
    detparindx <- parindx[!(names(parindx) %in% c('D'))]
    detlink <- link[!(names(link) %in% c('D'))]
    realparval  <- makerealparameters (design, beta, detparindx, detlink, fixedpar)
    realparval0 <- makerealparameters (design0, beta, detparindx, detlink, fixedpar)

    #--------------------------------------------------------------------
    # Density
    D.modelled <- !CL & is.null(fixedpar$D)
    if (!CL ) {
        D <- getD (designD, beta, mask, parindx, link, fixedpar,
                    levels(grp[[1]]), sessionlevels)
        if (sum(D) <= 0)
            warning ("invalid density <= 0")
    }
    #--------------------------------------------------------------------

   ###############################################################################################
   ###############################################################################################
   sessionLL <- function (sessnum) {

        ## in multi-session case must get session-specific data from lists
        if (ms(capthist)) {
            session.capthist <- capthist[[sessnum]]
            session.traps    <- traps(capthist)[[sessnum]]
            session.mask     <- mask[[sessnum]]
            session.grp      <- grp[[sessnum]]
            session.xy       <- 0
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
        sessg <- min (sessnum, design$R)
        cell <- attr(session.mask,'area')
        session.mask <- as.matrix(session.mask[,1:2])
        knownclass <- getknownclass(session.capthist, nmix, hcov)

        LL <- 0

        #----------------------------------------
        # 2012-10-20,21
        # mingled proximity & telemetry
        xylist <- attr(session.capthist,'xylist')
        if (!is.null(xylist)) {
            if (.localstuff$iter < 1) {
                outside <- outsidemask (session.capthist, session.mask)
                if (sum(outside)>0)
                    warning (sum(outside), " centres lie outside mask and",
                             " will be assigned to the nearest mask point")
            }
            telem <- row.names(session.capthist) %in% names(xylist)
            T.session.capthist <- subset(session.capthist, telem)
            session.capthist <- subset(session.capthist, !telem)
            nc <- nrow(session.capthist)
            ## prwi for 'known' centres
            LL <- LL + telemloglik(T.session.capthist,
                session.traps, session.mask, detectfn=detectfn, detectpar=realparval )

            ## realparval is temporary; does not allow for non-constant model or scaling!

            ## optional use of telemetry locations to inform sigma
            if (details$telemetrysigma) {
                index <- ifelse(nrow(realparval)==1, 1, sessnum)
                sigma <- realparval[index,'sigma']
                z <- ifelse (detectfn %in% c(1,3,5,6,7,8), realparval[index,'z'], 1)
                LL <- LL + telemetryloglik (T.session.capthist, detectfn, sigma, z)$value
            }
        }

        #----------------------------------------

        dettype <- detectorcode(session.traps)

        if (dettype %in% c(5,9,12)) {    # cue or signal strength
            session.signal <- signal(session.capthist)
            if (dettype==12)
                session.signal <- c (session.signal, noise(session.capthist))
            session.signal <- switch( details$tx,
                log = log(session.signal),
                logit = logit(session.signal),
                identity = session.signal
            )
            ## 2012-01-30 code missing values as negative for C code
            session.signal[is.na(session.signal)] <- -1
        }
        else
            session.signal <- 0

        ## miscparm is used to package beta parameters that are not modelled
        ## and hence do not have a beta index specified by parindx.
        ## This includes the signal threshold and the mean and sd of noise.

        miscparm <- c(0,0,0)
        if (dettype == 9)   ## cue rate
            miscparm[] <- exp(beta[max(unlist(parindx))+1])   ## fudge: last
        else if (detectfn %in% c(12,13))   ## experimental signal-noise
            miscparm[] <- c(details$cutval, beta[max(unlist(parindx))+1:2])   ## fudge: last 2
        else if (detectfn %in% c(10,11))  ## Dawson&Efford 2009 models
            miscparm[] <- details$cutval
        if (dettype %in% c(3,6,13)) {    # polygonX, polygon, telemetry
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

#        if (dettype == 8) {    # times -- phony use of 'signal'
#            session.signal <- times(session.capthist)
#        }

        trps <- unlist(session.traps, use.names=F)
        usge <- usage(session.traps)
        if (is.null(usge) | details$ignoreusage)
            usge <- matrix(1, nrow = K, ncol = s)

        #---------------------------------------------------
        binomN <- details$binomN
        #---------------------------------------------------

        ## differentiate so density & g do not both need to use sessions
        if (CL)
            density <- 0
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

        #------------------------------------------
        # allow for scaling of detection

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
        # For conditional likelihood, supply a value for each
        #  animal, not just groups

        K <- ifelse (detector(session.traps) %in% .localstuff$polydetectors, length(k)-1, k)
        if (CL) tempPIA0 <- design0$PIA[sessg,1:nc,1:s,1:K, ]
        else tempPIA0 <- design0$PIA[sessg,1:ngroup,1:s,1:K, ]     ## drop=FALSE unnecessary?

        if (is.numeric(details$distribution)) {
            if (details$distribution < nc)
                stop ("superpopulation (details$distribution) ",
                      "is less than number observed")
            distrib <- details$distribution
        }
        else
            distrib <- switch (tolower(details$distribution), poisson=0, binomial=1, 0)

        #-------------------------------------------
        # debug
        #-------------------------------------------
        if (details$debug) {

           arglist <- list(CL=CL, dettype=dettype,
                           details=details, distrib=distrib,
                           capthist=session.capthist,
                           xy=unlist(session.xy),
                           signal=session.signal, grp=session.grp,
                           nc=nc, s=s, k=k, m=m, ngroup=ngroup,
                           nmix=nmix, trps=trps, usge=usge,
                           mask=session.mask, density=density,
                           Xrealparval=Xrealparval,
                           Xrealparval0=Xrealparval0,
                           nrowXrealparval=nrow(Xrealparval),
                           nrowXrealparval0=nrow(Xrealparval0),
                           PIA=design$PIA[sessg,1:nc,1:s,1:K,],
                           tempPIA0=tempPIA0, cell=cell,
                           miscparm=miscparm, detectfn=detectfn,
                           binomN=binomN)

           save(arglist, file= paste('arglist',sessnum,'.RData',sep=''))
       }

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
            ## simple      (Royle-Nichols)
            ## integrated
            ## pairwise

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
            if ((.localstuff$iter < 1) & (details$presence %in% c('simple','pairwise')) &
                !(detectfn %in% c(1,4,7,8)))
                warning ("simple presence requires detectfn for which",
                    " sigma = radius (1,4,7,8)", .call = FALSE)
            type <- switch (details$presence, simple = 0, integrated = 1, pairwise = 2, 3)

            temp <- .C('presenceloglik', PACKAGE = 'secr',
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
        #-------------------------------------------
        # experimental 'telemetry' detector type
        #-------------------------------------------
        else if (dettype == 13) {
            ## xy locations of individuals
            ## does not allow within-session models
            index <- ifelse(nrow(Xrealparval)==1, 1, sessnum)
            sigma <- Xrealparval[index,'sigma']
            z <- ifelse (detectfn %in% c(1,3,5,6,7,8), Xrealparval[index,'z'], 1)
            temp <- telemetryloglik (session.capthist, detectfn, sigma, z)
        }
        #--------------------------------------------
        # typical call (not 'presence' or 'unmarked')
        #--------------------------------------------
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
                as.integer(k), # may be zero-terminated vector for parts of polygon
                               # or transect detector
                as.integer(m),
                as.integer(ngroup),
                as.integer(nmix),
                as.integer(knownclass), ## 2013-04-12
                as.double(trps),
                as.double(usge),
                as.double(session.mask),
                as.double(density),                          # density at each mask point x ngroup cols
                as.double(Xrealparval),
                as.double(Xrealparval0),
                as.integer(nrow(Xrealparval)),               # number of rows in lookup table
                as.integer(nrow(Xrealparval0)),              # ditto, naive
                as.integer(design$PIA[sessg,1:nc,1:s,1:K,]), # index of nc,S,K,mix to rows in Xrealparval
                as.integer(tempPIA0),                        # index of ngroup,S,K,mix to rows in Xrealparval0
                as.double(cell),                             # mask cell area
                as.double(miscparm),                         # miscellaneous parameter
                as.integer(detectfn),
                as.integer(binomN),
                as.double(details$minprob),
                a=double(nc),
                value=double(1),
                resultcode=integer(1))
        }

        LL <- ifelse (temp$resultcode == 0, LL + temp$value, NA)

        ####################################################
        ## unclear whether this is correct wrt groups
        if (logmult & (detector(session.traps) %in% .localstuff$simpledetectors)) {
            LL <- LL + logmultinom(session.capthist,
                                   group.factor(session.capthist, groups))
        }
        LL

    } ## end sessionLL
   ###############################################################################################

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
############################################################################################

