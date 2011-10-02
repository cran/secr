############################################################################################
## package 'secr'
## regionN.R
## population size in an arbitrary region
## 2011-08-18 (fixed expected.n)
## 2011-09-26 fixed multi-session bug in region.N
############################################################################################

predictD <- function (object, region, session) {
    ## not exported
    ## groups not thought through yet...
    ## if (!is.null(group))
    ##    newdata <- subset (newdata[newdata$g == group,]

    if (is.null(object$model$D)) {    ## object$CL == TRUE
        temp <- derived(object) ## inefficient as repeats for each sess
        if (!is.data.frame(temp))
            temp <- temp[[session]]
        D <- temp['D', 'estimate']
        return (D)
    }
    else {
        newdata <- D.designdata (region, object$model$D,
             group.levels(object$capthist,object$groups),
             session(object$capthist), sessioncov = object$sessioncov)
        ## select a single session
        if ("session" %in% names(newdata)) {
# following does not work 2011-09-26
#            newdata <- subset(newdata, subset = (newdata$session == session))
            newdata <- newdata[newdata$session == session,]
        }
        if ("Session" %in% names(newdata)) {
#            newdata <- subset(newdata, subset = newdata$Session ==
#               as.numeric(factor(session, levels =  session(object$capthist)))-1)
            newdata <- newdata[newdata$Session ==
               (as.numeric(factor(session, levels =  session(object$capthist)))-1),]
        }
        class(newdata) <- c('mask', 'data.frame')
        attr (newdata, 'area') <- attr(region, 'area')
        indx <- object$parindx$D
        beta <- object$fit$par[indx]
        if (object$model$D == ~1) {
            D <- untransform(beta, object$link$D)
            return ( rep(D, nrow(region) ))
        }
        else {
            vars <- all.vars(object$model$D)
            if (any(!(vars %in% names(newdata))))
                stop ("one or more model covariates not found")
            newdata <- as.data.frame(newdata)
            mat <- model.matrix(object$model$D, data=newdata)
            lpred <- mat %*% beta
            return ( untransform(lpred, object$link$D) )
        }
    }
}

############################################################################################

region.N <- function (object, region = NULL, spacing = NULL, session = NULL,
    group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,
    keep.region = FALSE, nlowerbound = TRUE) {

    ## Notes
    ## se.N = FALSE returns scalar N

    ###########################################################
    ## for gradient of E.N wrt density betas
    betaN <- function (betaD, object, region, session, group) {
        ## assume region is a mask (no need for spacing)
        ## assume single session
        ## indx identifies beta parameters for density D
        object$fit$par[indx] <- betaD
        region.N(object, region, spacing = NULL, session = session,
            group = group, se.N = FALSE, keep.region = FALSE)
    }
    ###########################################################
    ## for gradient of R.N wrt all betas
    betaRN <- function (beta, object, region) {
        ## assume region is a mask (no need for spacing)
        ## assume single session
        ## n, cellarea, sessnum global
        object$fit$par <- beta
        D <- predictD(object, region, session)
        n + sumDpdot (object, sessnum, region, D, cellarea,
                  constant = FALSE, oneminus = TRUE)[1]
    }
    ###########################################################

    if (is.null(region)) {
        region <- object$mask
        ## warning ("using entire mask as region")
    }

    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object ")

    if (!is.null(group))
        stop ("not yet working for groups")

    if (!all(group %in% interaction (object$groups)))
        stop ("unrecognised groups")

    if (is.null(session))
        session <- session(object$capthist)

    ####################################################################
    ## if N requested for multiple sessions,
    ## call region.N recursively for each session
    nsess <- length(session)
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in session) {
            if (ms(region))
                tempregion <- region[[sess]]
            else
                tempregion <- region
            out[[sess]] <- region.N (object, region = tempregion,
                spacing = spacing, session = sess, group = group,
                se.N = se.N, alpha = alpha, loginterval = loginterval,
                keep.region = keep.region)
        }
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {

        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask
        ########################################################
        ## if necessary, convert vector region to raster
        if (!inherits(region, 'mask')) {
            if (is.null(spacing)) {
                ## use mask spacing by default
                spacing <- spacing(mask)
            }
            if (inherits(region, 'SpatialPolygonsDataFrame')) {
                if (!require (sp))
                    stop ("package 'sp' required for SpatialPolygonsDataFrame ",
                          "in region.N")
                bbox <- bbox(region)
            }
            else {
                bbox <- apply(region, 2, range)
            }
            region <- make.mask(bbox, poly = region, buffer = 0,
                spacing = spacing, type = 'polygon',
                check.poly = FALSE)
        }
        #######################################################

        ## inserted to fix bug
        ## 2011-09-26
        if (ms(region))
            region <- region[[session]]

        ## region now inherits from mask, so has area attribute
        cellarea <- attr(region, 'area')
        regionarea <- nrow(region) * cellarea
        if (ms(object))
            n <- nrow(object$capthist[[session]])
        else
            n <- nrow(object$capthist)
        sessnum <- match (session, session(object$capthist))

        #######################################################
        ## for conditional likelihood fit,
        if (object$CL) {
            temp <- derived(object) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
            seD <- temp['D', 'SE.estimate']
            N <- D * regionarea
            if (!se.N) return (N)    ## and stop here
            seN <- seD * regionarea
        }

        #######################################################
        ## for full likelihood fit...
        else {
            if (is.null(object$model$D) | is.null(object$link$D))
                stop ("model or link function not found in object")

            if (object$model$D == ~1) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- predicted['D','estimate']
                N <- D * regionarea
                seN <- predicted['D','SE.estimate'] * regionarea
            }
            else {
                D <- predictD (object, region, session)
                N <- sum(D) * cellarea
                if (!se.N) return (N)    ## and stop here
                indx <- object$parindx$D
                dNdphi <- gradient (object$fit$par[indx],
                    betaN, object = object, region = region, session =
                    session, group = group)
                beta.vcv <- object$beta.vcv[indx,indx]
                seN <- (dNdphi %*% beta.vcv %*% dNdphi)^0.5
            }
        }

        #######################################################################
        ## realised N
        ## only makes sense for individual detectors (not unmarked or presence)
        ## amended 2011-09-26
        if (ms(object))
            det <- detector(traps(object$capthist)[[session]])
        else
            det <- detector(traps(object$capthist))
        if (det %in% .localstuff$individualdetectors) {
            notdetected <- sumDpdot (object, sessnum, region, D,
                cellarea, constant = FALSE, oneminus = TRUE)[1]
            dNdbeta <- gradient (object$fit$par, betaRN, object =
                object, region = region)
            RN <- n + notdetected
            pdotvar <- dNdbeta %*% object$beta.vcv %*% dNdbeta
            seRN <- (notdetected + pdotvar)^0.5
            En <- sumDpdot (object, sessnum, mask, D, attr(mask,'area'),
                 constant = FALSE, oneminus = FALSE)[1]
        }
        else { RN <- NA; seRN <- NA; En <- NA }
        #######################################################################

        temp <- data.frame(
            row.names = c('E.N','R.N'),
            estimate = c(N,RN),
            SE.estimate = c(seN,seRN))
        ## lower bound added 2011-07-15
        if (nlowerbound)
            temp <- add.cl (temp, alpha, loginterval, c(0, n))
        else
            temp <- add.cl (temp, alpha, loginterval, c(0, 0))
        temp$n <- rep(n, nrow(temp))
        temp$E.n <- rep(round(En,2), nrow(temp))
        if (keep.region)
            attr(temp, 'region') <- region
        temp
    }
}

############################################################################################
## 2011-05-05
############################################################################################

## see also Dpdot in pdot.R
## modelled on esa.R

sumDpdot <- function (object, sessnum = 1, mask, D, cellarea, constant = TRUE,
                      oneminus = FALSE)

# Return integral for given model and new mask, D
# 'sessnum' is integer index of session (factor level of the 'session' attribute in capthist)
# object must have at least capthist, detectfn
# D should be scalar or vector of length nrow(mask)

# if 'constant' a much simplified calculation is used, assuming
# constant detection and density, and full detector usage

{
    if (ms(object))
        capthists <- object$capthist[[sessnum]]
    else
        capthists <- object$capthist

    if (ms(mask))
        mask <- mask[[sessnum]]

    beta <- object$fit$par

    traps   <- attr(capthists, 'traps')  ## need session-specific traps
    ## 2011-09-26 check added
    if (!(detector(traps) %in% .localstuff$individualdetectors))
        stop ("require individual detector type for sumDpdot")
    dettype <- detectorcode(traps)
    n       <- max(nrow(capthists), 1)
    s       <- ncol(capthists)

    noccasions <- s

    nmix    <- object$details$nmix
    nmix    <- ifelse (is.null(nmix), 1, nmix)

    ##############################################
    ## adapt for marking occasions only 2009 10 24
    ## hangover from esa... leave for now
    q <- attr(capthists, 'q')
    if (!is.null(q))
        if (q<s) s <- q
    ##############################################

    if (dettype %in% c(3,6)) {
        k <- c(table(polyID(traps)),0)
        K <- length(k)-1
    }
    else if (dettype %in% c(4,7)) {
        k <- c(table(transectID(traps)),0)
        K <- length(k)-1
    }
    else {
        k <- nrow(traps)
        K <- k
    }
    m      <- length(mask$x)            ## need session-specific mask...
    if (constant) {
        if (is.null(beta))
            real <- detectpar(object)
        else {
            real <- makerealparameters (object$design0, beta,
                object$parindx, object$link, object$fixed)  # naive
            real <- as.list(real)
            names(real) <- parnames(object$detectfn)
        }
        a <- cellarea * sum(pdot(X = mask, traps = traps, detectfn = object$detectfn,
                             detectpar = real, noccasions = noccasions))
        return(a * D)
    }
    else {

        ## allow for old design object
        if (length(dim(object$design0$PIA))==4)
            dim(object$design0$PIA) <- c(dim(object$design0$PIA),1)
        PIA <- object$design0$PIA[sessnum,,1:s,,,drop=F]
        ncolPIA <- dim(object$design0$PIA)[2]
        #############################################
        ## trick to allow for changed data 2009 11 20
        ## nmix>1 needs further testing 2010 02 26
        ## NOTE 2010-11-26 THIS LOOKS WEAK
        if (dim(PIA)[2] != n) {
            PIA <- array(rep(PIA[1,1,,,],n), dim=c(s,K,nmix,n))
            PIA <- aperm(PIA, c(4,1,2,3))   ## n,s,K,nmix
            ncolPIA <- n     ## 2010 02 26
        }
        #############################################

        realparval0 <- makerealparameters (object$design0, beta,
            object$parindx, object$link, object$fixed)  # naive
        Xrealparval0 <- scaled.detection (realparval0, FALSE,
            object$details$scaleg0, NA)

        used <- usage(traps)
        if (any(used==0))
        PIA <- PIA * rep(rep(t(used),rep(n,s*K)),nmix)
        ncolPIA <- n

        param <- object$details$param
        if (is.null(param))
            param <- 0    ## default Borchers & Efford (vs Gardner & Royle)
        gamma <- 1  ## DUMMY

        ## add density as third column of mask
        if (!(length(D) %in% c(1,nrow(mask))))
            stop ("D does not match mask in sumDpdot")
        if (length(D) == 1)
            D <- rep(D[1], nrow(mask))
        mask <- cbind (mask, D)

        useD <- TRUE
        temp <- .C("integralprw1", PACKAGE = 'secr',
            as.integer(dettype),
            as.integer(param),
            as.double(Xrealparval0),
            as.integer(n),
            as.integer(s),
            as.integer(k),
            as.integer(m),
            as.integer(nmix),
            as.double(unlist(traps)),
            as.double(as.numeric(unlist(mask))),
            as.integer(nrow(Xrealparval0)), # rows in lookup
            as.integer(PIA),                # index of nc*,S,K to rows in realparval0
            as.integer(ncolPIA),            # ncol - if CL, ncolPIA = n, else ncolPIA = 1 or ngrp
            as.double(cellarea),
            as.double(gamma),
            as.integer(object$detectfn),
            as.integer(object$details$binomN),
            as.double(object$details$cutval),
            as.integer(useD),
            a=double(n),
            resultcode=integer(1)
       )
       if (oneminus) {
           temp$a <- sum(D) * cellarea - temp$a
       }
       if (temp$resultcode == 3)
           stop ("groups not implemented in external function 'integralprw1'")
       if (temp$resultcode != 0)
           stop ("error in external function 'integralprw1'")
       return(temp$a)
    }
}
############################################################################################


expected.n <- function (object, session = NULL, group = NULL, bycluster = FALSE,
                        splitmask = FALSE) {

    ## Note
    ## splitmask toggles between two methods for clustered detectors:
    ## 1. is to integrate over the whole mask, restricting detection to each cluster in turn
    ## 2. is to split mask into Dirichlet subregions by distance to detector centre
    ## and to integrate over all detectors, assuming those far away will never detect from
    ## within a subregion to which they do not belong
    ## Probably, 1. is more robust

    if (!all(session %in% session(object$capthist)))
        stop ("session incompatible with object")

    if (!is.null(group))
        stop ("not yet working for groups")

    if (!all(group %in% interaction (object$groups)))
        stop ("unrecognised groups")

    if (is.null(session))
        session <- session(object$capthist)

    ####################################################################
    ## if En requested for multiple sessions,
    ## call En recursively for each session
    nsess <- length(session)
    if (nsess > 1) {
        ## predict for each session
        out <- vector('list')
        for (sess in session) {
            out[[sess]] <- expected.n (object, sess, group, bycluster)
        }
        out
    }

    ####################################################################
    ## otherwise, this is a non-recursive call for one session...
    else {

        if (ms(object$mask))
            mask <- object$mask[[session]]
        else
            mask <- object$mask

        cellarea <- attr(mask, 'area')
        if (ms(object)) {
            n <- nrow(object$capthist[[session]])
            trps <- traps(object$capthist[[session]])
        }
        else {
            n <- nrow(object$capthist)
            trps <- traps(object$capthist)
        }
        sessnum <- match (session, session(object$capthist))

        #######################################################
        ## for conditional likelihood fit,
        if (object$CL) {
            temp <- derived(object) ## inefficient as repeats for each sess
            if (!is.data.frame(temp))
                temp <- temp[[session]]
            D <- temp['D', 'estimate']
        }

        #######################################################
        ## for full likelihood fit...
        else {
            if (is.null(object$model$D) | is.null(object$link$D))
                stop ("model or link function not found in object")

            if (object$model$D == ~1) {
                predicted <- predict(object)
                if (!is.data.frame(predicted))
                    predicted <- predicted[[1]]
                D <- predicted['D','estimate']
            }
            else {
                D <- predictD (object, mask, session)
            }
        }
        #################################################################
        if (length(D) == 1)
            D <- rep(D, nrow(mask))

        #################################################################
        if (bycluster) {
            centres <- cluster.centres(trps)
            nclust <- nrow(centres)
            out <- numeric (nclust)
            if (is.null(attr(trps, 'cluster'))) {
                clusterID(trps) <- 1:nclust
            }
            if (splitmask) {
                cluster <- nearesttrap (mask, centres)
                mask <- split (mask, cluster)
                D <- split(D, cluster)
            }
            for (i in 1:nclust) {
                if (splitmask) {
                    out[i] <- sumDpdot(object = object, sessnum = sessnum,
                        mask=mask[[i]], D = D[[i]], cellarea = cellarea,
                        constant = FALSE, oneminus = FALSE)[1]
                }
                else {
## replaced 2011-08-18
##                    traps(object$capthist) <- subset(trps, subset =
##                        as.numeric(clusterID(trps)) == i)

                    temptrap <- subset(trps, subset = as.numeric(clusterID(trps)) == i)
                    if (ms(object))
                        traps(object$capthist[[sessnum]]) <- temptrap
                    else
                        traps(object$capthist) <- temptrap

                    out[i] <- sumDpdot(object = object, sessnum = sessnum,
                        mask=mask, D = D, cellarea = cellarea,
                        constant = FALSE, oneminus = FALSE)[1]
                }
            }
            out
        }
        else sumDpdot (object, sessnum, mask, D, attr(mask,'area'),
             constant = FALSE, oneminus = FALSE)[1]
        #################################################################
    }
}

