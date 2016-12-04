############################################################################################
## package 'secr'
## empirical.R
## design-based variance of derived density
## 2011 05 10 - derived.systematic quarantined in 'under development' folder
## 2012 03 17 - derived.nj extended to single n
## 2012 12 24 - clust argument required in derived.external
## 2014-10-17 revised to require spsurvey
## 2016-06-05 started derived.stratified
##            to do - treat sessions resulting from mash() clustergroups as strata
############################################################################################

## source ('d:\\density secr 2.5\\secr\\r\\empirical.R')

############################################################################

localvar <- function (z, xy) {
    ## z is response variable
    ## xy is dataframe or list with x, y
    ## choice of output verified by comparing to
    ## sum ((nj - n/J)^2) * J / (J-1) for ovensong vector
    ## with vartype='SRS'
    if (requireNamespace('spsurvey', quietly = TRUE))
        temp1 <- spsurvey::total.est(z = z, wgt = rep(1,length(z)),
                                     x = xy$x, y = xy$y)
    else
        stop ("package 'spsurvey' required for local variance")

    temp1[1,4]^2
}

derived.nj <- function ( nj, esa, se.esa, method = c('SRS','local','poisson','binomial'), xy = NULL,
    alpha = 0.05, loginterval = TRUE, area = NULL ) {
    method <- match.arg(method)
    esa <- unlist(esa)
    if (length(esa) == 2) {
        se.esa <- esa[2]
        esa <- esa[1]
    }
    n <- sum (nj)                          ## total animals
    J <- length(nj)                        ## number of replicates

    ## 2012-03-17
    if (J == 1) {
        method <- tolower(method)
        if (!(method %in% c('poisson','binomial'))) {
            method <- 'poisson'
            warning ("nj is of length 1; n variance component defaults to Poisson")
        }
        if (method == 'binomial') {
            if (is.null(area))
                stop ("must specify 'area' for binomial variance")
            N <- area * n / esa
        }
    }

    if ((method=='local') & is.null(xy))
        stop ("'local' method requires x and y coordinates")
    varn <- switch (method,
        local = localvar (nj, xy),         ## from total.est in spsurvey
        SRS = sum ((nj - n/J)^2) * J / (J-1),
        poisson = n,                                      ## 2012-03-17
        binomial = N * (esa / area * (1 - esa / area)),   ## 2012-03-17
        NA
    )
    D <- n / esa / J                       ## average density
    varD <- D^2 * (varn/n^2 + (se.esa/esa)^2)
    temp <- data.frame(row.names = 'D', estimate = D, SE.estimate = varD^0.5)
    temp <- add.cl(temp, alpha, loginterval)
    temp$CVn <- varn^0.5/n
    temp$CVa <- se.esa/esa
    temp$CVD <- varD^0.5/D
    attr(temp, 'nj') <- nj
    attr(temp, 'esa') <- as.numeric(esa)
    attr(temp, 'se.esa') <- as.numeric(se.esa)
    temp
}
############################################################################

derived.session <- function ( object, method = c('SRS','local'), xy = NULL,
    alpha = 0.05, loginterval = TRUE ) {
    method <- match.arg(method)
    if (!inherits (object, 'secr') | !ms(object))
        stop ("requires fitted multi-session model")
    if (('session' %in% object$vars) | ('Session' %in% object$vars) |
        (!is.null(object$sessioncov)))
        stop ("derived.session assumes detection model constant over sessions")
    der <-  derived(object, se.esa = TRUE, se.D = FALSE)[[1]][,1:2]
    nj <- sapply(object$capthist, nrow)    ## animals per session
    derived.nj(nj, der['esa','estimate'], der['esa','SE.estimate'], method,
               xy, alpha, loginterval)
}

############################################################################

derived.mash <- function (object, sessnum = NULL, method = c('SRS',
                          'local'), alpha = 0.05, loginterval = TRUE) {
    if (ms(object) & (is.null(sessnum))) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                output[[i]] <- derived.mash(object, sessnum = i, method, alpha, loginterval)
            }
            names(output) <- sessnames
            output
    }
    else {
        method <- match.arg(method)
        if (is.null(sessnum)) {
            n.mash <- attr(object$capthist, 'n.mash')
            centres <- attr(object$capthist, 'centres')
        }
        else {
            if (length(sessnum)>1)
                stop ("requires single sessnum")
            n.mash <- attr(object$capthist[[sessnum]], 'n.mash')
            centres <- attr(object$capthist[[sessnum]], 'centres')
        }

        if (is.null(n.mash))
            stop ("requires capthist with 'n.mash' attribute")

        if (ms(object))
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[[sessnum]][,1:2]
        else
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[,1:2]

        derived.nj(n.mash, der['esa','estimate'], der['esa','SE.estimate'],
                   method, centres, alpha, loginterval)

    }
}
####################################################################################
derived.cluster <- function (object, sessnum = NULL, method = c('SRS','local'), alpha = 0.05,
    loginterval = TRUE) {

    if (ms(object)) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                output[[i]] <- derived.cluster(object, sessnum = i, alpha, loginterval)
            }
            names(output) <- sessnames
            output
    }
    else {
        method <- match.arg(method)
        if (is.null(sessnum)) {
            capthist <- object$capthist
            sessnum <- 1
        }
        else
            capthist <- object$capthist[[sessnum]]

        if (!is.null(attr(capthist, 'n.mash')))
            warning ("use derived.mash for mashed data")

        nj <- cluster.counts(capthist)
        centres <- cluster.centres(traps(capthist))

        if (ms(object))
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[[sessnum]][,1:2]
        else
            der <-  derived(object, se.esa = TRUE, se.D = FALSE)[,1:2]

        derived.nj(nj, der['esa','estimate'], der['esa','SE.estimate'], method,
                   centres, alpha, loginterval)

    }
}
####################################################################################

derived.external <- function (object, sessnum = NULL, nj, cluster,
    buffer = 100, mask = NULL, noccasions = NULL, method = c('SRS','local'), xy = NULL,
    alpha = 0.05, loginterval = TRUE) {

    se.deriveesa <- function (selection, asess, noccasions) {
        CLesa  <- esagradient (object, selection, asess, noccasions, clust = NULL)
        sqrt(CLesa %*% object$beta.vcv %*% CLesa)
    }

    if (ms(object) & is.null(sessnum)) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                output[[i]] <- derived.external(object, sessnum = i, nj, cluster,
                    buffer, mask, noccasions, method, xy, alpha, loginterval)
            }
            names(output) <- sessnames
            output
    }
    else {
        method <- match.arg(method)
        if (!(inherits(cluster, 'traps')) | ms(cluster))
            stop ("cluster should be a single-session traps object")

        if (is.null(sessnum)) {
            capthist <- object$capthist
            sessnum <- 1
        }
        else
            capthist <- object$capthist[[sessnum]]

        if (is.null(noccasions))
            noccasions <- ncol(capthist)

        if (is.null(mask)) {
            mask <- make.mask(cluster, buffer = buffer)
        }

        ## substitute new traps and mask into old object
        traps(object$capthist) <- cluster
        object$mask <- mask

        ## use CLmeanesa for harmonic mean in case of individual variation
        esa.local <- CLmeanesa (object$fit$par, object, 1:nrow(capthist),
            sessnum, noccasions)
        se.esa <- se.deriveesa(1, sessnum, noccasions = noccasions)

        derived.nj (nj, esa.local, se.esa, method, xy, alpha, loginterval)
    }
}
####################################################################################


derived.stratified <- function (object, strata, mask = NULL, sessnum = NULL, plt = FALSE, ...) {

    if (ms(object)) {
        ## recursive call if MS
        stop("derived.stratified not yet implemented for multi-session fit")
        sessnames <- session(object$capthist)
        nsess <- length(sessnames)
        output <- vector('list', nsess )
        for (i in 1:nsess) {
            output[[i]] <- derived.stratified(object, strata, sessnum = i, mask, ...)
        }
        names(output) <- sessnames
        output
    }
    else {

        if (is.null(sessnum)) {
            capthist <- object$capthist
            sessnum <- 1
        }
        else
            capthist <- object$capthist[[sessnum]]

        if (is.null(mask))  ## multisession case?
            mask <- object$mask

        CH <- reduce(capthist, outputdetector = 'proximity')
        trapstrata <- group.factor(traps(CH), groups=strata)

        ## step 1: get stratum-specific counts of individuals nh
        nh <- tapply(apply(abs(CH), 3, sum), trapstrata, sum)
        if (nrow(CH) != sum(nh))
            warning ("some animals detected in more than one stratum")

        ## step 2: get stratum-specific a and D-hat
        ## WARNING: HAVE NOT ALLOWED
        ## INDIVIDUAL-SPECIFIC OR TRAP-SPECIFIC MODELS
        ## IS SUBSETTING PIA REQUIRED?
        stratalevels <- levels(trapstrata)
        nstrata <- length(stratalevels)
        if (is.null(mask))
             mask <- object$mask

        which.trap <- nearesttrap(mask, traps(CH))
        out <- vector('list', nstrata+1)
        for (i in 1:nstrata) {
            stratumi <- stratalevels[i]
            object$capthist <- subset(CH, traps = trapstrata == stratumi)
            if (all (strata %in% names(covariates(object$mask)))) {
                ## derive mask from covariates
                maskstrata <- group.factor(mask, groups = strata) ## not yet
            }
            else {
                maskstrata <- trapstrata[which.trap]
            }

            object$mask <- subset(mask, maskstrata == stratumi)
            tmp <- derived(object, se.esa = TRUE, se.D = TRUE, ...)
            out[[i]] <- data.frame(esa = tmp[1,1], SE.esa = tmp[1,2],
                                   wh = nrow(object$mask), tmp[2,])
        }
        out <- do.call(rbind, out)
        out$wh <- out$wh/sum(out$wh)
        out2 <- rbind(out, rep(NA,ncol(out)))
        rownames(out2) <- c(stratalevels,'Total')
        out2[nstrata+1, 'wh'] <- sum(out$wh)
        out2[nstrata+1, 'estimate'] <- sum(out$wh * out$estimate)
        varD <- sum(out$wh^2 * out$SE.estimate^2) ## IGNORING ANY COVARIANCES
        out2[nstrata+1, 'SE.estimate'] <- varD^0.5
        out2[nstrata+1, 'CVD'] <- out2[nstrata+1, 'SE.estimate'] / out2[nstrata+1, 'estimate']

        if (plt) {
            covariates(mask) <- data.frame(stratum = )
        }

        out2

    }
}
####################################################################################

