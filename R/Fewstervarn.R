############################################################################################
## package 'secr'
## fewstervarn.R
## adapted in part from RMF
## simplified from derivedSystematic 2018-12-24-27
############################################################################################

Fewstervarn <- function (nj, xy, design, esa, detectfn, detectpar, nocc,
                         basenx, df, maskspacing, keep, extrapolate = TRUE) {

    # design is a list with components --
    #   cluster -- traps object for one cluster)
    #   region  -- SpatialPolygons
    #   spacing -- distance between clusters)
    # and possibly --
    #   exclude 
    #   edgemethod
    #   exclmethod
    
    minxy <- bbox(design$region)[,1]
    
    #########################################################################
    ## Discrete sample space (bases for systematic samples)

    rx <- design$spacing             ## spacing between cluster centres
    drx <- rx / basenx               ## spacing between base points
    bx <- seq(drx/2, rx, drx)
    b <- expand.grid(x = bx, y = bx) ## assume square grid
    B <- basenx^2                    ## number of base points

    #########################################################################
    ## Boxlet centres clipped to region
    
    boxlets <- make.mask(type = 'polygon', poly = design$region, 
                         spacing = maskspacing)
    ## optionally clip out 'lakes'
    if (!is.null(design$exclude)) {
        boxlets <- subset(boxlets, !pointsInPolygon(boxlets, design$exclude))
    }

    #########################################################################
    ## Model trend
    
    njxy <- data.frame (nj = nj, xy)  ## dataframe with one row per cluster
    if (is.null(df)) {
        lgam <- gam (nj ~ s(x, y), data = njxy,
                     family = poisson(link = log), offset = log(esa))
    }
    else {
        lgam <- gam (nj ~ s(x, y, fx = TRUE, k = df+1), data = njxy,
                     family = poisson(link = log), offset = log(esa))
    }
    
    #########################################################################
    ## Sum over boxlets for each starting position in b
    ## bxy is 2-vector of x- and y- displacement corresp to b
    ## pdot (X, trps, ...) gives p.(x) = 1 - (1-p(x))^nocc for animal at X
    ## where p(x) is agregate over all detectors
    ## msk is subset of region near detectors

    sumoverboxlets <- function (bxy) {
        ## place grid for base 'bxy' and clip to region
        design$origin <- bxy + minxy
        trps <- do.call(make.systematic, design)
        ## proceed cluster by cluster for speed, and to reduce memory demand
        onecluster <- function (msk, trps) {
            pxy <- covariates(msk)$pxy
            a <- pdot (msk, trps, detectfn, detectpar, nocc)
            c(sum(pxy * a), sum(a) * attr(msk, 'area'))
        }
        splittrps <- split(trps, clusterID(trps))
        splitbox  <- split(boxlets, clusters = splittrps)
        bycluster <- mapply(onecluster, splitbox, splittrps, SIMPLIFY = TRUE)
        ## return 2-vector of sum(pxy.a), sum(a) over clusters
        apply(bycluster, 1, sum)   
    }

    #########################################################################
    ## predicted p, saved as covariate of boxlet mask

    pxy <- predict.gam(lgam, newdata = boxlets, type = "response")
    covariates(boxlets) <- data.frame(pxy = pxy)
    if (!extrapolate) {
        inner <- xy[chull(xy),]
        inner <- pointsInPolygon(boxlets, inner)
        innerboxlets <- subset(boxlets, inner)
        outerboxlets <- subset(boxlets, !inner)
        nearest <- nearesttrap(outerboxlets, innerboxlets)
        covariates(boxlets)$pxy[!inner] <- covariates(boxlets)$pxy[nearest]
    }
    covariates(boxlets)$pxy <- covariates(boxlets)$pxy / 
        sum(covariates(boxlets)$pxy, na.rm=TRUE)
    
    #########################################################################
    ## Qb including pdot weight

    detectfn <- valid.detectfn(detectfn)  ## converts character code to numeric 
    tmp <- apply(b, 1, sumoverboxlets)    ## 2 x b matrix
    Qb <- tmp[1,]
    Ab <- tmp[2,]

    #########################################################################
    ## Finally...

    A <- sum(esa)
    regionarea <- masksize(boxlets)       ## hectares or km
    N <- sum(nj) / A  * regionarea        ## H-T estimate
    out <- 1/B * sum ((N * Qb * (1-Qb) + N^2 * Qb^2) / Ab^2) -
        sum((N * Qb / Ab) / B)^2
    out <- out * A^2                      ## scale var(n/A) to var(n)
    
    #########################################################################
    ## And optionally save lots of intermediate values
    
    if (keep) {
        attr(out, 'xy') <- xy
        attr(out, 'N') <- N
        attr(out, 'design') <- design
        attr(out, 'b') <- b
        attr(out, 'boxlets') <- boxlets  # pxy is covariate of mask
        attr(out, 'Qb') <- Qb
        attr(out, 'Ab') <- Ab
    }
    out  # var (n)
}
############################################################################

derivedSystematic <- function ( object, xy, design = list(),
    basenx = 10, df = 9, extrapolate = TRUE, alpha = 0.05, loginterval = TRUE, 
    independent.esa = FALSE, keep = FALSE) {
    ## 'design' should comprise, e.g.,
    ## list(cluster = hollowsquare, spacing = 600, region = possumarea)
    
    warning("derivedSystematic is experimental in secr 3.2.0", call. = FALSE)
    
    if (!inherits (object, 'secr') | !is.null(object$groups))
        stop ("requires fitted secr model without groups")
    else if (ms(object)) {
        ## multisession, cluster = session
        if (('session' %in% object$vars) | ('Session' %in% object$vars) |
            (!is.null(object$sessioncov))) {
            stop ("derivedSystematic assumes detection model constant over sessions")
        }
        der <-  derived(object, se.esa = TRUE, se.D = FALSE, bycluster = FALSE)
        nj <- sapply(object$capthist, nrow)    ## animals per session
        maskspacing <- spacing(object$mask[[1]])
        detectpar <- detectpar(object)[[1]]
        nocc <- dim(object$capthist[[1]])[2]
        if (is.null(design$region)) {
            design$region <- attr(object$mask[[1]], 'polygon')
        }
    }
    else if (!is.null(clusterID(object$capthist))) {
        ## single session, native clusters
        der <-  derived(object, se.esa = TRUE, se.D = FALSE, bycluster = TRUE)
        nj <- cluster.counts(object$capthist)    ## animals per cluster
        maskspacing <- spacing(object$mask)
        detectpar <- detectpar(object)
        nocc <- dim(object$capthist)[2]
        if (is.null(design$region)) {
            design$region <- attr(object$mask, 'polygon')
        }
    }
    else {
        stop("derivedSystematic requires clustered input")
    }
    
    der <- do.call(rbind, der)
    esa <- der[grepl('esa', rownames(der)),1]
    se.esa <- der[grepl('esa', rownames(der)),2]
    
    J <- length(nj)
    
    if (is.null(design$region))
        stop ("region not found")
    
    ## convert matrix to SpatialPolygons if not already
    design$region <- boundarytoSP(design$region)      ## see utility.R
    design$exclude <- boundarytoSP(design$exclude)    ## may be NULL

    #########
    ## var(n)
    varn <- Fewstervarn (nj, xy, design, esa, object$detectfn, detectpar, nocc,
                         basenx, df, maskspacing, keep)
    ##########
    ## var (A)
    A <- sum(esa)
    if (independent.esa)
        varA <- sum(se.esa^2)
    else
        varA <- J^2 * sum(esa/A * se.esa^2)
     
    ##################
    ## assemble output
    n <- sum (nj)                          ## total animals
    D <- n / A 
    varD <- D^2 * (varn/n^2 + varA/A^2)
    temp <- data.frame(row.names = c('esa','D'), estimate = c(A,D), SE.estimate = c(varA,varD)^0.5)
    temp <- add.cl(temp, alpha, loginterval)
    temp$CVn <- c(NA, varn^0.5/n)
    temp$CVa <- c(NA, sqrt(varA)/A)
    temp$CVD <- c(NA, varD^0.5/D)
    attr(temp, 'nj') <- nj
    attr(temp, 'esa') <- esa
    attr(temp, 'se.esa') <- se.esa
    
    ######################
    ## optional attributes
    if (keep) {
        attr(temp, 'xy') <- xy
        attr(temp, 'design') <- attr(varn, 'design')
        attr(temp, 'N') <- attr(varn, 'N')
        attr(temp, 'b') <- attr(varn, 'b')
        attr(temp, 'boxlets') <- attr(varn, 'boxlets')  # pxy is covariate of mask
        attr(temp, 'Qb') <- attr(varn, 'Qb')
        attr(temp, 'Ab') <- attr(varn, 'Ab')
    }
        
    temp
}

plotSystematic <- function (out, dec = 0, legend = TRUE) {
    nj <- attr(out, 'nj')
    esa <- attr(out, 'esa')
    b <- attr(out, 'b')
    xy <- attr(out, 'xy')
    N <- attr(out, 'N')
    design <- attr(out, 'design')
    boxlets <- attr(out, 'boxlets')
    Qb <- attr(out, 'Qb')
    
    par(mfrow = c(2,3), mar = c(2,1,2,1), cex = 0.8)
    
    plot(design$region)
    text (xy[,1], xy[,2], as.character(nj))
    mtext(side = 3, line = -0.5, 'nj')
    
    plot(design$region)
    text (xy[,1], xy[,2], as.character(round(nj/esa,1)))
    mtext(side = 3, line = -0.5, 'nj/esa')
    
    plot(design$region)
    points(b[,1] + bbox(design$region)['x','min'], 
           b[,2] + bbox(design$region)['y','min'], xpd = TRUE)
    mtext(side = 3, line = -0.5, 'base points b')
    
    design$origin <- b[1,] + bbox(design$region)[,'min']
    trps <- do.call(make.systematic, design)
    plot(design$region)        
    plot(trps, detpar = list(cex=0.7), add = TRUE)
    mtext(side = 3, line = -0.5, 'grid for b[1]')
    
    covariates(boxlets)$pxy4 <- covariates(boxlets)$pxy4 <- covariates(boxlets)$pxy * 1e4
    plot(boxlets)
    plot(boxlets, cov = 'pxy4', breaks = 5, add = TRUE, legend = legend)
    mtext(side = 3, line = -0.5, 'pxy * 1E4')
    plot(design$region, add = TRUE)

    opar <- par(pty = 's', mar = c(4,4,2,1))
    plot(b, type = 'n', xlab = '', ylab = '', axes = FALSE, 
         xlim=c(0,design$spacing), ylim=c(0,design$spacing))
    axis(1, at = unique(b$x))
    axis(2, at = unique(b$y))
    text (b[,1], b[,2], round(N*Qb, dec), cex = 0.7)
    mtext(side = 3, line = 1, 'E(n|b)')
    par(opar)
        
}

# plotSystematic(tmp)

