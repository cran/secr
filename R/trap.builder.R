###############################################################################
## package 'secr'
## trap.builder.R
## repeat trap layout systematically, by GRTS, or at simple random centres
## across a region
## Also make.systematic, mash(), cluster.counts(), cluster.centres()
## 2011-08-16 (full argument names)
## 2012-01-11 (cutval and signal in mash(); mxy[1:2])
## 2012-02-07 mash noise
###############################################################################

## spsurvey uses sp

boundarytoSPDF <- function (boundary) {
    ## build sp SpatialPolygonsDataFrame object
    ## input is 2-column matrix for a single polygon
    ## requires package sp
    Sr1 <- Polygon(boundary)
    Srs1 <- Polygons(list(Sr1), "s1")
    SpP <- SpatialPolygons(list(Srs1))
    attr <- data.frame(a = 1, row.names = "s1")
    SpatialPolygonsDataFrame(SpP, attr)
}

###############################################################################

trap.builder <- function (n = 10, cluster, region = NULL, frame =
    NULL, method = 'SRS', edgemethod = 'clip', samplefactor = 2,
    ranks = NULL, rotation = NULL, detector, plt = FALSE, add = FALSE) {

    ## region may be -
    ## matrix x,y
    ## sp SpatialPolygonsDataFrame object

    ## future: allow
    ##   minimum separation
    ##   shapefile

    #####################################################
    # 1. form polygon
    # 2. get n random origins SRS, GRTS
    # 3. translate n times; covariate for cluster of origin
    # 4. optionally clip reject edge clusters
    # 5. rbind
    #####################################################

    if (!require(sp))
        stop ("package 'sp' required in trap.builder")
    .local <- new.env()   ## for clusteri

    allinside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy))
        !any(is.na(overlay (xy, region)))
    }

    position <- function (i, cluster) {
        newtraps <- shift(cluster, origins[i,])
        if (!is.null(rotation)) {
            if (rotation<0)
                rotation <- runif(1) * 360
            newtraps <- rotate(newtraps, rotation, apply(newtraps,2,mean))
        }
        i <- .local$clusteri
        .local$clusteri <- .local$clusteri + 1    ## global update!
        clusterID(newtraps) <- factor(rep(i,nrow(newtraps)), levels=i)
        clustertrap(newtraps) <- as.numeric(polyID(newtraps))
        newtraps
    }

    ## option for single-trap clusters
    if (missing(cluster))
        cluster <- NULL
    if (is.null(cluster)) {
        if (missing(detector))
            detector <- 'multi'
        if (!(detector %in% .localstuff$pointdetectors))
            stop ("solitary detectors must be of a point detector type")
        cluster <- make.grid(nx = 1, ny = 1, detector = detector)
        edgemethod <- 'allowoverlap'
    }
    else {
        if ((attr(cluster,'detector') %in% .localstuff$polydetectors) &
            (ndetector(cluster) > 1))
            stop("clusters with multiple polygons or transects not supported")
    }

    if (method == 'all') {
        n <- nrow(frame)
    }
    if (is.null(region))
        edgemethod <- 'allowoverlap'

    if (is.null(frame)) {
        if (is.null(region)) {
            stop ("specify at least one of 'region' or 'frame'")
        }
        SPDF <- inherits(region, 'SpatialPolygonsDataFrame')
        if (!SPDF) {
            region <- matrix(unlist(region), ncol = 2)
            region <- rbind (region, region[1,])  # force closure of polygon
            region <- boundarytoSPDF(region)
        }
        if (plt & !add)
            plot(region)
    }
    else {
        if (plt & !add) {
            if (!is.null(region))
                plot(region)
            else {
                require (MASS)
                eqscplot (frame, axes=F, xlab='', ylab='', pch=1, cex=0.5)
            }
        }
    }

    ntrial <- max(n * samplefactor, 5)
    ####################################
    if (method == 'SRS') {
        if (is.null(frame)) {
            origins <- coordinates(spsample(region, ntrial, type='random'))
        }
        else {
            if (ntrial > nrow(frame))
                stop ("too few rows in frame for requested sample")
            OK <- sample.int(nrow(frame), ntrial, replace = FALSE)
            origins <- as.matrix(frame[OK, ])
        }
    }
    ####################################
    else if (method == 'GRTS') {
        if (!require (spsurvey))
            stop ("package 'spsurvey' required for grts in trap.builder")
        ## make a list in the format needed by grts()
        design <- list(None = list(panel = c(Panel1 = n),
            seltype = "Equal", over = ntrial))
        src <-ifelse (is.null(frame), 'sp.object', 'att.frame')
        typ <-ifelse (is.null(frame), 'area', 'finite')
        GRTS.sites <- grts (design = design, type.frame = typ,
            src.frame = src, sp.object = region, att.frame = frame,
            shapefile = FALSE)
        origins <- coordinates(GRTS.sites)
    }
    ####################################
    else if (method == 'all') {
        if (is.null(frame)) {
            stop ("'all' requires finite frame")
        }
        origins <- as.matrix(frame)
    }
    ####################################
    else if (method == 'rank') {
        if (is.null(frame)) {
            stop ("'rank' requires finite frame")
        }
        if (is.null(ranks)) {
            stop ("'rank' requires ranks")
        }
        nframe <- nrow(frame)
        if (nframe<n)
            stop ("not enough rows in frame for requested n")
        ranks <- ranks + runif(nframe)/(nframe+1)
        frameorder <- rev((1:nframe)[order(ranks)])
        frame <- frame[frameorder,]
        origins <- as.matrix(frame)
    }
    else
        stop ("method not recognised")

    #######################################################
    ## centre cluster on (0,0)
    if (nrow(cluster)>1)
        cxy <- apply(cluster,2,mean)
    else
        cxy <- unlist(cluster)    ## assume one detector
    cluster[,] <- sweep(cluster, MARGIN=2, FUN='-', STATS=cxy)
    #######################################################

    .local$clusteri <- 1    ## updated within position()
    if (method %in% c('all','rank')) {
        ## position all, even if we will later reject some on ranks
        traps <- lapply(1:nrow(frame), position, cluster)
    }
    else {
        traps <- lapply(1:(ntrial), position, cluster)
    }
    if (edgemethod %in% c('clip', 'allowoverlap')) {
        if (n==1)
            traps <- traps[[1]]
        else {
            traps <- traps[1:n]
            traps$renumber <- FALSE
            traps <- do.call(rbind, traps)
        }
        if (edgemethod == 'clip') {
            xy <- SpatialPoints(as.matrix(traps))
            OK <- overlay (xy, region)
            traps <- subset(traps, subset = !is.na(OK))
        }
    }
    else if (edgemethod %in% c('allinside')) {
        if (is.null(region))
            stop ("allinside requires 'region'")
        OK <- sapply(traps, allinside)
        if (method == 'all')
            n <- sum(OK)
        if (sum(OK) < n)
            stop ("not enough clusters inside polygon")
        traps <- traps[OK][1:n]   ## first n usable clusters
        if (n==1)
            traps <- traps[[1]]
        else {
            traps$renumber <- FALSE
            traps <- do.call(rbind, traps)
        }
    }
    else {
        stop ("edgemethod not recognised")
    }

    ## renumber clusters
    if (attr(cluster,'detector') %in% .localstuff$polydetectors) {
        npoly <- ndetector(traps)
        npercluster <- nrow(cluster)
        polyID(traps) <- factor(rep(1:npoly, rep(npercluster, npoly)))
        clustertrap(traps) <- rep(1, nrow(traps))
        clusterID(traps) <- polyID(traps)
        vertexpart <- rep(rownames(cluster), npoly)
        row.names(traps) <- paste(polyID(traps), vertexpart, sep = '.')
    }
    else {
        clusterID(traps) <- factor(as.numeric(clusterID(traps)))
        if (nrow(cluster) == 1)
            newnames <- clusterID(traps)
        else
            newnames <- paste(clusterID(traps),
                row.names(cluster)[clustertrap(traps)], sep='.')
        row.names(traps) <- newnames
    }

    ####################################
    ## optional plot
    if (plt) {
        plot(traps, add=TRUE)
        invisible(traps)
    }
    else
        traps
    ####################################
}
###############################################################################

make.systematic <- function (n, cluster, region, spacing = NULL,
    origin = NULL, ...) {

## 'cluster' is a traps object for one module
## 'region' is a rectangular survey region
## ... arguments passed to trap.builder (rotate, detector)

    if (!require(sp))
        stop ("package 'sp' required in make.systematic")

    SPDF <- inherits(region, "SpatialPolygonsDataFrame")
    if (!SPDF) {
        ## convert to SpatialPolygonsDataFrame
        ## future: recognise & import shapefile
        region <- matrix(unlist(region), ncol = 2)
        region <- rbind (region, region[1,])  # force closure of polygon
        region <- boundarytoSPDF(region)
    }
    wd <- diff(bbox(region)[1,])
    ht <- diff(bbox(region)[2,])

    if (missing(cluster)) {
        ## this case is passed to trap builder for single detector placement
        ## if ... does not include detector, detector defaults to 'multi'
        cluster <- NULL
        clwd <- 0
        clht <- 0
    }
    else {
        clwd <- diff(range(cluster$x))
        clht <- diff(range(cluster$y))
    }

    wx <- clwd/2
    wy <- clht/2

    if (!is.null(spacing)) {
        rx <- spacing[1]
        ry <- ifelse(length(spacing)>1, spacing[2], rx)
        nx <- round ((wd-2*wx)/rx)
        ny <- round ((ht-2*wy)/ry)
    }
    else {
        if (length(n)>1) {
            nx <- n[1]
            ny <- n[2]
        }
        else {
            area <- sum(sapply(region@polygons, function(x) x@area))
            cell <- sqrt(area / n)
            nx <- round ((wd - 2*wx) / cell)
            ny <- round ((ht - 2*wy) / cell)
        }
        rx <- (wd - 2*wx) / nx
        ry <- (ht - 2*wy) / ny
    }
    rxy <- c(rx,ry)
    if (is.null(origin))
        origin <- runif(2) * rxy + bbox(region)[,1] + c(wx,wy)
    else {
        origin <- origin + rxy * trunc((bbox(region)[,1] - origin) / rxy)
    }
    centres <- expand.grid (
        x = seq(0, by = rx, len = nx) + origin[1],
        y = seq(0, by = ry, len = ny) + origin[2])
    centres <- SpatialPoints(as.matrix(centres))
    OK <- !is.na(overlay (centres, region))
    centres <- coordinates(centres[OK,])
    trap.builder (cluster = cluster, frame = centres, region = region,
        method = 'all', ...)
}

###############################################################################

mash <- function(object, origin = c(0,0), clustergroup = NULL, ...) {

## mash() recasts a capthist object in which the detectors belong to
## multiple clusters as a capthist with multiple detections at one cluster
## This assumes independence of clusters: if any individuals were detected
## on multiple clusters their new detection histories will be misleading

    if (is.list(clustergroup) & (length(clustergroup) > 1)) {
        if (ms(object))
            stop ("cannot regroup multisession capthist")
        out <- vector('list')
        for (i in 1:length(clustergroup)) {
            out[[i]] <- mash (object, origin, clustergroup[[i]])
        }
        names(out) <- names(clustergroup)
        class(out) <- c('list', 'capthist')
        if (length(out) == 1) out <- out[[1]]
        return(out)
    }
    else if (ms(object)) {
        out <- lapply(object, mash, origin, clustergroup)
        names(out) <- names(clustergroup)
        class(out) <- c('list', 'capthist')
        if (length(out) == 1) out <- out[[1]]
        return(out)
    }
    else {
        if (!is.null(clustergroup)) {
            trapsi <- clusterID(traps(object)) %in% clustergroup
            object <- subset(object, traps = trapsi)
        }
        trps <- traps(object)
        if (!is.null(covariates(trps)))
            warning ("detector covariates are discarded by mash()")
        if (!is.null(usage(trps)))
            warning ("usage discarded by mash()")
        cluster <- clusterID(trps)
        centres <- cluster.centres(trps)

        ## how many individuals per cluster?
        ## assign each to the first cluster in which it appears
        cl <- cluster[trap(object, names = FALSE)]
        ID <- animalID(object, names = FALSE)
        n.mash <- table (cl[match(unique(ID),ID)])

        if (is.null(cluster))
            stop ("requires cluster covariate")
        tmp <- split(trps, cluster)
        if (length(unique(sapply(tmp, nrow))) != 1)
            warning ("unequal number of detectors per cluster")

        trapnum <- clustertrap(trps)
        if (is.null(trapnum)) {
            tmp <- lapply(tmp, function(x) {x$trapnum <- 1:nrow(x); x})
            trapnum <- unlist(sapply(tmp, function(x) x$trapnum))
        }

        ## take first cluster for new traps
        newtraps <- tmp[[1]]
        rownames(newtraps) <- 1:nrow(newtraps)
        mxy <- apply(newtraps, 2, min)
        newtraps <- shift(newtraps, origin-mxy[1:2])

        sigcov <- NULL
        if ( length(animalID(object)) == 0) {
            tempdf <- data.frame(
                session = session(object),
                ID = 'NONE',
                occ = ncol(object),
                trap = 1)
        }
        else {
            tempdf <- data.frame(
                session = rep(session(object), length(animalID(object))),
                ID = animalID(object),
                occ = occasion(object),
                trap = trapnum[trap(object, names=FALSE)]
            )
            if (!is.null(attr(object, 'signalframe'))) {
                tempdf <- cbind(tempdf, attr(object, 'signalframe'))
                sigcov <- names(tempdf)[!(names(tempdf) %in% c('signal','noise'))]
            }
        }
        tempcapt <- make.capthist(tempdf, newtraps, cutval = attr(object, "cutval"),
                                  signalcovariates = sigcov, ...)
        attr(tempcapt, 'n.mash') <- as.numeric(n.mash)
        attr(tempcapt, 'centres') <- centres
        tempcapt
    }
}

###############################################################################

cluster.counts <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires capthist object")
    clust <- clusterID(traps(object))
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(traps(object)))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    cl <- clust[trap(object, names = FALSE)]
    tmp <- data.frame(ID=animalID(object), cluster = cl)
    sapply(split(tmp,tmp$cluster), function(x) length(unique(x$ID)))
}
###############################################################################

cluster.centres <- function (object) {
    if (!inherits(object, 'traps'))
        stop ("requires traps object")
    clust <- clusterID(object)
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(object))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    data.frame(x = tapply(object$x,clust,mean),
               y = tapply(object$y,clust,mean))
}
###############################################################################
