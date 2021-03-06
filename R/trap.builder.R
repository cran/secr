###############################################################################
## package 'secr'
## trap.builder.R
## repeat trap layout systematically, by GRTS, or at simple random centres
## across a region
## Also make.systematic, mash(), cluster.counts(), cluster.centres()
## 2011-08-16 (full argument names)
## 2012-01-11 (cutval and signal in mash(); mxy[1:2])
## 2012-02-07 mash noise
## 2012-07-26 mash cleans out trap attributes; names sessions for list input
## 2012-09-17 mash moved to mash.r
## 2013-03-02 exclude, exclmethod arguments added
## 2013-04-20 replace deprecated overlay with over
## 2014-10-25 revamp polygon requirement for grts
## 2014-12-11 set proj4string to NA
## 2014-12-11 method = "GRTS" changed to method == "GRTS" !
## 2018-09-27 coerce CRS of region to CRS()
## 2018-12-18 clipping left to trap.builder()
## 2018-12-19 make.systematic() grid of centres expanded
## 2018-12-29 boundarytoSPDF and boundarytoSP moved to utility.R
## 2018-12-29 make.systematic() argument 'chequerboard'
## 2019-12-20 make.systematic() arguments 'rotate', 'centrexy'
## 2020-01-27 fix bug (saved n = NULL when n missing)
###############################################################################

## spsurvey uses sp

trap.builder <- function (n = 10, cluster, region = NULL, frame =
    NULL, method = c("SRS", "GRTS", "all", "rank"), edgemethod =
    c("clip", "allowoverlap", "allinside", "anyinside", "centreinside"), 
    samplefactor = 2, ranks = NULL, rotation = NULL, detector, 
    exclude = NULL, exclmethod = c("clip", "alloutside", "anyoutside", "centreoutside"), 
    plt = FALSE, add = FALSE) {
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

    .local <- new.env()   ## for clusteri
    method <- match.arg(method)
    edgemethod <- match.arg(edgemethod)
    exclmethod <- match.arg(exclmethod)
    allinside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy), proj4string = CRS())
        ## 2014-10-25 polygons() works with both SP and SPDF
        !any(is.na(sp::over (xy, polygons(region))))
    }

    anyinside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy), proj4string = CRS())
        any(!is.na(sp::over (xy, polygons(region))))
    }

    centreinside <- function (xy) {
        xy <- apply(as.matrix(xy),2,mean)
        xy <- SpatialPoints(matrix(xy, ncol=2), proj4string = CRS())
        !is.na(sp::over (xy, polygons(region)))
    }

    alloutside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy), proj4string = CRS())
        ## 2014-10-25 polygons() works with both SP and SPDF
        all(is.na(sp::over (xy, polygons(exclude))))
    }

    anyoutside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy), proj4string = CRS())
        any(is.na(sp::over (xy, polygons(exclude))))
    }

    centreoutside <- function (xy) {
        xy <- apply(as.matrix(xy),2,mean)
        xy <- SpatialPoints(matrix(xy, ncol=2), proj4string = CRS())
        is.na(sp::over (xy, polygons(exclude)))
    }

    position <- function (i, cluster) {
        #newtraps <- secr::shift(cluster, origins[i,])
        newtraps <- shift(cluster, origins[i,])
        if (!is.null(rotation)) {
            if (rotation<0)
                rotation <- runif(1) * 360
            #bnewtraps <- secr::rotate(newtraps, rotation, apply(newtraps,2,mean))
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
        ## 2019-09-25 block this 
        ## edgemethod <- 'allowoverlap'
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

    SP <- inherits(region, 'SpatialPolygons')
    SPx <- inherits(exclude, 'SpatialPolygons')

    if (is.null(frame)) {
        if (is.null(region)) {
            stop ("specify at least one of 'region' or 'frame'")
        }

        if (SP) {
            ## spsurvey requires SPDF
            if ((method == 'GRTS') & (!inherits(region, 'SpatialPolygonsDataFrame'))) {
                attr <- data.frame(a = 1, row.names = "s1")
                region <- SpatialPolygonsDataFrame(region, attr)
            }
        }
        else {
            region <- matrix(unlist(region), ncol = 2)
            region <- rbind (region, region[1,])  # force closure of polygon
            region <- boundarytoSPDF(region)
        }
        
        if (!is.null(exclude)) {
            exclude <- boundarytoSP(exclude)
        }
        
        if (plt & !add) {
            plot(region)
            if (!is.null(exclude))
                plot(exclude, col='lightgrey', add=TRUE, border='lightgrey')
        }
    }
    else {
        if (plt & !add) {
            if (!is.null(region)) {
                plot(region)
            }
            else {
                MASS::eqscplot (frame, axes = F, xlab = '', ylab = '', pch = 1, cex = 0.5)
            }
            if (!is.null(exclude))
                plot(exclude, col = 'lightgrey', add = T)
        }
    }

    ntrial <- max(n * samplefactor, 5)
    ## 2018-09-27
    if (!is.null(region))  proj4string(region) <- CRS()
    if (!is.null(exclude)) proj4string(exclude) <- CRS()
    
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
        if (!requireNamespace ('spsurvey', quietly = TRUE))
            stop ("package 'spsurvey' required for grts in trap.builder")
        ## make a list in the format needed by grts()
        design <- list(None = list(panel = c(Panel1 = n),
            seltype = "Equal", over = ntrial))
        src <-ifelse (is.null(frame), 'sp.object', 'att.frame')
        typ <-ifelse (is.null(frame), 'area', 'finite')
        GRTS.sites <- spsurvey::grts (design = design, type.frame = typ,
            src.frame = src, sp.object = region, att.frame = frame,
            shapefile = FALSE)
        origins <- coordinates(GRTS.sites)
    }
    ####################################
    else if (method == 'all') {
        if (is.null(frame) || (nrow(frame)==0)) {
            stop ("method = 'all' requires finite frame")
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
    else {
        stop ("method not recognised")
    }
    rownames(origins) <- 1:nrow(origins)
    
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
    ## subset whole clusters
    if ((edgemethod %in% c('allooutside', 'allinside', 'anyinside', 'centreinside')) | 
        (exclmethod %in% c('alloutside', 'anyoutside', 'centreoutside'))) {
        
        if (edgemethod %in% c('allinside', 'anyinside', 'centreinside')) {
            if (is.null(region)) stop ("'edgemethod' requires 'region'")
            OK <- sapply(traps, edgemethod) 
        }
        else OK <- rep(TRUE, length(traps)) ## 2018-12-18
        
        if (!is.null(exclude) & (exclmethod %in% c('alloutside', 'anyoutside', 'centreoutside')))
            OK <- OK & sapply(traps, exclmethod)

        if (method == 'all')
            n <- sum(OK)
        if (sum(OK) < n)
            stop ("not enough clusters inside polygon")
        
        traps <- traps[OK]  ## subset whole clusters
    }
    traps <- traps[1:n]   ## first n usable clusters
    ## convert list of clusters to flat traps object
    if (n == 1)
        traps <- traps[[1]]
    else {
        traps$renumber <- FALSE
        traps <- do.call(rbind, traps)
    }
    ## drop excluded sites, if requested
    if (edgemethod %in% c('clip', 'centreinside')) {
        xy <- SpatialPoints(as.matrix(traps), proj4string = CRS())
        OK <- sp::over (xy, polygons(region))
        traps <- subset(traps, subset = !is.na(OK))
    }
    if (!is.null(exclude) & (exclmethod %in% c('clip', 'centreoutside'))) {
        xy <- SpatialPoints(as.matrix(traps), proj4string = CRS())
        notOK <- sp::over (xy, polygons(exclude))
        traps <- subset(traps, subset = is.na(notOK))
    }

    ## renumber clusters
    oldnames <- unique(clusterID(traps))
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
    attr(traps, 'centres') <- origins[oldnames,]

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
    origin = NULL, originoffset = c(0,0), chequerboard = c('all','black','white'), 
    order = c('x', 'y', 'xb', 'yb'), 
    rotate = 0, centrexy = apply(bbox(region),1,mean), 
    keep.design = TRUE, ...) {
    
    # NOT GOOD IF REGION IS MATRIX NOT SPDF

    ## 'cluster' is a traps object for one module
    ## 'region' is a rectangular survey region
    ## ... arguments passed to trap.builder (rotate, detector)
    temporigin <- origin
    chequerboard <- match.arg(chequerboard)
    order <- match.arg(order)
    SP <- inherits(region, "SpatialPolygons")
    if (SP) {
        region <- polygons(region)
    }
    else{
        ## convert to SpatialPolygons
        ## future: recognise & import shapefile
        region <- matrix(unlist(region), ncol = 2)
        region <- rbind (region, region[1,])  # force closure of polygon
        region <- boundarytoSP(region)
    }
    
    if (rotate != 0) {
        if (requireNamespace('maptools', quietly = TRUE)) {
            region <- maptools::elide(region, rotate = -rotate, center = centrexy)
        }
        else {
            warning("install maptools for rotate option in make.systematic")
        }
    }
    
    ## 2018-09-27
    proj4string(region) <- CRS()

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
        nx <- round ((wd-2*wx)/rx + max(spacing)/rx)  ## extra to ensure coverage 2018-10-13
        ny <- round ((ht-2*wy)/ry + max(spacing)/ry)  ## extra to ensure coverage 2018-10-13
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
    
    minxy <- bbox(region)[,1]
    if (is.null(origin)) {
        origin <- runif(2) * rxy + minxy + originoffset
        originbox <- data.frame(
            x = minxy[1] + originoffset[1] + c(0,0,rx,rx),
            y = minxy[2] + originoffset[2] + c(0,ry,ry,0))
        ## originoffset = c(wx,wy) was standard before 3.1.8 2018-10-13
    }
    else {
        origin <- unlist(origin)  ## 2018-12-28
        ## unnecessary 2018-12-29
        ## origin <- origin + rxy * floor((minxy - origin) / rxy)
        originbox <- NULL
    }

    ## 2018-12-29
    rowcol <- expand.grid(row = 0:(ny+1), col = 0:(nx+1))
    # centres <- data.frame(x = rowcol$col * ry + origin[1],
    #                       y = rowcol$row * rx + origin[2])
    ## 2019-12-17 apply rx to columns, ry to rows
    centres <- data.frame(x = rowcol$col * rx + origin[1],
                          y = rowcol$row * ry + origin[2])
    ##-----------------------------------------------------------------
    ## 2019-09-28  alternative cluster sequence
    nx2 <- nx+2; ny2 <- ny+2
    if (order == 'y') temp <- 1:nrow(centres)
    if (order == 'x') temp <- t(matrix(1:nrow(centres), ncol = nx2))
    if (order == 'yb') {
        temp <- matrix(1:(nx2*ny2), ncol = ny2)
            for (i in seq(2,ny2,2)) temp[,i] <- rev(temp[,i])
    }
    if (order == 'xb') {  
        temp <- t(matrix(1:(nx2*ny2), ncol = nx2))
        for (i in seq(2,nx2,2)) temp[i,] <- rev(temp[i,])
    }
    row.names(centres) <- as.numeric(temp)
    centres <- centres[order(temp),]   
    ##-----------------------------------------------------------------
    
    if (chequerboard != 'all') {
        whitesquares <- trunc(rowcol$row + rowcol$col + 0.1) %% 2L == 1L
        if (chequerboard == 'white')
            centres <- centres[whitesquares,]
        else 
            centres <- centres[!whitesquares,]
    }
    
    # blocked 2018-12-18; restored 2018-12-24 because NULL clusters not clipped in trap.builder
    args <- list(...)
    if (!is.null(args$edgemethod)) {
        if (args$edgemethod %in% c('allinside', 'centreinside')) {
            spcentres <- SpatialPoints(as.matrix(centres), proj4string = CRS())
            OK <- !is.na(sp::over (spcentres, region))
            centres <- coordinates(spcentres[OK,])
        }
    }
    
    traps <- trap.builder (cluster = cluster, frame = centres, region = region,
        method = 'all', ...)
    
    if (keep.design) {
        design <- list (
            'function' = 'make.systematic',
            n = if (missing(n)) NULL else n, 
            cluster = cluster, 
            region = region, 
            spacing = spacing,
            origin = temporigin, 
            originoffset = originoffset, 
            chequerboard = chequerboard,
            order = order, 
            rotate = rotate, 
            centrexy = centrexy)
        design <- c(design, list(...))
        attr(traps, 'design') <- design
    }
    else {
        attr(traps, 'origin') <- origin
        attr(traps, 'originbox') <- originbox
    }
    
    if (rotate != 0 && requireNamespace('maptools', quietly = TRUE)) {
        ## reverse rotation applied to region polygon
        ## silently as any maptools warning issued when rotate first encountered
        traps <- rotate(traps, rotate, centrexy)
    }
    
    traps
}

###############################################################################
make.lacework <- function (region, spacing = c(100, 20), times = NULL, 
                           origin = NULL, rotate = 0, 
                           radius = NULL, detector = "multi", keep.design = TRUE) {
    spacing <- as.numeric(spacing)
    if (length(spacing) == 1) {
        if (is.null(times)) stop("make.lacework requires times if spacing length 1")
        spacing <- c(spacing * times, spacing)
    }
    temporigin <- origin
    fraction <- 1.0  ## suspended code
    if (length(spacing) != 2) {
        stop("invalid spacing 2-vector in make.lacework")
    }
    if (inherits(region, 'SpatialPolygons')) {
        region <- sp::polygons(region)
    }
    else {
        region <- matrix(unlist(region), ncol = 2)
        region <- rbind (region, region[1,])  # force closure of polygon
        region <- boundarytoSP(region)
    }
    if (rotate != 0) {
        if (requireNamespace('maptools', quietly = TRUE)) {
            centrexy <- apply(sp::bbox(region),1,mean)
            region <- maptools::elide(region, rotate = -rotate, center = centrexy)
        }
        else {
            warning("rotation requires package maptools; not rotated")
        }
    }    
    bbox <- sp::bbox(region)                    ## after rotation
    if (is.null(origin)) {
        origin <- bbox[,1]
        origin <- origin - runif(2)*spacing[1]
    }
    rx <- diff(bbox[1,])
    ry <- diff(bbox[2,])
    n1 <- ceiling((rx + spacing[1]) / spacing[1]) + 2
    n2 <- ceiling((ry + spacing[1]) / spacing[2]) + 2
    gridx <- make.grid(nx=n1, ny=n2, spacex = spacing[1], spacey = spacing[2], originxy = origin, 
                       detector = detector, ID = 'numxb')
    n1 <- ceiling((ry + spacing[1]) / spacing[1]) + 2
    n2 <- ceiling((rx + spacing[1]) / spacing[2]) + 2
    gridy <- make.grid(ny=n1, nx=n2, spacey = spacing[1], spacex = spacing[2], originxy = origin, 
                       detector = detector, ID = 'numyb')
    if (fraction < 1) {
        OKx <- ((gridx$y-origin[2]) %% spacing[1]) < fraction * spacing[1]
        OKy <- ((gridy$x-origin[1]) %% spacing[1]) < fraction * spacing[1]
        gridx <- subset(gridx,OKx)
        gridy <- subset(gridy,OKy)
    }
    grid <- rbind(gridx, gridy, renumber = FALSE, checkdetector = FALSE)
    dupl <- duplicated(round(grid,5))
    crossings <- subset(grid, dupl)
    grid <- subset(grid, !dupl)
    grid <- subset(grid, pointsInPolygon(grid, region))
    crossings <- subset(crossings, pointsInPolygon(crossings, region))
    if (!is.null(radius)) {
        grid <- subset(grid, distancetotrap(grid, crossings)<=radius)
    }
    rownames(grid) <- sapply(lapply(strsplit(rownames(grid), ".", TRUE), rev), paste, collapse='.')
    if (rotate != 0) {
        grid <- rotate(grid, degrees = rotate, centrexy = centrexy)
        crossings <- rotate(crossings, degrees = rotate, centrexy = centrexy)
    }
    attr(grid, 'crossings') <- as.matrix(crossings)
    if (keep.design) {
        design <- list (
            'function' = 'make.lacework',
            region = region, 
            spacing = spacing,
            origin = temporigin, 
            rotate = rotate, 
            radius = radius,
            detector = 'multi')
        attr(grid, 'design') <- design
    }
    
    grid
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

# lacework test
# library(secr)
# bx <- data.frame(x=c(0,0,270,270,0), y = c(0,270,270,0,0))
# plot(make.lacework(bx, c(20,5)), border=10)
# lines(bx)