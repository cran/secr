###############################################################################
## package 'secr'
## methods.R
## Methods for classes traps, capthist and mask
## Last changed
## 2010 04 27 version 1.4
## 2010 04 30 scrapped read.captures
## 2010 05 02 usage may be separated by commas in read.traps
## 2010 05 03 bugfix in subset.capthist for 'multi' data
## 2010 06 08 sim.popn moved to own file sim.popn.R
## 2010 07 01 alphanumeric detectfn option
## 2010 06 08 fixed bug in insidepoly (wrong call of C fn)
## 2010 07 03 rbind.mask handles covariates
## 2010 08 30 print.mask removed
## 2010 10 10 enabled compound halfnormal in secr.fit
## 2010 10 10 added function 'detectpar'
## 2010 10 11 secr.fit default binomN changed to Poisson if detector == count
## 2010 10 12 covariates<- modified for multi-session data
## 2010 10 12 plot.mask now handles multi-session mask
## 2010 10 15 function ms to recognise multisession objects
## 2010 10 15 make.mask poly clipping extended to all types
## 2010 10 19 version 1.5
## 2010-10-21 AIC.secrlist
## 2010-10-22 predict.secrlist
## 2010-10-24 multisession detectpar
## 2010-10-25 revision of all error and warning messages in R code
## 2010-11-04 logLik.secr
## 2011-01-04 G&R parameterisation
## 2011-01-12 cue detector
## 2011-02-01 remove binomN from traps object attribute list
## 2011-02-07 major revisions for polygonX, transectX
## 2011-03-18 unmarked detector
##            for counts of unidentified animals at points
## 2011-06-07 rbind.capthist allows for xy, signal
## 2011-06-17 moved insidepoly to utility.R (renamed to pointsInPolygon 2011-10-21)
## 2011-06-23 polygon and transect names retain ordering of input dataframe
##            in make.poly and make.transect
## 2011-06-24 fixed transectlength, searcharea
## 2011-07-08 adjust 'class' to 'inherits' in make.mask
## 2011-08-18 check for null rownames in animalID
## 2011-09-11 subset.popn commissioned
## 2011-10-10 make.mask moved to make.mask.R
## 2011-10-10 make.grid, make.poly etc. moved to make.grid.R
## 2011-10-20 predict.secr modified for user-supplied D function
## 2011-10-20 subset.popn gets poly argument
## 2012-05-14 predict.secr patch for missing 'session' in newdata (cf Deb's problem 7/2/12)
###############################################################################

# Generic methods for extracting attributes etc

usage      <- function (object, ...) UseMethod("usage")
clusterID  <- function (object, ...) UseMethod("clusterID")
clustertrap <- function (object, ...) UseMethod("clustertrap")
covariates <- function (object, ...) UseMethod("covariates")
traps      <- function (object, ...) UseMethod("traps")
detector   <- function (object, ...) UseMethod("detector")
spacing    <- function (object, ...) UseMethod("spacing")
session    <- function (object, ...) UseMethod("session")
trim       <- function (object, drop, keep) UseMethod("trim")

reduce     <- function (object, columns, ...) UseMethod("reduce")
rotate     <- function (object, degrees, centrexy=NULL, ...) UseMethod("rotate")
shift      <- function (object, shiftxy, ...) UseMethod("shift")
flip       <- function (object, lr=F, tb=F, ...) UseMethod("flip")

ms         <- function (object, ...) UseMethod("ms")

# Default methods for specialised functions

ms.default <- function (object, ...)       {
    inherits(object, 'list')
}
ms.mask <- function (object, ...)       {
    !is.data.frame(object)
}
ms.secr <- function (object, ...)       {
    ms(object$capthist)
}

usage.default <- function (object, ...)       {
    if (ms(object)) lapply(object, usage.default, ...)
    else attr(object,'usage')
}

clusterID.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clusterID.default, ...)
    else attr(object,'cluster')
}

clustertrap.default <- function (object, ...)       {
    if (ms(object)) lapply(object, clustertrap.default, ...)
    else attr(object,'clustertrap')
}

covariates.default <- function (object, ...)  {
    if (ms(object)) lapply(object, covariates.default, ...)
    else attr(object,'covariates')
}

traps.default <- function (object, ...)       {
    if (ms(object)) {
        temp <- lapply(object, traps.default, ...)
        class(temp) <- c('list', 'traps')
        temp
    }
    else{
        attr(object,'traps')
    }
}

detector.default <- function (object, ...)    {
## assumed constant across MS
    if (ms(object)) {
        detector.default (object[[1]], ...)
    }
    else {
        if (is.null(object)) NULL
        else attr(object,'detector')
    }
}

spacing.default <- function (object, ...)    {
    if (is.null(object))
        NULL
    else {
        attr(object,'spacing')
    }
}

spacing.traps <- function (object, ...)    {
    if (ms(object)) {
        sapply(object, spacing.traps, ...)
    }
    else {
        if (is.null(object)) {
            NULL
        }
        else {
            temp <- attr(object,'spacing')
            if (is.null(temp) & (nrow(object)>1)) {
                spacing <- as.matrix(dist(object))
                sp <- apply(spacing,1,function(x) min(x[x>0]))
                mean(sp)
            }
            else
                temp
        }
    }
}

spacing.mask <- function (object, ...)    {
    if (ms(object)) {
        sapply(object, spacing.mask, ...)
    }
    else {
        if (is.null(object)) NULL
        else {
            temp <- attr(object,'spacing')
            if (is.null(temp) & (nrow(object)>1) ) {
                spacing <- as.matrix(dist(object))
                sp <- apply(spacing,1,function(x) min(x[x>0]))
                mean(sp)
            }
            else
                temp
        }
    }
}

## traps object
polyID <- function (object)    {
    if (ms(object)) {
        polyID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            temp <- attr(object,'polyID')
#            if (is.null(temp)) temp <- factor(rep(1,nrow(object)))
            if (is.null(temp)) temp <- factor(1:nrow(object))   ## all different
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract polyID from 'capthist' object")
        }
        else stop ("polyID requires 'traps' object")
    }
}

## traps object
transectID <- function (object)    {
    if (ms(object)) {
        transectID (object[[1]])
    }
    else {
        if (inherits(object,'traps')) {
            if (!detector(object) %in% c('transect','transectX'))
                stop ("requires transect detector")
            temp <- attr(object,'polyID')
#            if (is.null(temp)) temp <- factor(rep(1,nrow(object)))
            if (is.null(temp)) temp <- factor(1:nrow(object))
            temp
        }
        else
        if (inherits(object,'capthist')) {
            stop ("use trap() to extract transectID from 'capthist' object")
        }
        else stop ("transectID requires 'traps' object")
    }
}

xy <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, xy)
    }
    else {
        if (detector(traps(object)) %in%
            c('polygonX', 'transectX', 'polygon','transect')) {
            attr(object, 'detectedXY')
        }
        else
            NULL
    }
}

alongtransect <- function (object, tol = 0.01) {
    ptalongtransect <- function (i) {
        ## where is point i on its transect k?
        k <- trans[i]
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        temp <- .C('alongtransect',  PACKAGE = 'secr',
            as.double (xyi[i,]),
            as.integer (0),
            as.integer (nr-1),
            as.integer (nr),
            as.double (transectxy),
            as.double (tol),
            result = double(1))
        temp$result
    }
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, alongtransect)
    }
    else {
        trps <- traps(object)
        if (detector(trps) %in% c('transectX', 'transect')) {
            trans <- trap(object, names = TRUE)
            xyi <- xy(object)
            lxy <- split (trps, levels(transectID(trps)))
            sapply(1:nrow(xyi), ptalongtransect)
        }
        else
            NULL
    }
}

clusterID <- function (object) {
    if (ms(object)) {
        lapply(object, clusterID)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clusterID(trps)[trap(object, names=FALSE)]
        }
        else
            attr(object, 'cluster')
    }
}

clustertrap <- function (object) {
    if (ms(object)) {
        lapply(object, clustertrap)
    }
    else {
        if (inherits(object, 'capthist')) {
            trps <- traps(object)
            clustertrap(trps)[trap(object, names=FALSE)]
        }
        else
        attr(object, 'clustertrap')
    }
}

signal <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, signal)
    }
    else {
        if (detector(traps(object)) %in% c('cue','signal','signalnoise')) {
            attr(object, 'signalframe')$signal
        }
        else
            NULL
    }
}

noise <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, noise)
    }
    else {
        if (detector(traps(object)) %in% c('signalnoise')) {
            attr(object, 'signalframe')$noise
        }
        else
            NULL
    }
}

times <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")

    if (ms(object)) {
        lapply(object, times)
    }
    else {
        if (detector(traps(object)) %in% c('times')) {
            attr(object, 'times')
        }
        else
            NULL
    }
}

occasion <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, occasion)
    }
    else {
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            col(object)[abs(object)>0]
        }
        else {
            apo <- aperm(object,c(2,1,3))
            temp <- matrix(apo, nrow = dim(apo)[1])
            temp <- array(row(temp), dim = dim(apo))
            rep(aperm(temp,c(2,1,3)), abs(object))
        }
    }
}

alive <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, alive)
    }
    else
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            sign(object)[abs(object)>0] > 0
        }
        else {
            temp <- sign(object)[abs(object)>0] > 0
            rep(temp, abs(object[abs(object)>0]))
        }
}

trap <- function (object, names = TRUE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, trap, names=names)
    }
    else {
        if (names)
            values <- row.names(traps(object))
        else
            values <- 1:nrow(traps(object))
        if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
            values[abs(object[abs(object)>0])]
        }
        else {
            apo <- aperm(object,c(3,1,2))
            temp <- matrix(apo, nrow = dim(apo)[1])
            temp <- array(row(temp), dim = dim(apo))
            k <- aperm(temp,c(2,3,1))
            rep(values[k], abs(object))
        }
    }
}

animalID <- function (object, names = TRUE) {
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, animalID, names=names)
    }
    else {
        if (nrow(object) == 0)
           ## ''
           character(0)  ## 2011-04-08
        else {
            if (names & !is.null(row.names(object)))  ## 2011-08-18 null check
                values <- row.names(object)
            else
                values <- 1:nrow(object)
            if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
                detrow <- row(object)[abs(object)>0]
                values[detrow]
            }
            else {
                temp <- matrix(object, nrow = dim(object)[1])
                n <- array(row(temp), dim=dim(object))
                rep(values[n], abs(object))
            }
        }
    }
}

detectionindex <- function (object) {
## detectionindex is non-exported function 2012-02-11
## to which original cell in dim3 capthist object does a detection relate?
    if (!inherits(object, 'capthist'))
        stop ("requires 'capthist' object")
    if (ms(object)) {
        lapply(object, detectionindex)
    }
    else {
        if (nrow(object) == 0)
           numeric(0)
        else {
            if (detector(traps(object)) %in% .localstuff$exclusivedetectors) {
                stop("indices is not intended  for exclusive detectors")
            }
            else {
                x <- object
                x[] <- 1:length(x)
                x[abs(object)==0] <- 0
                x <- aperm(x, c(3,2,1))
                x <- x[x>0]
                count <- aperm(object, c(3,2,1))
                count <- count[abs(count)>0]
                rep(x,count)
            }
        }
    }
}

polyarea <- function (xy, ha = TRUE) {
    if (!require(sp))
        stop ("package 'sp' required for polyarea")
    if (inherits(xy, 'SpatialPolygons')) {
        if (!require(rgeos))
            stop ("package rgeos is required for area of SpatialPolygons")
        temparea <- gArea(xy)
    }
    else {
        nr <- length(xy$x)
        ## beat problem with integer overflow 2010-11-21
        xy$x <- as.double(xy$x)
        xy$y <- as.double(xy$y)
        sc <- sum (xy$x[-nr] * xy$y[-1] - xy$x[-1] * xy$y[-nr])
        temparea <- abs(sc/2)
    }
    ifelse (ha, temparea/10000, temparea)
}

searcharea <- function (object)    {
## use search cell area if attribute present 2009 11 08
## searchcell, or spacex.spacey, or spacing^2, or 1
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, searcharea)
    }
    else {
        if (detector(object) %in% c('polygon','polygonX')) {
            ## assume poly closed
            sapply(split(object, polyID(object)), polyarea)
        }
        else {
            sc <- attr(object,'searchcell')
            if (is.null(sc)) {
                spx <- attr(object,'spacex')
                spy <- attr(object,'spacey')
                if (is.null(spx) | is.null(spy))
                    sp2 <- attr(object,'spacing')^2
                else
                    sp2 <- spx * spy
                ifelse (is.null(sp2), 1, sp2 / 10000)
            }
            else sc
        }
    }
}

transectlength <- function (object)    {
## requires traps object
    if (ms(object)) {
        ## 2011-06-24
        lapply(object, transectlength)
    }
    else {
        if (detector(object) %in% c('transect','transectX')) {
            calclength <- function (xy) {
                nr <- nrow(xy)    ## number of vertices
                segments <- (xy$x[-nr] - xy$x[-1])^2 + (xy$y[-nr] - xy$y[-1])^2
                sum (segments^0.5)
            }
            sapply(split(object, polyID(object)), calclength)
        }
        else NA
    }
}

session.default <- function (object, ...)     {
## bypass session attribute for multi-session objects 2009 12 22
## use names(object) for lists i.e. multi-session objects

    if (ms(object)) {
       temp <- names(object)
    }
    else {
        temp <- attr(object,'session')
        if (is.null(temp)) temp <- 1    ## added 2010-02-03
    }
    names(temp) <- NULL
    as.character(temp)      ## 2010 02 25
}

trim.default <- function (object, drop, keep)     {
## drop unwanted named components of a list
## conservative resolution of conflicts between drop & keep
    objnames <- names(object)
    indices <- 1:length(object)
    if (missing(drop)) drop <- indices  # all! but wait...
    if (missing(keep)) keep <- 0
    if (is.character(keep)) keep <- match(keep, objnames)
    if (is.character(drop)) drop <- match(drop, objnames)
    drop <- drop[drop %in% indices[! (indices %in% keep)]]
    ## by index, so have to work from end
    for (i in sort(drop, decreasing = T)) object[[i]] <- NULL
    object
}

reduce.default <- function (object, columns, ...) {
  object <- as.matrix(object)
  if (any(is.na(object)))
      warning ("NAs in input converted to zero")
  firsttrap <- function (y) y[abs(y)>0][1]    # first non-zero
  fnmulti   <- function (occ) apply (object[,occ,drop=F], 1, firsttrap)
  nrow <- nrow(object)
  nnew <- length(columns)
  temp <- sapply (columns, fnmulti)
  temp[is.na(temp)] <- 0
  temp
}
###############################################################################

rotate.default <- function (object, degrees, centrexy=NULL, ...) {

    rotatefn <- function (xy) {
        # about centre
        x <- xy[1] - centrexy[1]
        y <- xy[2] - centrexy[2]
        x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
        y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
        c(x2,y2)
      }

##    if (ms(object)) lapply(object, rotate.default, degrees=degrees, centrexy=centrexy, ...)
##    else

    object <- as.matrix(object)
    if (dim(object)[2] <2)
        stop ("requires at least 2 columns")
    if (abs(degrees)>0) {
        if (is.null(centrexy)) centrexy <- c(0,0)
        theta <- 2*pi*degrees/360 # convert to radians
        temp <- t(apply (object,1,rotatefn))
        if (!is.null(dimnames(object)[[2]]))
            dimnames(temp)[[2]] <- dimnames(object)[[2]][1:2]
        temp
    } else object[,1:2]
}
###############################################################################

shift.default <- function (object, shiftxy, ...) {
##    if (ms(object)) lapply(object, shift.default, shiftxy=shiftxy, ...)
##    else

    object <- as.matrix(object[,1:2])
    object[,1] <- object[,1] + shiftxy[1]
    object[,2] <- object[,2] + shiftxy[2]
    object
}
###############################################################################

flip.default <- function (object, lr=F, tb=F, ...) {
##    if (ms(object)) lapply(object, flip.default, lr=lr, tb=tb, ...)
##    else
    object <- as.matrix(object[,1:2])
    if (is.logical(lr)) {
        if (lr) object[,1] <- 2 * mean(object[,1]) - object[,1]  ## flip about mean
    } else
        if (is.numeric(lr)) object[,1] <- 2*lr-object[,1]  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object[,2] <- 2 * mean(object[,2]) - object[,2]  ## flip about mean
    } else
        if (is.numeric(tb)) object[,2] <- 2*tb-object[,2]  ## flip about tb

    object
}
###############################################################################

# Generic methods for replacing values

'usage<-' <- function (object, value) structure (object, usage = value)
'clusterID<-' <- function (object, value) structure (object, cluster = value)
'clustertrap<-' <- function (object, value) structure (object, clustertrap = value)

'covariates<-' <- function (object, value) {
## modified for multi-session data 2010-10-15
    if (is.null(value))
        structure (object, covariates = NULL)
    else {
        if (ms(object)) {
            if (length(object) != length(value))
                stop ("mismatch between multisession object and covariates")
            value <- lapply(value, as.data.frame)
        }
        else {
            value <- as.data.frame(value)
            nrequired <- nrow(object)
            if (inherits(object, 'traps'))
                if (detector(object) %in% .localstuff$polydetectors)
                    nrequired <- length(levels(polyID(object)))
            if (nrow(value) != nrequired)
                stop ("length of covariate does not match")
        }
        structure (object, covariates = value)
    }
}

'detector<-' <- function (object, value) {
    if (!(value %in% .localstuff$validdetectors))
        stop ("invalid detector type")
    structure (object, detector = value)
}

'spacing<-' <- function (object, value) {
    if (!(is.numeric(value)))
        stop ("non-numeric spacing")
    if (ms(object)) {
        stop ("not sure how to replace spacing of ms object")
    }
    else {
        structure (object, spacing = value)
    }
}

'polyID<-' <- function (object, value) {
    if (!inherits(object,'traps'))
        warning ("polyID requires 'traps' object")
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'transectID<-' <- function (object, value) {
    if (length(value)==1) value <- rep(value, nrow(object))
    value <- factor(value)
    structure (object, polyID = value)
}

'xy<-' <- function (object, value) {
    if (!is.null(value)) {
        if (nrow(value) != sum(abs(object)))
            stop ("requires one location per detection")
        if (!(detector(traps(object)) %in% .localstuff$polydetectors) |
                !(inherits(object,'capthist')))
            stop ("requires 'capthist' object with ",
                  "'polygon' or 'transect' detector")
        if (ms(object))
            stop ("requires single-session 'capthist' object")
    }
    structure (object, detectedXY = value)
}

'clusterID<-' <- function (object, value) {

    if (ms(object))
        stop ("cluster requires single-session 'traps' object")

    if (length(value)==1) value <- rep(value, nrow(object))

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one cluster per detector or detector vertex")
    }
    else
        value <- NULL

    structure (object, cluster = factor(value))
}

'clustertrap<-' <- function (object, value) {

    if (ms(object))
        stop ("clustertrap requires single-session 'traps' object")

    if (!(inherits(object, 'traps')))
        stop ("requires clustered 'traps' object")

    if (length(value) > 0) {
        if (length(value) != nrow(object))
            stop ("requires one clustertrap per detector or detector vertex")
    }
    else
        value <- NULL

    structure (object, clustertrap = value)
}

'signal<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one signal per detection")
    if (!(detector(traps(object)) %in% c('signal','signalnoise')) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signal' or 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe')
    if (is.null(sf))
        sf <- data.frame(signal = value)
    else
        sf$signal <- value
    structure (object, signalframe = sf)
}

'noise<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("requires one noise per detection")
    if (!(detector(traps(object)) %in% c('signalnoise')) |
        !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'signalnoise' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    sf <- attr(object, 'signalframe')
    if (is.null(sf))
        sf <- data.frame(noise = value)
    else
        sf$noise <- value
    structure (object, signalframe = sf)
}

'times<-' <- function (object, value) {
    if (length(value) != sum(abs(object)))
        stop ("require one time per detection")
    if (!(detector(traps(object)) == 'times') | !(inherits(object,'capthist')))
        stop ("requires 'capthist' object with 'times' detector")
    if (ms(object))
        stop ("requires single-session 'capthist' object")
    structure (object, times = value)
}

'traps<-' <- function (object, value) {
    if (!is(value,'traps'))
        stop ("'traps' object required for replacement")

    ## MODIFIED 2010 04 27
    if (ms(object)) {
        nsess <- length(object)
        temp <- vector(mode='list', nsess)
        if (nsess != length(value))
            stop ("replacement value has wrong length")
        for (i in 1:nsess) temp[[i]] <- `traps<-`(object[[i]], value[[i]])
        class(temp) <- c('list', 'traps')
        temp
    }
    else {
        structure (object, traps = value)
    }
}

'session<-' <- function (object, value) {
    if (ms(object)) {
       if (length(value) != length(object))
           stop ("invalid replacement value")
       for (i in 1:length(object)) session(object[[i]]) <- value[i]   ## 2010 03 26
       structure (object, names = as.character(value))
    }
    else {
        if (length(value) > 1)
            stop ("requires only one session name per session")
        structure (object, session = as.character(value))
    }
}
###############################################################################

######################################
## Class : traps
## defines an array of detectors
## detector may be 'single', 'multi', 'proximity' etc.
######################################

## 2011-10-10 make.grid, makepoly etc. moved to make.grid.R

rotate.traps <- function (object, degrees, centrexy=NULL, ...)
{
##    if (ms(object)) lapply(object, rotate.traps, degrees, centrexy, ...)
##    else

  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }

  if (abs(degrees)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(object$x), mean(object$y))
    theta <- 2*pi*degrees/360 # convert to radians

    traps2 <- data.frame(t(apply (object,1,rotatefn)))
    names(traps2) <- c('x','y')
    attr(traps2,'class')  <- c('traps', 'data.frame')
    detector(traps2)      <- detector(object)
    if (!is.null(usage(object)))
        usage(traps2)         <- usage(object)
    if (!is.null(covariates(object)))
        covariates(traps2)    <- covariates(object)
    if (!is.null(polyID(object)))   ## includes transectID
        polyID(traps2)    <- polyID(object)
  }
  else traps2 <- object
  traps2
}
###############################################################################

shift.traps <- function (object, shiftxy, ...)
{
##    if (ms(object)) lapply(object, shift.traps, shiftxy, ...)
##    else

  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  object
}
###############################################################################

flip.traps <- function (object, lr=F, tb=F, ...) {

##    if (ms(object)) lapply(object, flip.traps, lr, tb, ...)
##    else

    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(object$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(object$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb

    object
}
###############################################################################

rotate.popn <- function (object, degrees, centrexy=NULL, ...)
{
  rotatefn <- function (xy) {
    # about centre
    x <- xy[1] - centrexy[1]
    y <- xy[2] - centrexy[2]
    x2 <- x * cos(theta) + y * sin(theta) + centrexy[1]
    y2 <- - x * sin(theta) + y * cos(theta) + centrexy[2]
    c(x2,y2)
  }
  bbox <-  attr(object, 'boundingbox')
  if (abs(degrees)>0) {
    if (is.null(centrexy)) centrexy <- c(mean(bbox$x), mean(bbox$y))
    theta <- 2*pi*degrees/360 # convert to radians

    popn2 <- data.frame(t(apply (object,1,rotatefn)))

    object[,] <- popn2[,]
    attr(object, 'boundingbox') <- data.frame(rotate (bbox, degrees, centrexy))
  }
  object
}
###############################################################################

shift.popn <- function (object, shiftxy, ...)
{
  object$x <- object$x + shiftxy[1]
  object$y <- object$y + shiftxy[2]
  bbox <-  attr(object, 'boundingbox')
  attr(object, 'boundingbox') <- data.frame(shift(bbox, shiftxy))
  object
}
###############################################################################

flip.popn <- function (object, lr=F, tb=F, ...) {
    bbox <- attr(object, 'boundingbox')
    if (is.logical(lr)) {
        if (lr) object$x <- 2 * mean(bbox$x) - object$x  ## flip about centre
    } else
        if (is.numeric(lr)) object$x <- 2*lr - object$x  ## flip about lr

    if (is.logical(tb)) {
        if (tb) object$y <- 2 * mean(bbox$y) - object$y  ## flip about centre
    } else
        if (is.numeric(tb)) object$y <- 2*tb - object$y  ## flip about tb
    attr(object, 'boundingbox') <- data.frame(flip(bbox, lr=lr, tb=tb, ...))
    object
}
###############################################################################

print.traps <- function(x, ...) {
    if (ms(x)) {
        for (i in 1:length(x)) {
            cat('\n')
            if (!is.null(names(x)))
                cat(names(x)[i], '\n')
            print (x[[i]], ...)
        }
        invisible()
    }
    else {
        temp <- data.frame(row.names=attr(x,'row.names'), x=x$x, y=x$y)
        print(temp, ...)
    }
}
###############################################################################

subset.traps <- function (x, subset = NULL, occasions = NULL, ...) {
    # subset may be numeric index or logical

    if (ms(x)) {
        temp <- lapply(x, subset.traps, subset=subset, occasions = occasions, ...)
        class(temp) <- c('list', 'traps')
        temp
    }
    else {
        ## 2011-03-29
        if (is.null(subset)) {
            if (detector(x) %in% .localstuff$polydetectors)
                subset <- 1:length(levels(polyID(x)))
            else
                subset <- 1:nrow(x)
        }

        ## polygon & transect objects subset by whole polygons or transects
        ## 2011-01-24
        rowsubset <- subset  ## default
        if (detector(x) %in% c('polygon', 'polygonX'))
            rowsubset <- as.character(polyID(x)) %in% levels(polyID(x))[subset]
        if (detector(x) %in% c('transect', 'transectX'))
            rowsubset <- as.character(transectID(x)) %in% levels(transectID(x))[subset]
        ## apply subsetting
        temp <- x[rowsubset,,drop=F]
        class(temp) <- c('traps','data.frame')
        detector(temp) <- detector(x)

        ## 2011-05-09
        clusterID(temp) <- factor(clusterID(x)[rowsubset])
        clustertrap(temp) <- factor(clustertrap(x)[rowsubset])

        ## restore polyiD, transectID, usage, covariates
        if (detector(x) %in% c('polygon', 'polygonX'))
            polyID(temp) <- factor(polyID(x)[rowsubset])
        if (detector(x) %in% c('transect', 'transectX'))
            transectID(temp) <- factor(transectID(x)[rowsubset])

        if (!is.null(usage(x))) {
            if (is.null(occasions))
                occasions <- 1:ncol(usage(x))
            usage(temp) <- usage(x)[subset,occasions,drop=F]
        }
        if (!is.null(covariates(x)))
            covariates(temp) <- covariates(x)[subset,,drop=F]
        temp
    }
}
###############################################################################

split.traps <- function (x, f, drop = FALSE, prefix='S', ...) {
  if (!inherits(x, 'traps'))
      stop ("argument to split.traps should have class 'traps'")
  if (ms(x))
      stop ("'split.traps' not suitable for multi-session traps")
  options(warn=-1)
  f <- factor(f)

## changed 2011-04-13
##  if (any(!is.na(as.numeric(levels(f))))) {
  if (any(is.na(as.numeric(levels(f))))) {
      f <- factor(paste (prefix,f,sep=''))
      sp <- paste(prefix, levels(polyID(x)), sep='')
  }
  else {
      sp <- levels(polyID(x))
  }

  if (detector(x) %in% .localstuff$polydetectors) {
      if (length(f) > length(levels(polyID(x))))
          warning ("split factor does not match traps object")
  }

  options(warn=0)
  out <- list()
  for (i in levels(f)) {
      if (detector(x) %in% .localstuff$polydetectors) {
          temp <- subset (x, subset = (sp == i), ...)
      }
      else
          temp <- subset (x, subset = (f == i), ...)
    if (!drop | (nrow(temp)>0))
      out[[i]] <- temp
  }
  class (out) <- c('list', 'traps')
  out
}

###############################################################################

rbind.traps <- function (..., renumber = TRUE) {
# combine 2 or more traps objects
# what if multi-session?
    allargs <- list(...)
    check <- function (x) {
        if (!is(x,'traps'))
            stop ("all arguments must be 'traps' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
        if (!is.null(usage(x)))
            if (ncol(usage(x)) != ncol(usage(allargs[[1]])))
                stop ("'usage' is defined for varying number of occasions")
        if (detector(x) != detector(allargs[[1]]))
            warning ("detector types vary; using first")

    }

    if (length(allargs) <= 1)
        stop ("requires more than one traps object")
    sapply (allargs, check)
    temp <- rbind.data.frame(...)
    class(temp) <- c('traps', 'data.frame')

    ## 2010 07 05
    tempdet <- lapply(allargs, detector)
    if (length(unique(sapply(tempdet, detector))) >1 )
        ## stop ("cannot combine detector types - change with reduce()")
        ## 2011-09-26
        warning ("combining detector types - maybe change first with reduce()?")
    detector(temp) <- tempdet[[1]]

    ## 2011-02-07
    ## code for originating traps object
    ## uses 'polyID = transectID' behaviour of polyID()
    if (detector(allargs[[1]]) %in% .localstuff$polydetectors) {
        oldtrapsID <- rep(1:length(allargs), sapply(allargs, nrow))
        newpolyID <- unlist(lapply(allargs, function(x) as.character(polyID(x))))
        newpolyID <- factor(paste (oldtrapsID,newpolyID, sep='.'))
        polyID(temp) <- newpolyID
    }

    ## covariates 2010 07 05, 2010-08-28
    tempcov <- lapply(allargs, covariates)
    covariates(temp) <- do.call(rbind, tempcov)

    ## clusters  2011-04-12
    tempclus <- lapply(allargs, clusterID)
    if (!is.null(tempclus[[1]])) {
        tempclus <- lapply(tempclus, function(x) as.numeric(as.character(x)) )
        clusterID(temp) <- do.call(c, tempclus)
        clustertrap(temp) <- unlist(lapply(allargs, clustertrap))
    }
    else {
        clusterID(temp) <- NULL
        clustertrap(temp) <- NULL
    }

    ## usage
    tempusage <- lapply(allargs, usage)
    if (any(!sapply(tempusage, is.null))) {
        nocc <- unique(unlist(sapply(tempusage, ncol)))
        nr <- unlist(sapply(allargs, nrow))
        if (length(nocc)>1)
            warning ("varying number of occasions; using maximum")
        nocc <- max(nocc)
        fillmissing <- function(x, nr) {
            flush.console()
            if(is.null(x))
                matrix(1, nrow = nr, ncol = nocc)
            else {
                if (ncol(x) < nocc)
                    cbind(x, matrix(0, nrow = nr, ncol = nocc))
                else
                    x
            }
        }
        for (i in 1: length(temp))
            tempusage[[i]] <- fillmissing(tempusage[[i]], nr[i])
        usage(temp) <- do.call(rbind, tempusage)
    }

    tn <- sapply(allargs, row.names, simplify=F)


 #   if (length(unique(unlist(tn))) != length(unlist(tn))) {  # renumber
        if (renumber) {
            if (detector(temp) %in% .localstuff$polydetectors) {
                if (!is.null(polyID(temp)))  ## polyID also stores transectID
                    polyID(temp) <- factor(as.numeric(polyID(temp)))
                temp <- renamepolyrows(temp)    ## see read.traps
            }
            else
                row.names(temp) <- 1:nrow(temp)
        }
        else {
            for (i in 1:length(tn)) tn[[i]] <- paste(tn[[i]],i,sep='.')
            row.names(temp) <- unlist(tn)
        }
#    }
    if (!is.null(usage(temp)))
        if (nrow(usage(temp))>0)
            row.names(usage(temp)) <- row.names(temp)

    if (!is.null(covariates(temp))) {
        if (nrow(covariates(temp))>0)
            row.names(covariates(temp)) <- row.names(temp)
    }

    temp
}
###############################################################################

plot.traps <- function(x,
    border = 100,
    label = FALSE,
    offset=c(6,6),
    add = FALSE,
    hidetr = FALSE,
    detpar=list(),
    txtpar=list(),
    bg='white',
    gridlines=TRUE,
    gridspace=100,
    gridcol='grey',
    markused=FALSE,
    markvarying=FALSE,
    labelclusters = FALSE,
    ... )
{
#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12

    if (ms(x)) {
        lapply(x, plot.traps,
            border, label, offset,
            add, hidetr, detpar, txtpar,
            bg, gridlines, gridspace, gridcol,
            markvarying, ...)
    }
    else {

        buff <- c(-border,+border)
        offsety <- ifelse (length(offset)==2, offset[2], offset[1])
        dcol <- 'red'
        detpar <- replacedefaults (list(col=dcol, pch=3, cex=0.8), detpar)
        txtpar <- replacedefaults (list(col='blue', cex=0.7), txtpar)

        if (!is.null(usage(x))) {
            used <- apply(attr(x,'usage'),1,function(z) any(z>0))
            varying <- used * apply(attr(x,'usage'),1,function(z) any(z==0))
        }
        else {
            used  <- rep(TRUE, nrow(x))
            varying <- rep(FALSE, nrow(x))
        }
        initialpar <- par(detpar)

        if (!add) {
            par(bg=bg)
            require(MASS)
            ## axes = FALSE blocks bty = 'o' 2011-05-08
            eqscplot (x$x, x$y, xlim=range(x$x)+buff, ylim=range(x$y)+buff,
                xlab='', ylab='', type='n', axes=F, ...)
            if (gridlines) {
                xl <- range(x$x)+buff
                yl <- range(x$y)+buff
                strtx <- gridspace * floor(xl[1]/gridspace)
                strty <- gridspace * floor(yl[1]/gridspace)
                finx  <- gridspace * (floor(xl[2]/gridspace) + 1)
                finy  <- gridspace * (floor(yl[2]/gridspace) + 1)
                for (xi in seq(strtx, finx, gridspace))
                    segments(xi, strty, xi, finy, col=gridcol)
                for (yi in seq(strty, finy, gridspace))
                    segments(strtx, yi, finx, yi, col=gridcol)
            }
        }

        if (hidetr==F) {
            if (detector(x) %in% c('polygon','polygonX')) {
                templist <- split (x, levels(polyID(x)), prefix='')
                lapply(templist, function (y)
                       polygon (y$x, y$y, col=detpar$col, density=0))
                if (label) for (k in 1:length(templist)) {
                    xbar <- mean(range(templist[[k]]$x))
                    ybar <- mean(range(templist[[k]]$y))
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else if (detector(x) %in% c('transect','transectX')) {
                templist <- split (x, levels(transectID(x)), prefix='')
                lapply(templist, function (y) lines (y$x, y$y, col=detpar$col))
                if (label) for (k in 1:length(templist)) {
                    xbar <- mean(range(templist[[k]]$x))
                    ybar <- mean(range(templist[[k]]$y))
                    text (xbar+offset[1], ybar+offsety, names(templist)[k])
                }
            }
            else {
                points (x$x, x$y)
                if (markused) {
                    points (x$x[used], x$y[used], pch = 1, cex = 0.8)
                }
                if (markvarying & any(varying)) {
                    points (x$x[varying], x$y[varying], pch = 16, cex = 0.8)
                }
            }
            par(txtpar)
            if (label && !(detector(x) %in% .localstuff$polydetectors)) {
                text (x$x+offset[1], x$y+offsety, rownames(x))
            }
            if (labelclusters && !(detector(x) %in% .localstuff$polydetectors)) {
                if (is.null(clusterID(x)) | is.null(clustertrap(x)))
                    stop ("require clustered traps to label with clusterID")
                cl1 <- clustertrap(x) == 1
                text (x$x[cl1]+offset[1], x$y[cl1]+offsety, clusterID(x)[cl1])
            }
            par(initialpar)   # restore
        }
        invisible()
    }
}
###############################################################################

summary.traps <- function(object, getspacing = TRUE, ...) {

#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12

    if (ms(object)) lapply(object, summary.traps, getspacing = getspacing, ...)
    else {
        if (is.null(object$x) | is.null(object$y))
            stop ("not a valid 'traps' object")
        nd <- nrow(object)
        np <- NA
        if (detector(object) %in% .localstuff$polydetectors) {
            spacing <- NA
        }
        else {
            spacing <- spacing(object)
            if (is.factor(covariates(object))) {
                susage <- by(usage(object), covariates(object), function(y) apply(y,2,sum))
                sumusage <- matrix(unlist(susage), byrow = T, nrow = length(susage))
                dimnames(sumusage) <- list(levels(covariates(object)), names(susage[[1]]))
            }
            else if (!is.null(usage(object)))  sumusage <- apply(usage(object),2,sum)
                 else sumusage <- NULL
            sumcovar <- NULL
                tempcovar <- covariates(object)
            if (!is.null(tempcovar))
                if ((nrow(tempcovar)>0) & (ncol(tempcovar)>0)) ## amended to check ncol 2009 09 18
                    sumcovar <- summary(tempcovar, ...)
        }

        ## defaults
        area <- NA
        totallength <- NA
        np = ndetector(object)
        spacex <- NA
        spacey <- NA
        ## unblock these 2010-12-17
        sumusage <- NULL
        sumcovar <- NULL
        xrange <- range(object$x)
        yrange <- range(object$y)
        if (detector(object) %in% c('polygon', 'polygonX'))
            area <- searcharea(object)
        if (detector(object) %in% c('transect', 'transectX'))
            totallength <- transectlength(object)

        temp <- list (
          detector = detector(object),
          ndetector = nd,
          npart = np,
          xrange = xrange,
          yrange = yrange,
          spacing = spacing,
          area  = area,
          totallength = totallength,
          usage = sumusage,
          covar = sumcovar
        )

        class(temp) <- 'summary.traps'
        temp
    }
}
###############################################################################

print.summary.traps <- function (x, terse = FALSE, ...) {

#### NEED TO HANDLE CLUSTER, CLUSTERTRAP 2011-04-12

    if (!terse)
    cat ('Object class     ', 'traps', '\n')
    cat ('Detector type    ', x$detector, '\n')
    if (x$detector %in% c('polygon', 'polygonX')) {
        cat ('Number vertices  ', x$ndetector-x$npart, '\n')  ## assume each polygon closed
        cat ('Number polygons  ', x$npart, '\n')
        cat ('Total area       ', sum(x$area), 'ha \n')
    }
    else if (x$detector %in% c('transect', 'transectX')) {
        cat ('Number vertices  ', x$ndetector, '\n')
        cat ('Number transects ', x$npart, '\n')
        cat ('Total length     ', x$totallength, 'm \n')
    }
    else {
        cat ('Detector number  ', x$ndetector, '\n')
        cat ('Average spacing  ', x$spacing, 'm \n')
    }
    cat ('x-range          ', x$xrange, 'm \n')
    cat ('y-range          ', x$yrange, 'm \n')
    if (!terse) {
        if (!is.null(x$covar)) {
            cat ('\n')
            cat ('Summary of covariates', '\n')
            print(x$covar, ...)
        }
        if (!is.null(x$usage)) {
            cat ('Usage by occasion', '\n')
            print(x$usage, ...)
        }
    }
}

###############################################################################

####################################
## Class : capthist
## capture data
####################################

###############################################################################

plot.popn <- function (x, add = FALSE, frame = TRUE, circles = NULL, ...) {

    if (ms(x)) {
        ## force shared frame
        temp <- do.call(rbind, lapply(x, function(y) attr(y,'boundingbox')))
        vertices <- apply(temp,2,range)
        for (i in 1:length(x)) attr(x,'boundingbox') <- vertices
        lapply (x, plot, add, frame, circles, ...)
        invisible()
    }
    else {
        vertices <- attr(x,'boundingbox')

        if (add==FALSE)
        {
            require(MASS)
            if (frame)
                eqscplot (x$x, x$y, xlab='', ylab='', xlim=range(vertices$x),
                    ylim=range(vertices$y), type='n', axes = FALSE, ...)
            else
                eqscplot (x$x, x$y, xlab='', ylab='', type='n', axes = FALSE,
                          ...)
        }
        if (is.null(circles) | (nrow(x) == 0))    ## second condition 2011-09-14
            points (x$x, x$y, ...)
        else {
            if (length(circles) == 1)
                circles <- rep(circles, nrow(x))
            symbols (x$x, x$y, circles = circles, inches = FALSE,
                add = TRUE, ...)
        }
        if (frame) {
            if (!is.null(attr(x,'polygon')))
                polygon (attr(x,'polygon'))
            else
                polygon (vertices)
        }
    }
}

###############################################################################

rbind.popn <- function (..., renumber = TRUE) {
## combine 2 or more popn objects
## ... may be a single list object

    dots <- match.call(expand.dots = FALSE)$...
    allargs <- list(...)

    names(allargs) <- lapply(dots, as.character)

    if ((length(dots)==1) & (!inherits(allargs[[1]],'popn'))) allargs <- allargs[[1]]

    if (length(allargs)==1) return(allargs[[1]])

    ## check input
    check <- function (x) {
        if (!is(x,'popn'))
            stop ("all arguments must be 'popn' objects")
        if (is.null(covariates(x)) != is.null(covariates(allargs[[1]]) ))
            stop ("covariates must be provided for all or none")
    }
    sapply (allargs, check)

    ## row names
    an <- unlist(sapply(allargs, row.names, simplify=F))
    names(an) <- NULL
    if (any(duplicated(an))) # renumber
    {
        if (renumber) rn <- 1:length(an)
        else {
            for (i in 1:length(an)) an[[i]] <- paste(an[[i]],i,sep='.')
            rn <- unlist(an)
        }
    }
    else rn <- an

    ## construct output
    animals <- data.frame(abind (allargs, along=1), row.names = rn)
    names(animals) <- c('x','y')
    class(animals) <- c('popn', 'data.frame')
    ## following 2 lines modified 2010-06-13
    attr(animals, 'Ndist') <- 'user'
    attr(animals, 'model2D') <- attr(allargs[[1]], 'model2D')
    xl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$x))
    yl <- range(sapply(allargs, function(x) attr(x,'boundingbox')$y))
    attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
    if (!is.null(covariates(allargs[[1]]))) {
        cov <- lapply(allargs, function(x) covariates(x))
        covariates(animals) <- data.frame(abind(cov, along=1), row.names=rn)
    }
    animals
}
###############################################################################

subset.popn <- function (x, subset = NULL, sessions = NULL, poly = NULL,
    poly.habitat = TRUE, keep.poly = TRUE, renumber = FALSE, ...)
## x - popn object
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## sessions - vector (integer or logical) to subscript sessions
{
    if (ms(x)) {
        if (!is.null(sessions))
            x <- x[sessions]
        out <- vector('list')
        for (i in 1:length(x)) {
            if (is.list(subset))
                sset <- subset[[i]]
            else
                sset <- subset
            out[[i]] <- subset(x[[i]], subset[[i]], NULL, renumber, ...)
        }
        class(out) <- c('list','popn')
        out
    }
    else {
        #-------------------------
        # default subset is all
        if (is.null(subset))
            subset <- 1:nrow(x)
        #-------------------------
        # restrict to a polygon
        # added 2011-10-20

        if (!is.null(poly)) {
            OK <- pointsInPolygon(x, poly)
            if (!poly.habitat)
                OK <- !OK
            subset <- subset[OK]
        }
        #-------------------------
        # apply subsetting
        pop <- x[subset,]
        #-------------------------
        if (renumber)
            rownames(pop) <- 1 : nrow(pop)

        class(pop) <- c('popn', 'data.frame')
        attr(pop, 'Ndist') <- NULL     ## no longer known
        attr(pop, 'model2D') <- NULL   ## no longer known
        attr(pop, 'boundingbox') <- attr(x, 'boundingbox')
        if (!is.null(poly) & keep.poly) {
            attr(pop, 'polygon') <- poly
            attr(pop, 'poly.habitat') <- poly.habitat
        }
        if (!is.null(covariates(x))) {
            covariates(pop) <- covariates(x)[subset,,drop=FALSE]
        }
        pop
    }
}

###############################################################################
subset.capthist <- function (x, subset=NULL, occasions=NULL, traps=NULL,
    sessions=NULL, cutval=NULL, dropnullCH=TRUE, dropnullocc=FALSE,
    dropunused = TRUE, droplowsignals = TRUE, renumber=FALSE, ...)  {

## modified 2011-11-14; 2012-02-11

## x - capthist object (array with 2 or 3 dimensions)
## subset - vector (character, logical or integer) to subscript rows (dim(1)) of x
## occasions - vector (integer or logical) to subscript occasions
## traps - vector (character, integer or logical) to subscript rows of traps object
## sessions - vector (integer or logical) to subscript sessions

    if (ms(x)) {
        if (is.null(sessions)) sessions <- 1:length(x)
        temp <- lapply (x[sessions], subset,
            subset = subset,
            occasions = occasions,
            traps = traps,
            sessions = sessions,   ## inserted 2009 10 01
            cutval = cutval,
            dropnullCH = dropnullCH,
            dropnullocc = dropnullocc,
            dropunused = dropunused,
            droplowsignals = droplowsignals,
            renumber = renumber, ...)
        class(temp) <- c('list', 'capthist')
        if (length(temp) == 1) temp <- temp[[1]]  ## 2009 09 25
        return(temp)
    }
    else {
        detector <- detector(traps(x))
        dim3 <- length(dim(x)) == 3
        nk <- ndetector (traps(x))
        ## detectionindex is non-exported function 2012-02-11
        if (is.logical(subset) & (length(subset) != nrow(x)))
            stop ("if 'subset' is logical its length must match number of animals")
        if (is.logical(occasions) & (length(occasions) != ncol(x)))
            stop ("if 'occasions' is logical its length must match number of occasions")
        if (is.logical(traps) & (length(traps) != nk))
            stop ("if 'traps' is logical its length must match number of detectors")
        if (is.null(occasions)) occasions <- 1:ncol(x)
        if (is.null(traps))  traps <- 1:nk
        if (is.null(subset)) subset <- 1:nrow(x)

        #############################################
        ## coerce subset, traps, occasions to logical
        if (is.character(subset))
            subset <- dimnames(x)[[1]] %in% subset
        else
            if (!is.logical(subset))
                subset <- (1:nrow(x)) %in% subset
        if (is.character(traps))
            traps <- rownames(traps(x)) %in% traps
        else
            if (!is.logical(traps))
                traps <- (1:nk) %in% traps
        if (!is.logical(occasions))
            occasions <- (1:ncol(x)) %in% occasions

        #####################################
        ## preliminaries for signal detectors
        if (detector %in% c('cue','signal','signalnoise')) {
            detectionindices <- detectionindex(x)
            if (is.null(cutval)) cutval <- attr(x, 'cutval')
            if (cutval < attr(x, 'cutval'))
                stop ("cannot decrease 'cutval'")
            if (detector == 'signalnoise')
                signalOK <- (signal(x)-noise(x)) >= cutval
            else
                signalOK <- signal(x) >= cutval

            newcount <- table(
                factor(animalID(x,names=FALSE), levels=1:nrow(x))[signalOK],
                factor(occasion(x))[signalOK],
                factor(trap(x, names=FALSE), levels=1:nk)[signalOK])
            x[] <- newcount * sign(x)  ## retain deads, in principle
        }

        ###########################
        ## condition missing values
        x[is.na(x)] <- 0

        #################################
        ## optionally drop traps never used on the specified occasions
        if (dropunused && !is.null(usage(traps(x)))) {
            used <- apply(usage(traps(x))[,occasions, drop=F],1,sum) > 0
            traps <- traps & used
        }

        #############################
        ## for signalframe 2012-02-11
        if (dim3) {
            i <- x
            i[] <- 1:length(i)
            i <- i[subset, occasions, traps, drop = FALSE]
        }

        #################################
        ## prepare to drop null histories
        if (dropnullCH) {
            if (dim3) {
                nonnull <- apply(abs(x[,occasions, traps, drop=FALSE]),1,sum) > 0
            }
            else {
                x[!(abs(x) %in% (1:nk)[traps])] <- 0
                nonnull <- apply(abs(x[,occasions, drop=F]),1,sum) > 0
            }
        }
        else
            nonnull <- rep(TRUE, nrow(x))
        subset <- subset & !is.na(subset) & nonnull
        if (nrow(x)==0) subset <- 0

        #################################
        ## perform main subset operation
        if (dim3) {
            temp <- x[subset, occasions, traps, drop = FALSE]
            ## parallel operation to track indices for dim3 case
        }
        else {
            temp <- x[subset, occasions, drop = FALSE]
            ## drop rejected trap sites; 'abs' added 2011-11-14
            temp[!(abs(temp) %in% (1:nk)[traps])] <- 0
        }

        nocc <- ncol(temp)

        #################################
        ## drop null occasions
        OK2 <- rep(T,nocc)
        if (nrow(temp)>0)
       if ((!(detector %in% c('cue','signal','signalnoise'))) &&
            any( apply(abs(temp),2,sum) ==0)) {
            if (dropnullocc) {
                OK2 <- apply(abs(temp),2,sum) > 0
                if (dim3) {
                    temp <- temp[,OK2,, drop = FALSE]
                    i <- i[,OK2,, drop = FALSE]
                }
                else
                    temp <- temp[,OK2, drop = FALSE]
            }
            else
                warning ("no detections on occasion(s) ",
                    paste((1:nocc)[apply(abs(temp),2,sum) ==0], collapse=', '), "\n")
            nocc <- dim(temp)[2]  # refresh
        }
        #################################
        ## attributes
        class(temp) <- 'capthist'
        traps(temp) <- subset (traps(x), traps)
        covariates(temp) <- covariates(x)[subset,,drop = FALSE]
        usage(traps(temp)) <- usage(traps(x))[traps, occasions,
            drop = FALSE][,OK2, drop = FALSE]  ## drop null occasions
        session(temp) <- session(x)
        attr(temp, 'n.mash') <- attr(x, 'n.mash')
        attr(temp, 'centres') <- attr(x, 'centres')

        ###################################
        ## subset signal of signal capthist
        if (detector %in% c('cue','signal','signalnoise')) {
            if (!droplowsignals) {
                if (any(x<=0))
                    stop ("droplowsignals = FALSE cannot be applied to CH",
                          " objects with incomplete detection")
                signalOK <- TRUE
            }
            # otherwise signalOK remains a logical vector with length equal
            # to the original number of detections
            retained <- detectionindices %in% i
            OK <- retained & signalOK
            df <- attr(x, 'signalframe')  ## unchanged, so far
            attr(temp, 'signalframe') <- df[OK,, drop=FALSE]
            attr(temp, 'cutval') <- cutval
        }

        ############################################
        ## subset xy of polygon or transect capthist
        if (!is.null(xy(x))) {
            df <- data.frame(trap=trap(x, names = F),
                             occ=occasion(x),
                             ID=animalID(x,names = F),
                             x=xy(x)[,1], y=xy(x)[,2])
            ## subset, traps and occasions are logical vectors 2011-01-21, 2011-11-14
            OK <- occasions[df$occ] & subset[df$ID] & traps[df$trap]
            df <- df[OK,, drop = FALSE]
            df <- df[order(df$trap, df$occ, df$ID),, drop=FALSE]
            attr(temp, 'detectedXY') <- df[,c('x','y')]
        }

        ################################################
        # 2011-03-29, 2011-11-14
        # reassign trap numbers to allow for dropunused
        if (!dim3) {
            oldtrapnum <- 1:nk
            newtrapnum <- match(oldtrapnum, oldtrapnum[traps])
            trapsites <- as.numeric(abs(temp))
            temp[temp!=0] <- as.numeric(sign(temp[temp!=0])) * newtrapnum[trapsites]
        }

        ## renumber if desired
        if (renumber) {
            if (length(dim(temp))==3)
                dimnames(temp) <- list(1:nrow(temp),1:nocc,NULL)   # renew numbering
            else
                dimnames(temp) <- list(1:nrow(temp),1:nocc)
        }
        temp
    }
}

###############################################################################
##
## MS.capthist and rbind.capthist removed to rbind.capthist.R 13/9/2011
##
###############################################################################

sort.capthist <- function (x, decreasing = FALSE, by = '', byrowname = TRUE, ...) {
    if (ms(x)) {
        newx <- vector(mode='list')
        for (i in 1:length(x)) {
            xi <- subset(x, session=i)
            newx[[i]] <- sort(xi, decreasing=decreasing, by=by)
        }
        temp <- do.call(MS.capthist, newx)
        names(temp) <- names(x)
        session(temp) <- session(x)
        temp
    }
    else {
        if (is.character(by)) {
            if (by == '') {
                by <- vector(mode='list')   ## zero-length
            }
            else {
                if (!all(by %in% names(covariates(x))))
                    stop ("unrecognised sort field(s)")
                by <- as.list(covariates(x)[,by, drop=FALSE])
            }
        }
        else {
            by <- as.list(data.frame(by))
        }
        if (byrowname) by$rownames <- row.names(x)
        by$decreasing <- decreasing
        roworder <- do.call(order, by)
        newrownames <- rownames(x)[roworder]
        tempx <- x   ## for detectionorder
        tempx[tempx!=0] <- 1:length(animalID(x))
        if (length(dim(x))==2) {
            tempx[,] <- tempx[roworder,]
            x[,] <- x[roworder,]
        }
        else {
            tempx[,,] <- tempx[roworder,,]
            x[,,] <- x[roworder,,]
        }
        rownames(x) <- newrownames
        ## covariates
        if (!is.null(covariates(x))) {
            temp <-covariates(x)[roworder,,drop=F]
            rownames(temp) <- newrownames
            covariates(x) <- temp
        }

        detectionorder <- tempx[tempx!=0]
        ## signal
        sf <- attr(x, 'signalframe')
        if (!is.null(sf))
            attr(x, 'signalframe') <- sf[detectionorder,,drop=FALSE]
        ## xy
        if (!is.null(xy(x)))
            xy(x) <- xy(x)[detectionorder,,drop=FALSE]
        ## return
        x
    }
}

print.capthist <- function (x,..., condense = FALSE, sortrows = FALSE)
{
    ## recursive if list of capthist
    if (ms(x)) lapply (x, print.capthist, ..., condense = condense,
        sortrows = sortrows)
    else { # strip attributes, but why bother?
        cat('Session = ', session(x), '\n')
        if (is.null(detector(traps(x)))) {
            print.default(x, ...)
            invisible (x)
        }
        else
        if (detector(traps(x)) == 'unmarked') {
            temp <- apply(x>0, 2:3, sum)
            colnames(temp) <- rownames(traps(x))
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (detector(traps(x)) == 'presence') {
            temp <- apply(apply(x>0, 2:3, sum) > 0,2,any)*1.0
            print.default(temp, ...)
            invisible (temp)
        }
        else
        if (condense & (detector(traps(x)) %in% c('proximity', 'count',
                                                  'cue','signal','signalnoise'))) {
            temp <- apply(x, 3, function(y) y[apply(abs(y),1,sum)>0,, drop=F])
            trapnames <- rownames(traps(x))
            traps <- trapnames[rep(1:length(temp), sapply(temp, nrow))]
            detections <- as.matrix(abind(temp, along=1))
            temp <- data.frame(AnimalID = rownames(detections), Trap = traps, detections,
                stringsAsFactors = FALSE, row.names=NULL)
            names(temp)[3:ncol(temp)] <- 1:(ncol(temp)-2)
            if (sortrows) {
                lab <- temp$AnimalID
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                temp <- temp[order(lab, temp$Trap ), ]
            }
            rownames(temp) <- 1:nrow(temp)
            print(temp, ...)
            invisible (temp)
        }
        else {
            temp <- array (x, dim=dim(x))
            dimnames(temp) <- dimnames(x)
            if (sortrows) {
                lab <- dimnames(temp)[[1]]
                if (suppressWarnings( all(!is.na(as.numeric(lab)))))
                    lab <- as.numeric(lab)
                if (length(dim(x)) == 3) temp <- temp[order(lab),,]
                else temp <- temp[order(lab),]
            }
            print.default(temp,...)
            invisible(temp)
        }
    }
}
############################################################################################

plot.capthist <- function(x, rad = 5,
   hidetraps = FALSE, tracks = FALSE,
   title = TRUE, subtitle = TRUE,
   add = FALSE,
   varycol = TRUE, icolours = NULL, randcol = FALSE,
   lab1cap = FALSE, laboffset = 4,
   ncap = FALSE,
   splitocc = NULL, col2 = 'green',
   type = 'petal',
   cappar = list(cex=1.3, pch=16, col='blue'),
   trkpar = list(col='blue', lwd=1),
   labpar = list(cex=0.7, col='black'),
   ...)

    # reorganised 2010-03-30, 2011-05-09
    # see also version in d:\single sample with stems=F, mst=F 2009 02 22

{
## recursive if list of capthist
    if (ms(x)) {
        sapply (x, plot.capthist,
            rad = rad, hidetraps = hidetraps, tracks = tracks,
            title = title, subtitle = subtitle, add = add, varycol = varycol, icolours =
            icolours, randcol = randcol, lab1cap = lab1cap, laboffset =
            laboffset, ncap = ncap, splitocc = splitocc, col2 = col2,
            type = type, cappar = cappar, trkpar = trkpar, labpar = labpar, ...)
    }
    else {

        plotprox <- function (x) {
            x <- abs(x)   # do not distinguish deads for now
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day
            dx  <- rep((cos((1:nocc) * 2 * pi / nocc) * rad), ncol(x))[x>0]
            dy  <- rep((sin((1:nocc) * 2 * pi / nocc) * rad), ncol(x))[x>0]
            if (varycol) icol <<- icol+1  ## 2009 10 02
            par(trkpar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            if (tracks) lines (traps$x[x]+dx, traps$y[x]-dy)
            par(cappar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            points (traps$x[x]+dx, traps$y[x]-dy)
        }
        plotpolygoncapt <- function (xy) {
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day
            if (varycol) icol <<- icol+1  ## 2009 10 02
            par(trkpar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            if (tracks) lines (xy$x, xy$y)
            par(cappar)
            if (varycol) par(col=icol)    ## override 2009 10 02
            points (xy$x, xy$y)
        }
        plotcapt <- function (x) {
            x <- abs(x)   # do not distinguish deads for now
            dx <- (cos((1:nocc) * 2 * pi / nocc) * rad)
            dy <- (sin((1:nocc) * 2 * pi / nocc) * rad)
            occ2 <- (1:nocc) %in% splitocc
            x0 <- x>0
            x1 <- x0 & !(occ2)
            x2 <- x0 & occ2
            x[!x0] <- 1   # fool into keeping length..
            usesplit <- !is.null(splitocc) & (sum(x2)>0)
            if (varycol) icol <<- icol+1
            par(trkpar)
            if (varycol) par(col=icol)   # override

            if (tracks) {
                if (usesplit) {
                    par(col=col2)
                    lines ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])  # all
                    par(trkpar)
                    if (varycol) par(col=icol)   # override
                    lines ((traps$x[x]+dx)[x1], (traps$y[x]-dy)[x1])  # pre-split
                }
                else
                lines ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])  # all
            }

            par(cappar)
            if (varycol) par(col=icol)   # override
            points ((traps$x[x]+dx)[x0], (traps$y[x]-dy)[x0])
            if (usesplit) {
              par(col=col2)
              points ((traps$x[x]+dx)[x2], (traps$y[x]-dy)[x2])
            }
        }
        labcapt <- function (n) {
            if ( detectr %in% c('proximity', 'count', 'polygonX',
                'transectX', 'cue', 'signal', 'signalnoise', 'polygon',
                'transect', 'unmarked', 'presence') ) {
                warning ("labels not implemented for this detector type")
            }
            else {
                t1 <- abs(x[n,])
                o1 <- sum(cumsum(t1)==0)+1   # first occasion
                t1 <- t1[o1]                 # first trap site
                dx  <- (cos((1:nocc) * 2 * pi / nocc) * rad)[o1]
                dy  <- (sin((1:nocc) * 2 * pi / nocc) * rad)[o1]
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (traps$x[t1]+dx+laboffset[1], traps$y[t1]-dy+laboffsety, row.names(x)[n])
            }
        }
        labhead <- function (n, df) {
                par(labpar)
                if (varycol) par(col=n)   # override
                laboffsety <- ifelse (length(laboffset)==2, laboffset[2], laboffset[1])
                text (head(df[[n]],1)$x++laboffset[1], head(df[[n]],1)$y+laboffsety,
                      row.names(x)[n])
        }
        ncapt <- function (x) {
            ## if (detectr %in% c('proximity', 'count','cue')) {
            if (detectr %in% .localstuff$detectors3D){
               temp <- t(apply (abs(x),c(2,3),sum)) # capts/trap/day)
            }
            else if (detectr %in% .localstuff$polydetectors){
               stop ("ncap does not work with polygon and similar detectors")
            }
            else {
               fx   <- factor(x, levels=0:nrow(traps))
               temp <- table (fx, col(x)) # capts/trap/day
               temp <- temp[-1,,drop=F] # drop zeros
            }
            dx  <- rep(cos((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            dy  <- rep(sin((1:nocc) * 2 * pi / nocc) * rad, rep(nrow(traps),nocc))
            par(labpar)
            par(adj=0.5)
            OK <- temp>0
            text ((traps$x[row(temp)]+dx)[OK], (traps$y[row(temp)]-dy)[OK], as.character(temp[OK]))
        }

        plotsignal <- function (df, minsignal, maxsignal,n) {
            # df is a dataframe for one animal
            # some dupl points will be over plotted - could increase rad for
            # captures after first at a trap on a given day

            ## function (signal, occasion, trap, minsignal, maxsignal, n)
            .localstuff$i <- .localstuff$i+1
            dx <- rep( (cos((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            dy <- rep( (sin((.localstuff$i) * 2 * pi / n) * rad), nrow(df))
            sq <- order(df$signal)     # plot darkest points last
            df <- df[sq,]
            df$trap <- as.character(df$trap)
            if (maxsignal>minsignal)
                greycol <- grey(0.7 * (1 - (df$signal-minsignal)/(maxsignal-minsignal)))
            else
                greycol <- grey(0.5)
            if (tracks) lines (traps$x[df$trap]+dx, traps$y[df$trap]-dy, col=greycol)
            par(cappar)
            points (traps[df$trap,'x']+dx, traps[df$trap,'y']-dy, col = greycol)
        }

        ###########
        ## MAINLINE

        traps <- traps(x)
        detectr <- detector(traps)
        nocc <- ncol(x)
        nanimal <- nrow(x)

        if (type == 'petal')
            cappar <- replacedefaults (list(cex=1.3, pch=16, col='blue'), cappar)
        if (type %in% c('n.per.cluster','n.per.detector'))
            cappar <- replacedefaults (list(cex = 3, pch = 21), cappar)

        trkpar <- replacedefaults (list(col='blue', lwd=1), trkpar)
        labpar <- replacedefaults (list(cex=0.7, col='black'), labpar)
        initialpar <- par(cappar)

        if (!add) plot(traps, hidetr=hidetraps, ...)

        if (type == 'petal') {
            if (is.null(icolours)) icolours <- topo.colors((nanimal+1)*1.5)
            if (varycol) {
                if (randcol) icolours <- sample(icolours)
                test <- try (palette(icolours))  ## too many?
                if (inherits(test, 'try-error'))
                    stop ("requested too many colours; ",
                          "try with varycol = FALSE")
                icol <- 0
            }
            if ((nocc == 1) & ! (detectr %in% c('signal','signalnoise'))) rad <- 0

            if ( detectr %in% c('proximity', 'count', 'cue', 'unmarked', 'presence') )
            {
                w <- apply(x,1:2,function(x) (abs(x)>0) * (1:length(x)))
                w <- aperm(w, c(2,3,1))
                apply( w, 1, plotprox )
            }
            else
            if ( detectr %in% .localstuff$polydetectors ) {
                ## occasions not distinguished
                lxy <- split (xy(x), animalID(x, names=FALSE))
                lapply (lxy,plotpolygoncapt)
            }
            else
            if ( detectr %in% c('signal','signalnoise') )
            {
                .localstuff$i <- 0
                temp <- data.frame(ID = animalID(x), occ = occasion(x),
                                   trap=trap(x), signal = signal(x))
                lsignal <- split(temp, animalID(x, names = FALSE))
                lapply(lsignal, plotsignal, minsignal = min(temp$signal),
                    maxsignal = max(temp$signal), n=nanimal)
            }
            else  {   ## single, multi-catch traps
                apply( x, 1, plotcapt )
            }

            if (lab1cap) {
                if ( detectr %in% .localstuff$polydetectors ) {
                    lxy <- split (xy(x), animalID(x, names = FALSE))
                    sapply(1:nanimal, labhead, df=lxy)
                }
                else sapply(1:nanimal, labcapt)
            }

            if (ncap) { ncapt(x)}

        }
        else if (type %in% c('n.per.cluster','n.per.detector')) {
            if (type == 'n.per.detector') {
                ## never yields zeros
                temp <- table(trap(x), animalID(x))>0
                nj <- apply(temp,1,sum)
                centres <- traps(x)[names(nj),]
            }
            else if (type == 'n.per.cluster') {
                nj <- cluster.counts(x)
                centres <- cluster.centres(traps)
                # hide zeros, if present
                centres <- centres[nj>0,]
                nj <- nj[nj>0]
            }
            else stop ("unrecognised type")

            if (is.null(icolours)) {
                icolours <- topo.colors(max(nj)*1.5)
            }
            opal <- palette()
            palette (icolours)
            npal <- length(icolours)
            if (max(nj) < npal) {
                cols <- npal-nj
            }
            else {
                cols <- round(npal * (1-nj/max(nj)))
            }
            if (cappar$pch == 21)
                fg <- 'black'
            else
                fg <- cols
            par(cappar)
            points(centres, col = fg, bg = cols, pch = cappar$pch, cex = cappar$cex)
            palette(opal)

            if (ncap) {
                par(labpar)
                par(adj=0.5)
                vadj <- diff(par()$usr[3:4])/500  ## better centring!
                text(centres$x, centres$y + vadj, nj)
            }

            ## should export data for legend:
            ## count classes 0, 1-, 2-,...,-max
            ## and corresponding colours
            tempcol <- npal- (1:max(nj))
            output <- data.frame(
                legend = 1:max(nj),
                col = tempcol,
                colour = icolours[tempcol],
                stringsAsFactors = FALSE
            )
        }
        else
            stop ("type not recognised")


        ####################################################
        ## Titles
        if (is.logical(title)) {
            txt <- ifelse (is.null(session(x)), paste(deparse(substitute(x)),
                       collapse=''), session(x))
            title <- ifelse(title, txt, '')
        }
        if (title != '') {
            par(col='black')
            mtext(side=3,line=1.2, text = title, cex=0.7)
        }
        if (is.logical(subtitle)) {
            if(detectr %in% .localstuff$exclusivedetectors) nd <- sum(abs(x)>0)
            else nd <- sum(abs(x))
            if (subtitle) {
                if (detectr == 'cue')
                    subtitle <- paste(
                        nocc, 'occasion,' ,
                        nd, 'detections,',
                        nanimal, 'cues',
                        length(levels(covariates(x)$animal)), 'animals')
                else
                    subtitle <- paste(
                        nocc, 'occasions,' ,
                        nd, 'detections,',
                        nanimal, 'animals')
            }
            else subtitle <- ''
        }
        if (subtitle != '') {
            par(col='black')
            mtext(text = subtitle, side=3,line=0.2, cex=0.7)
        }
        ####################################################

        par(initialpar)   # restore
        if (type %in% c('n.per.detector','n.per.cluster'))
            invisible(output)
        else
            invisible(sum(abs(x)>0))
    }
}
############################################################################################

summary.capthist <- function(object, terse = FALSE, ...) {

    ## recursive if list of capthist
    if (ms(object)) {
        if (terse) {
            n     <- sapply(object, nrow)        # number caught
            nocc  <- sapply(object, ncol)        # number occasions
            ncapt <- sapply(object, function (xx) sum(abs(xx)>0))
            ndet  <- sapply(traps(object), ndetector) # number traps
            temp  <- as.data.frame(rbind(nocc, ncapt, n, ndet))
            names(temp) <- names(object)
            rownames(temp) <- c('Occasions','Detections','Animals','Detectors')
            temp
        }
        else
            lapply (object, summary.capthist, ...)
    }
    else {

        traps <- traps(object)
        detector <- detector(traps)
        cutval <- attr(object, 'cutval')   # signal strength only

        nd <- length(traps$x)

        # ni, ui, fi, M, losses etc.
        nocc <- ncol(object)
        counts <- matrix(nrow = 8, ncol = nocc)
        counts[,] <- 0
        signalsummary <- NULL
        if (nrow(object) > 0) {
            if (length(dim(object)) > 2) {
                tempx <- apply( object[,,,drop=F], c(1,2), function(x) sum(abs(x))>0)
                if (nocc>1) {  # distinction may not be needed...
                    counts [1,] <- apply(tempx, 2, function(x) sum(abs(x)>0) )
                    tempx2 <- apply(tempx, 1, function(x) cumsum(abs(x))>0)
                    counts [4,] <- apply(tempx2,1,sum)
                    counts [2,] <- c(counts[4,1],diff(counts[4,]))
                    counts [3,] <- tabulate(apply(tempx,1, function(x) sum(abs(x)>0)),nbins = nocc)
                    counts [5,] <- apply(tempx,2, function(x) sum(x<0))
                }
                else {
                    counts [1,1] <- sum(abs(tempx)>0)
                    counts [2,1] <- counts [1,1]
                    counts [3,1] <- counts [1,1]
                    counts [4,1] <- counts [1,1]
                    counts [5,1] <- sum(tempx<0)
                }
                if (detector %in% c('proximity','signal','signalnoise','cue'))
                   counts [6,] <- apply(object[,,, drop=F], 2, function(x) sum(abs(x)>0))
                ## abs (x) added following 2 statements 2011-02-09
                if (detector %in% c('count', 'polygon','transect','unmarked','presence'))
                   counts [6,] <- apply(object[,,, drop=F], 2, function(x) sum(abs(x)))
                tempt <- apply(object[,,,drop=F],c(2,3), function(x) sum(abs(x))>0)
                counts [7,] <- apply(tempt,1,sum)
            }
            else  {     ## .localstuff$exclusivedetectors
                if (nocc == 1) {
                    counts [1,1] <- sum(abs(object)>0)
                    counts [2,1] <- counts [1,1]
                    counts [3,1] <- counts [1,1]
                    counts [4,1] <- counts [1,1]
                    counts [5,1] <- sum(object<0)
                }
                else {
                    counts [1,] <- apply(object[,,drop=F], 2, function(x) sum(abs(x)>0))
                    if (nrow(object) > 0) {
                        tempM <- apply(object[,,drop=F], 1, function(x) cumsum(abs(x))>0)
                        counts [4,] <- apply(tempM,1,sum)
                    }
                    counts [2,] <- c(counts[4,1], diff(counts[4,]))
                    counts [3,] <- tabulate(apply(object[,,drop=F],1, function(x) sum(abs(x)>0)),
                                            nbins = nocc)
                    counts [5,] <- apply(object[,,drop=F],2, function(x) sum(x<0))
                }
                counts [6,] <- apply(object[,,drop=F],2, function(x) sum(abs(x)>0))
                counts [7,] <- apply(object[,,drop=F],2, function(x) length(unique(x[x!=0])))
            }
        }
        if (!is.null(traps)) {
            if (!is.null(usage(traps)))
                counts[8,] <- apply(usage(traps),2,function(x) sum(x>0))
            else
                counts[8,] <- rep(ndetector(traps),nocc)
        }

        counts <- as.data.frame(counts)
        dimnames(counts) <- list(c('n','u','f','M(t+1)','losses','detections',
                                   'detectors visited','detectors used'), 1:nocc)
        counts$Total <- apply(counts, 1, sum)
        counts$Total[4] <- counts[4, nocc]
        PSV <-  NULL
        dbar <- NULL
        if (is.null(traps)) {
            trapsum <- NULL
            signalsummary <- NULL
        }
        else {

            if (detector(traps) %in% .localstuff$individualdetectors) {
                if (-diff(counts$Total[1:2]) > 1)
                    PSV <- RPSV(object)
                if (length(dim(object)) == 2)
                    dbar <- dbar(object)
            }
            trapsum <- summary(traps)
            if (detector == 'signal')
                signalsummary <- summary(signal(object))
            if (detector == 'signalnoise')
                signalsummary <- list(signal = summary(signal(object)),
                                      noise = summary(noise(object)),
                                      diffSN = summary(signal(object)-noise(object)))
        }

        temp <- list (
            detector = detector,
            ndetector = nd,
            trapsum = trapsum,
            counts = counts,
            dbar = dbar,
            RPSV = PSV,
            cutval = cutval,        # signal, signalnoise only
            signalsummary = signalsummary
        )
        class(temp) <- 'summary.capthist'
        temp
    }
}
############################################################################################

counts <- function (CHlist, counts = 'M(t+1)') {
    if (!inherits(CHlist, 'capthist'))
        stop ("require capthist object")
    getc <- function (cnt) {
        getcnt <- function(x, maxocc) {
            temp <- x$counts[cnt,]
            lt <- length(temp)
            matrix(c(temp[-lt], rep(NA, maxocc-lt+1), temp[lt]), nrow = 1)
        }
        if (!is.list(CHlist))
            summary(CHlist)$counts[cnt,]
        else {
            maxocc <- max(sapply(CHlist,ncol))
            abind(lapply(summary(CHlist), getcnt, maxocc), along=1,
                new.names=list(session(CHlist), c(1:maxocc, 'Total')))
        }
    }
    temp <- lapply (counts, getc)
    names(temp) <- counts
    temp
}

print.summary.capthist <- function (x, ...) {
    cat ('Object class     ', 'capthist', '\n')
    print(x$trapsum, terse=TRUE)
    cat ('Counts by occasion \n')
    print(x$counts, ...)
    cat ('\n')
    if (x$detector %in% c('signal', 'signalnoise')) {
        cat ('Signal threshold ', x$cutval, '\n')
        print (x$signalsummary)
    }
}
############################################################################################

###############################
## Class : mask
## defines a habitat mask
###############################

subset.mask <- function (x, subset, ...) {

    if (ms(x))
        stop ("subset of multi-session mask not implemented")

    # subset may be numeric index or logical
    temp <- x[subset,,drop=F]
    spacing <- attr(x,'spacing')
    attr(temp,'type')        <- 'subset'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- attr(x,'area')
    attr(temp,'spacing')     <- spacing
    if (!is.null(covariates(x))) covariates(temp) <- covariates(x)[subset,,drop=F]
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)
    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    if (!is.null(attr(temp,'OK')))
        attr(temp,'OK') <- attr(temp,'OK')[subset]
    class(temp) <- c('mask', 'data.frame')
    temp
}
############################################################################################

rbind.mask <- function (...) {
# combine 2 or more mask objects

##    no check for multi-session masks at present
##        stop ('rbind of multi-session mask not implemented')
##
    dropduplicates <- TRUE   ## always
    allargs <- list(...)
    spacing <- attr(allargs[[1]],'spacing')
    area    <- attr(allargs[[1]],'area')
    check <- function (x) {
        if (!is(x,'mask'))
            stop ("arguments must be mask objects")
        if (attr(x,'spacing') != spacing)
            stop ("arguments must have same 'spacing' attribute")
        if (attr(x,'area') != area)
            stop ("arguments must have same area attribute")
    }
    sapply (allargs, check)
    temp <- rbind.data.frame(...)
    class(temp) <- c('mask', 'data.frame')
    tempcov <- lapply(allargs, covariates)
    covariates(temp) <- do.call(rbind, tempcov)  ## pass list of dataframes

    if (dropduplicates) {
        dupl <- duplicated(temp)
        droppedrows <- sum(dupl)
        if (droppedrows>0) {
            covariates(temp) <- covariates(temp)[!dupl,]
            temp <- temp[!dupl,]
            warning (droppedrows, " duplicate points dropped from mask")
        }
    }

    attr(temp,'type')        <- 'rbind'
    attr(temp,'meanSD')      <- getMeanSD(temp)
    attr(temp,'area')        <- area
    attr(temp,'spacing')     <- spacing
    xl <- range(temp$x) + spacing/2 * c(-1,1)
    yl <- range(temp$y) + spacing/2 * c(-1,1)

    ##  xl <- range(temp[,1])
    ##  yl <- range(temp[,2])

    attr(temp,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    temp
}
############################################################################################


read.mask <- function (file = NULL, data = NULL, spacing = NULL, columns = NULL, ...)
## 2010-04-11 SS for 'state-space' from SPACECAP
## 2011-05-11 data argument added
## 2011-11-01 re-modelled to add covariates
{
    if (is.null(data) & !is.null(file)) {
        fl <- nchar(file)
        SS <- tolower(substring(file,fl-3,fl)) == '.csv'
        if (SS) {
            data <- read.csv (file)
            if ('HABITAT' %in% names(data))
                data <- data[data$HABITAT == 1,]
        }
        else
            data <- read.table (file, ...)
    }
    else if (is.null(data))
       stop("require one of 'file' or 'data'")
    if (length(dim(data))!=2)
        stop ("require dataframe or matrix for 'data' input to read.mask")

    coln <- colnames(data)
    ixy <- match(c('x','y'), coln)
    if (any(is.na(ixy))) ixy <- 1:2
    mask <- as.data.frame(data[,ixy])
    names(mask) <- c('x', 'y')
    if (any(!apply(mask, 2, is.numeric)))
        stop ("non-numeric x or y coordinates")
    if (any(is.na(mask)))
        stop ("missing value(s) in x or y coordinates")

    class(mask) <- c('mask', 'data.frame')

    ## add covariates
    if (ncol(data) > 2) {
        df <- as.data.frame(data[,-ixy, drop = FALSE])
        if (!is.null(columns)) {
            if (!all(columns %in% names(mask)))
                stop ("columns missing from input")
            df <- df[,columns, drop=FALSE]
        }
        if (ncol(df)>0)
            covariates(mask) <- df
    }

    if (is.null(spacing))
    {
      sp      <- as.matrix(dist(as.matrix(mask)))
      spacing <- apply(sp,1,function(x) min(x[x>0]))
      spacing <- mean (spacing, na.rm=T)
    }

    area    <- spacing^2 / 10000

    xl <- range(mask$x) + spacing/2 * c(-1,1)
    yl <- range(mask$y) + spacing/2 * c(-1,1)

    attr(mask,'type')    <- 'user'
    attr(mask,'meanSD')  <- getMeanSD(mask)
    attr(mask,'area')    <- area
    attr(mask,'spacing') <- spacing
    attr(mask,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
    attr(mask,'polygon') <- NULL
    mask
}
###############################################################################

plot.mask <- function(x, border = 20, add = FALSE, covariate = NULL,
      axes = FALSE, dots = TRUE, col='grey', breaks = 12, meshcol = NA,
      ppoly = TRUE, polycol='red', ...)
{
    if (ms(x)) {
        lapply (x, plot.mask, ...)
    }
    else {

        buff <- c(-border,+border)
        if (!add) {
            require(MASS)
            eqscplot (x$x, x$y,
            xlim=range(x$x)+buff, ylim=range(x$y)+buff,
            xlab='', ylab='',
            axes=axes, type='n', ...)
        }

        if (!is.null(attr(x,'polygon')) & ppoly) {
            poly <- attr(x,'polygon')
            if (class(poly) == "SpatialPolygonsDataFrame") {
                if (!require(sp))
                    stop ("package 'sp' required to plot polygon in plot.mask")
# plot(poly, col = polycol, add = TRUE)
# poor control of colours
                plot(poly, add = TRUE)
            }
            else
                polygon (poly, col = polycol, density = 0)
        }

        if (is.null(covariate))
            covfactor <- factor(1)
        else {
            if (is.factor(covariates(x)[,covariate]))
                covfactor <- covariates(x)[,covariate]
            else
                covfactor <- cut ( covariates(x)[,covariate], breaks = breaks)
        }
        ncolour <- length(levels(covfactor))
        if (length(col) < ncolour)
            col <- heat.colors(ncolour)   # default set
        cols <- col[as.numeric(covfactor)]

        if (dots) {
            points (x$x, x$y, col = cols, pch = 16, cex = 0.4)
        }
        else {
            pixelsize <- attr(x,'spacing')
            dx <- c(-0.5, -0.5, +0.5, +0.5) * pixelsize
            dy <- c(-0.5, +0.5, +0.5, -0.5) * pixelsize
            plotpixel <- function (xy) {
                polygon (xy[1]+dx, xy[2]+dy, col=col[xy[3]],density=-1, border = meshcol)
            }
            apply(cbind(x,as.numeric(covfactor)),1,plotpixel)
        }

        if (!is.null(covariate))
            invisible(levels(covfactor))
    }
}
###############################################################################

summary.mask <- function(object, ...) {

  if (ms(object)) {
      temp <- lapply(object, summary.mask)
      class(temp) <- c('summary.mask', 'list')
      temp
  }
  else {
      if (is.null(object$x) | is.null(object$y))
          stop ("not a valid mask")
      nd <- length(object$x)
      if (length(object$x) != length(object$y))
          stop  ("not a valid mask")

      if (!is.null(covariates(object))) {
          sumcovar <- summary(covariates(object), ...)
      } else sumcovar <- NULL

      temp <- list (
        detector = attr(object,'detector'),
        type = attr(object,'type'),
        nmaskpoints = nrow(object),
        xrange = range(object$x),
        yrange = range(object$y),
        meanSD = attr(object,'meanSD'),
        spacing = attr(object,'spacing'),
        cellarea = attr(object,'area'),
        boundingbox = attr(object,'boundingbox'),
        covar = sumcovar
      )
      class(temp) <- 'summary.mask'
      temp
  }

}
############################################################################################

print.summary.mask <- function (x, ...) {
    if (ms(x)) {
        lapply (x, print.summary.mask)
    }
    else {
      cat ('Object class     ', 'mask', '\n')
      cat ('Mask type        ', x$type, '\n')
      cat ('Number of points ', x$nmaskpoints, '\n')
      cat ('Spacing m        ', x$spacing, '\n')
      cat ('Cell area ha     ', x$cellarea, '\n')
      cat ('Total area ha    ', x$cellarea * x$nmaskpoints, '\n')
      cat ('x-range m        ', x$xrange, '\n')
      cat ('y-range m        ', x$yrange, '\n')
      cat ('Bounding box     ','\n')
      print (x$boundingbox, ...)
      cat ('\n')
      if (!is.null(x$covar)) {
          cat ('Summary of covariates', '\n')
          print(x$covar, ...)
      }
  }
}
############################################################################################

####################################################
## Class : secr
## spatially explicit capture-recapture model fit
####################################################

############################################################################################

trim.secr <- function (object, drop = c('mask','design','design0'), keep = NULL) {
    trim.default(object, drop = drop, keep = keep)
}
############################################################################################


predict.secr <- function (object, newdata = NULL, se.fit = TRUE, alpha = 0.05,
    savenew = FALSE, scaled = FALSE, ...) {

    if (is.null(object$fit)) {
        warning ("empty (NULL) object")
        return(NULL)
    }
    if (is.null(newdata)) newdata <- secr.make.newdata (object)

    parindices <- object$parindx
    models <- object$model
    if ('cuerate' %in% object$realnames) {
        parindices$cuerate <- max(unlist(parindices)) + 1
        models$cuerate <- ~1
    }
    if (object$detectfn %in% c(12,13)) {
        ## experimental parameters not fitted
        ## construct dummies
        parindices$muN <- max(unlist(parindices)) + 1
        parindices$sdN <- max(unlist(parindices)) + 1
        models$muN <- ~1
        models$sdN <- ~1
        object$link$muN <- 'identity'
        object$link$sdN <- 'identity'
    }

    ## allow for fixed beta parameters 2009 10 19
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(NA, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
        fb[is.na(fb)] <- object$fit$par
        beta <- fb    ## complete
    }
    else {
        beta <- object$fit$par
        beta.vcv <- object$beta.vcv
    }

    getfield <- function (x) {
        if ((x == 'D') & userD(object)) {
            ## user-supplied density function
            ## return only intercept
            lpred <- matrix(ncol = 2, nrow = nrow(newdata),
               dimnames=list(NULL,c('estimate','se')))
            D0 <- parindices[[x]][1]
            lpred[,1] <- beta[D0]
            lpred[,2] <- beta.vcv[D0,D0]^0.5
            return(lpred)
        }
        else {
            secr.lpredictor (newdata = newdata, model = models[[x]],
                indx = parindices[[x]], beta = beta, field = x,
                beta.vcv = beta.vcv)
        }
    }

    predict <- sapply (object$realnames, getfield, simplify = FALSE)

    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!
    if (se.fit)  out <- list(nrow(newdata))
    else {
        out <- newdata
        ## add columns for real parameter estimates
        for (varname in object$realnames)
            out[,varname] <- rep(NA,nrow(out))
    }
    if (!is.null(predict$pmix)) {
        nmix <- object$details$nmix

        ######################
        ## replaced 2010 03 10
        # predict$pmix$estimate <- logit(mlogit.untransform(predict$pmix$estimate, 1:nmix))

        # assuming mixture is always last dimension...
        temp <- matrix(predict$pmix$estimate, ncol = nmix)
        temp2 <- apply(temp, 1, function(est) logit(mlogit.untransform(est, 1:nmix)))
        predict$pmix$estimate <- as.numeric(t(temp2))
        ######################

        predict$pmix$se <- NA    ## uncertain
    }
    for (new in 1:nrow(newdata)) {
        lpred  <- sapply (predict, function(x) x[new,'estimate'])
        Xlpred <- Xuntransform(lpred, object$link, object$realnames)

        if (ms(object)) {
            if (!('session' %in% names(newdata))) {
                n.mash <- NULL
                ## before 2012-05-14 was:
                ## stop ("could not discern session in predict... ",
                ##       "see the package author!")
            }
            else {
                sess <- newdata[new, 'session']
                n.mash <- attr (object$capthist[[sess]], 'n.mash')
            }
            n.clust <- length(n.mash)
            if (new==1)
                oldnclust <- n.clust
            else if (n.clust != oldnclust)
                warning ("number of mashed clusters varies between sessions")
        }
        else {
            n.mash <- attr (object$capthist, 'n.mash')
            n.clust <- length(n.mash)
        }

        #####################
        ## 2010 02 14 rescale
        if (object$details$scalesigma & scaled) {
            Xlpred['sigma'] <- Xlpred['sigma'] / Xlpred['D']^0.5
            lpred['sigma'] <- NA   ## disable further operations for SE
        }
        if (object$details$scaleg0 & scaled) {
            Xlpred['g0'] <- Xlpred['g0'] / Xlpred['sigma']^2
            lpred['g0'] <- NA   ## disable further operations for SE
        }
        #####################

        if (se.fit) {
            selpred <- sapply (predict,function(x) x[new,'se'])
            temp <- data.frame (
              row.names = object$realnames,
              link = unlist(object$link[object$realnames]),
              estimate = Xlpred,
              SE.estimate = se.Xuntransform (lpred, selpred, object$link, object$realnames),
              lcl = Xuntransform(lpred-z*selpred, object$link, object$realnames),
              ucl = Xuntransform(lpred+z*selpred, object$link, object$realnames)
              )
            # truncate density at zero; adjust for mash()
            if ('D' %in% row.names(temp)) {
                temp['D', -1][temp['D',-1]<0] <- 0
                if (!is.null(n.mash)) {
                    temp['D', -1] <- temp['D', -1] / n.clust
                }
            }

            if (nrow(newdata)==1) out <- temp
            else {
                out[[new]] <- temp
                names(out)[new] <- paste (
                        paste(names(newdata),'=', unlist(lapply(newdata[new,],as.character)),
                        sep=' ',collapse=', '),
                    sep=',')
            }
        }
        else { # no SE; terse format
            if ('D' %in% names(Xlpred)) {
                Xlpred['D'] <- ifelse (Xlpred['D']<0, 0, Xlpred['D'])
                if (!is.null(n.mash)) {
                    Xlpred['D'] <- Xlpred['D'] / n.clust
                }
            }
            out[new, (ncol(newdata)+1) : ncol(out)] <- Xlpred
        }
    }
    if (savenew) attr(out, 'newdata') <- newdata
    out
}
############################################################################################

############################################################################################
## 2010-10-22

predict.secrlist <- function (object, newdata = NULL, se.fit = TRUE, alpha = 0.05,
    savenew = FALSE, scaled = FALSE, ...) {
    lapply(object, predict, newdata, se.fit, alpha, savenew, scaled, ...)
}

############################################################################################

coef.secr <- function (object, alpha=0.05, ...) {
    beta   <- object$fit$par
    if (!is.null(object$beta.vcv))
        sebeta <- suppressWarnings(sqrt(diag(object$beta.vcv)))
    else sebeta <- rep(NA, length(beta))
    z <- abs(qnorm(1-alpha/2))
    temp <- data.frame(
        row.names = object$betanames,
        beta    = beta,
        SE.beta = sebeta,
        lcl = beta - z*sebeta,
        ucl = beta + z*sebeta
        )
    attr(temp, 'alpha') <- alpha
    temp
}
############################################################################################

detectpar <- function(object, ...) {
    extractpar <- function (temp) {

        if (!is.data.frame(temp))   ## assume list
            temp <- temp[[1]]
        if (!is.data.frame(temp) |
            (nrow(temp) > length(object$link)))
            stop ("unexpected input to detectpar()")

        temp <- temp[, 'estimate', drop = F]
        temp <- split(temp[,1], rownames(temp))
        temp <- c(temp, object$fixed)
        temp <- temp[parnames(object$detectfn)]
        if (object$detectfn > 9)
            temp <- c(temp, list(cutval = object$details$cutval))
        temp
    }
    if (!inherits(object,'secr'))
        stop ("requires 'secr' object")
    temppred <- predict (object, ...)
    if (ms(object)) {
        temp <- lapply(temppred, extractpar)
        names(temp) <- session(object$capthist)
        temp
    }
    else
        extractpar(temppred)
}
############################################################################################

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

print.secr <- function (x, newdata = NULL, alpha = 0.05, deriv = FALSE, ...) {

    cat ('\n')

    if (is.character(x$call))     ## pre secr 1.5 object
        cl <- x$call
    else {
        cl <- paste(names(x$call)[-1],x$call[-1], sep=' = ', collapse=', ' )
        cl <- paste('secr.fit(', cl, ')')
    }

    cat(strwrap(cl, getOption('width')), sep='\n  ')

    if (!is.null(x$version)) {
        cat ('secr ', x$version, ', ', x$starttime, '\n', sep='')
    }
    else {   ## for backward compatibility
        cat (x$fitted,'\n')
    }
    cat ('\n')

    print(summary(traps(x$capthist)), terse=TRUE)


    ###################
    ## Data description

    if (ms(x$capthist)) {
        print (summary(x$capthist, terse = TRUE))
        q <- sapply(x$capthist, function(y) attr(y,'q'))
        det <- detector(traps(x$capthist)[[1]])
    }
    else {
        det <- detector(traps(x$capthist))
        n  <- nrow(x$capthist)     # number caught
        if (length(dim(x$capthist))>2)
            ncapt <- sum(abs(x$capthist))
        else
            ncapt <- sum(abs(x$capthist)>0)

        q <- attr(x$capthist, 'q')

        if ('g' %in% x$vars) {
            Groups  <- table(group.factor(x$capthist, x$groups))
            temp <- paste (names(Groups), Groups, collapse=', ', sep='=')
            temp <- paste('(',temp,')', sep='')
        }
        else temp <- ''

        cat ('N animals       : ', n, temp, '\n')
        cat ('N detections    : ', ncapt, '\n')
        if (!is.null(q)) {
            cat ('N marking  occn : ', q, '\n')
            cat ('N sighting occn : ', ncol(x$capthist)-q, '\n')
        }
        else
            cat ('N occasions     : ', ncol(x$capthist), '\n')
        ## cat ('N detectors     : ', ndetector(traps(x$capthist)), '\n')
    }   # end of single-session

    if (det %in% .localstuff$countdetectors) {
        cat ('Count model     :  ')
        if (x$details$binomN == 0) cat ('Poisson \n')
        if (x$details$binomN == 1) cat ('Bernoulli \n')
        if (x$details$binomN < 0) cat ('Negative binomial k = ', abs(x$details$binomN), '\n')
        if (x$details$binomN > 1) cat('Binomial', x$details$binomN, '\n')
    }

    if (!ms(x$capthist))
    cat ('Mask area       : ', attr(x$mask,'area') * nrow(x$mask), 'ha \n')

    ####################
    ## Model description

    Npar <- max(unlist(x$parindx))
    ## allow for fixed beta parameters 2009 10 19
    if (!is.null(x$details$fixedbeta))
        Npar <- Npar - sum(!is.na(x$details$fixedbeta))
    AICval <- 2*(x$fit$value + Npar)
    n <- ifelse (ms(x$capthist), sum(sapply(x$capthist, nrow)), nrow(x$capthist))
    AICcval <- ifelse ((n - Npar - 1) > 0,
        2*(x$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)
    cat ('\n')
    cat ('Model           : ', model.string(x$model, x$details$userDfn), '\n')
    if (detector(traps(x$capthist))=='multi') {
        if (x$details$param == 1)
            cat ('Gardner, Royle & Wegan parameterisation for multi-catch traps','\n')
    }
    cat ('Fixed (real)    : ', fixed.string(x$fixed), '\n')
    cat ('Detection fn    : ', detectionfunctionname(x$detectfn), '\n')
    if (!x$CL)
    cat ('Distribution    : ', x$details$distribution, '\n')

    cat ('N parameters    : ', Npar, '\n')
    cat ('Log likelihood  : ', -x$fit$value, '\n')
    cat ('AIC             : ', AICval, '\n')
    cat ('AICc            : ', AICcval, '\n')

    cat ('\n')
    cat ('Beta parameters (coefficients)', '\n')

    print(coef(x), ...)

    if (!is.null(x$fit$hessian)) {
      cat ('\n')
      cat ('Variance-covariance matrix of beta parameters', '\n')
      print (x$beta.vcv, ...)
    }

    # scale newdata covariates... NOT FINISHED 10 05 08
    meanSD <- attr(x$mask,'meanSD')
    if (!is.null(newdata)) {
         for (i in 1:length(newdata)) {
           ind <- match (names(newdata[i]),names(meanSD))
           if (ind>0 & !is.na(meanSD[1,ind]))
             newdata[[i]] <- (newdata[[i]] - meanSD[1,ind]) / meanSD[2,ind]
         }
     }

    cat ('\n')
    cat ('Fitted (real) parameters evaluated at base levels of covariates', '\n')

    if (!is.null(x$realpar))
        print( x$realpar )
    else {
        temp <- predict (x, newdata, alpha)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n', names(temp)[new],'\n')
                print(temp[[new]], ...)
            }
    }

    #################################
    # Derived parameters (CL)
    #################################
    if (x$CL & deriv) {

        cat ('\n')
        cat ('Derived parameters (CL only)', '\n')

        temp <- derived(x, alpha=alpha, se.esa = TRUE)
        nd <- length(temp)
        if (is.data.frame(temp)) print(temp, ...)
        else for (new in 1:nd) {
                cat('\n',names(temp)[new],'\n')
                print(temp[[new]], ...)
            }

    }
    cat ('\n')
}
############################################################################################

oneline.secr <- function (secr) {

    if (ms(secr$capthist)) {
        n <- sum(sapply (secr$capthist, nrow))
        ncapt <- sum(sapply (secr$capthist, function (x) sum(abs(x>0))))
    }
    else  {
        n  <- nrow(secr$capthist)     # number caught
        ncapt <- sum( abs( secr$capthist)>0)
    }

    Npar <- max(unlist(secr$parindx))
    ## allow for fixed beta parameters 2009 10 19
    if (!is.null(secr$details$fixedbeta))
        Npar <- Npar - sum(!is.na(secr$details$fixedbeta))

    AICval <- 2*(secr$fit$value + Npar)
    AICcval <- ifelse ((n - Npar - 1)>0,
        2*(secr$fit$value + Npar) + 2 * Npar * (Npar+1) / (n - Npar - 1),
        NA)

    c (
       model  = model.string(secr$model, secr$details$userDfn),
       detectfn = detectionfunctionname(secr$detectfn),
       npar   = Npar,
       logLik = -secr$fit$value,
       AIC    = round(AICval, 3),
       AICc   = round(AICcval, 3),
       fitted = secr$fitted
    )
}
############################################################################################

logLik.secr <- function(object, ...) {
    npar <- length(object$fit$par)
    structure (-object$fit$value, df = npar, class = 'logLik')
}

############################################################################################

AIC.secr <- function (object, ..., sort = TRUE, k = 2, dmax = 10) {

    if (k != 2)
        stop ("'AIC.secr' defined only for k = 2")

    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
          as.character(match.call(expand.dots=FALSE)$...) ))
    if (any(sapply(allargs,class) != 'secr'))
        stop ("arguments must be secr objects")
    allargs <- c(list(object), allargs)

#    if (!all.equal(sapply (allargs, function(x) detector(traps(x$capthist)))))
#    stop('Models not compatible')
#    if (!all.equal(sapply (allargs, function(x) x$CL)))
#    stop('Models not compatible')

    output <- data.frame(t(sapply(allargs, oneline.secr)), stringsAsFactors=F)
    for (i in 3:6)
    output[,i] <- as.numeric(output[,i])
    output$dAICc <- output$AICc - min(output$AICc)

    OK <- abs(output$dAICc) < abs(dmax)
    sumdAICc <- sum(exp(-output$dAICc[OK]/2))
    output$AICwt <- ifelse ( OK, round(exp(-output$dAICc/2) / sumdAICc,4), 0)

    row.names(output) <- modelnames
    if (sort) output <- output [order(output$AICc),]

    if (nrow(output)==1) { output$dAICc <- NULL; output$AICwt <- NULL}
    output
}
############################################################################################
############################################################################################

    AIC.secrlist <- function (object, ..., sort = TRUE, k = 2, dmax = 10) {

    if (k != 2)
        stop ("AIC.secr defined only for k = 2")

    if (length(list(...)) > 0)
        warning ("... argument ignored in 'AIC.secrlist'")

    modelnames <- names(object)
    allargs <- object
    if (any(sapply(allargs,class) != 'secr'))
        stop ("components of 'object' must be 'secr' objects")

#    if (!all.equal(sapply (allargs, function(x) detector(traps(x$capthist)))))
#    stop('Models not compatible')
#    if (!all.equal(sapply (allargs, function(x) x$CL)))
#    stop('Models not compatible')

    output <- data.frame(t(sapply(allargs, oneline.secr)), stringsAsFactors=F)
    for (i in 3:6)
    output[,i] <- as.numeric(output[,i])
    output$dAICc <- output$AICc - min(output$AICc)

    OK <- abs(output$dAICc) < abs(dmax)
    sumdAICc <- sum(exp(-output$dAICc[OK]/2))
    output$AICwt <- ifelse ( OK, round(exp(-output$dAICc/2) / sumdAICc,4), 0)

    row.names(output) <- modelnames
    if (sort) output <- output [order(output$AICc),]

    if (nrow(output)==1) { output$dAICc <- NULL; output$AICwt <- NULL}
    output
}
############################################################################################

vcov.secr <- function (object, realnames = NULL, newdata = NULL, byrow = FALSE, ...) {
## return either the beta-parameter variance-covariance matrix
## or vcv each real parameters between points given by newdata (byrow = TRUE)
## or vcv for real parameters at points given by newdata (byrow = TRUE)

    if (is.null(dimnames(object$beta.vcv)))
        dimnames(object$beta.vcv) <- list(object$betanames, object$betanames)

    if (is.null(realnames))
        ## average beta parameters
        return( object$beta.vcv )
    else {
        ## average real parameters
        ## vcv among multiple rows

        if (byrow) {
            ## need delta-method variance of reals given object$beta.vcv & newdata
            if (is.null(newdata))
                newdata <- secr.make.newdata (object)
            nreal <- length(realnames)
            nbeta <- length(object$fit$par)

            rowi <- function (newdatai) {
                reali <- function (beta, rn) {
                    ## real from all beta pars eval at newdata[i,]
                    par.rn <- object$parindx[[rn]]
                    mat <- model.matrix(object$model[[rn]], data=newdatai)
                    lp <- mat %*% matrix(beta[par.rn], ncol = 1)
                    untransform (lp, object$link[[rn]])
                }
                grad <- matrix(nrow = nreal, ncol = nbeta)
                dimnames(grad) <- list(realnames, object$betanames)
                for (rn in realnames)
                    grad[rn,] <- fdHess (pars = object$fit$par, fun = reali, rn = rn)$gradient
                vcv <- grad %*% object$beta.vcv %*% t(grad)
                vcv
            }

            vcvlist <- list(nrow(newdata))
            for (i in 1:nrow(newdata)) vcvlist[[i]] <- rowi(newdata[i,])
            if (length(vcvlist) == 1) vcvlist <- vcvlist[[1]]
            return(vcvlist)
        }
        else {
            newdata <- as.data.frame(newdata)
            rownames <- apply(newdata, 1, function(x) paste(names(newdata), '=', x, sep='',
                                                            collapse=','))
            vcvlist <- list()
            for (rn in realnames) {
                par.rn <- object$parindx[[rn]]
                mat <- model.matrix(object$model[[rn]], data=newdata)
                lp <- mat %*% matrix(object$fit$par[par.rn], ncol = 1)
                real <- untransform (lp, object$link[[rn]])
                real <- as.vector(real)
                ## from Jeff Laake's 'compute.real' in RMark...
                deriv.real <- switch(object$link[[rn]],
                    logit = mat * real * (1-real),
                    log = mat * real,
                    identity = mat,
                    sin = mat * cos(asin(2*real-1))/2)
                vcvlist[[rn]] <- deriv.real %*% object$beta.vcv[par.rn, par.rn] %*% t(deriv.real)
                dimnames(vcvlist[[rn]]) <- list(rownames, rownames)
            }
            names (vcvlist) <- realnames
            return (vcvlist)
        }
        ## DIFFERENT VARIANCE TO secr.lpredictor for sigma because there use se.Xuntransfom
    }
}
############################################################################################

