############################################################################################
## package 'secr'
## reduce.capthist.R
## last changed 2009 12 02, 2010-12-01 ms()
## 2010-12-01 check for overlapping columns
## 2011-02-08 rewritten to include polygonX and transectX detectors, and simplified
## 2011-03-18 output to unmarked
## 2011-03-21 'by' argument
############################################################################################

#----------------------------------------------------------------------------------------------------
# Dimensions   2       2      3	         3      2        2          3       3+      3+
#----------------------------------------------------------------------------------------------------
#              single  multi  proximity  count  polygonX transectX  signal  polygon transect
# single	&#	&	*	   *      NA       NA         NA      NA      NA
# multi		&#	&	*	   *      NA       NA         NA      NA      NA
# proximity	&#	&	*	   *      NA       NA         NA      NA      NA
# count		&#@	&@	@	   *      NA       NA         NA      NA      NA
# polygonX	&#	&	*	   *      &$       NA         NA      NA      NA
# transectX	&#	&	*	   *      NA       &$         NA      NA      NA
# signal        &#	&	*	   @      @        @          $       NA      NA
# polygon     	&#@~	&@~	@~	   *~     @        NA         NA       *      NA
# transect     	&#@~	&@~	@~	   *~     NA       @          NA      NA       *

#----------------------------------------------------------------------------------------------------
#  * no loss of data
#  # must choose among animals (more than one animal)
#  & must choose among traps (animal caught more than once)
#  @ reduce to binary presence
#  $ apply rule for combining signals or locations (first, last, random, min, max, mean)
#  ~ form new point detectors from mean of vertices (assumes symmetry)
#  NA not feasible
############################################################################################

    poly2point <- function (object, detector = 'count') {
        if (!detector(object) %in% c('polygon','polygonX'))
            stop ("requires 'polygon' input")
        if (detector %in% .localstuff$polydetectors)
            stop ("requires non-polygon, non-transect output")
        temp <- split(object, polyID(object))
        temp <- lapply(temp, function(df) apply(df,2,mean))
        temp1 <- t(abind(temp, along=2))
        dimnames(temp1) <- list(levels(polyID(object)), c('x','y'))
        temp <- data.frame(temp1, row.names=NULL)
        class (temp)   <- c('traps', 'data.frame')
        detector(temp) <- detector
        usage(temp)    <- usage(object)
        covariates(temp) <- covariates(object)
        attr(temp,'spacex') <- 100 * (searcharea(object)/nrow(temp))^0.5
        attr(temp,'spacey') <- attr(temp,'spacex')
        temp
    }

    transect2point <- function (object, detector = 'count') {
        if (!detector(object) %in% c('transect','transectX'))
            stop ("requires 'transect' input")
        if (detector %in% c('transect', 'transectX'))
            stop ("requires non-transect output")
        temp <- split(object, transectID(object))
        temp <- lapply(temp, function(df) apply(df,2,mean))
        temp1 <- t(abind(temp, along=2))
        dimnames(temp1) <- list(levels(transectID(object)), c('x','y'))
        temp <- data.frame(temp1, row.names=NULL)
        class (temp)   <- c('traps', 'data.frame')
        detector(temp) <- detector
        usage(temp)    <- usage(object)
        covariates(temp) <- covariates(object)
        attr(temp,'spacex') <- mean(transectlength(object))/2   ## arbitrary
        attr(temp,'spacey') <- attr(temp,'spacex')
        temp
    }

## function to make list in which each component is a
## subset of occasions (for use in reduce.capthist)
## MGE 2011-03-10

split.by <- function (x, by) {
    if ((length(x) == 1) & (x[1] > 1))
        x <- 1:x
    if (by < 1)
        stop ("invalid 'by' argument")
    index <- 1:length(x)
    gp <- trunc((index-1)/by) + 1
    split (index, gp)
}

reduce.capthist <- function (object, columns = NULL, outputdetector =
    detector(traps(object)), select='last', dropunused = TRUE, verify = TRUE, sessions =
    NULL, by = 1, ...) {

    # columns - list, each element gives occasions to include in new capthist

    seltrap <- function (y) {
        y <- t(y)                  # allow for occasion x trap matrix in proximity, count data
        y <- y[abs(y)>0]
        if (length(y)<1) y <- 0
        if (length(y) == 1) y
        else switch (select,
            first = head(y,1),    # first non-null
            last = tail(y,1),     # last non-null
            random = sample (size=1, y)    # random non-null, weighted by frequency in sample
        )
    }

    selused <- function (y) {
        y <- y[abs(y)>0]
        if (length(y)<1) y <- 0
        if (length(y) == 1) y
        else switch (select,
            first = head(y,1),    # first non-null
            last = tail(y,1),     # last non-null
            random = sample (size=1, y)    # random non-null, weighted by frequency in sample
        )
    }

    #----------------------------------------------------------------------------
    # functions applied to collapse a set of occasions 'occ' to a single occasion
    # result is a vector (single, multi detectors)
    fnused <- function (occ, fn) {
        if (length(occ)>0) {
            temp <- usage(traps(object))[,occ,drop=F]
            apply (temp, 1, fn)
        }
        else NULL
    }
    #----------------------------------------------------------------------------

    collapse <- function (df) {
        ## reduce data frame to a single row
        if (nrow(df)>1) {
            df$alive <- rep(all(df$alive), nrow(df))
            allowedCriteria <- c('first','last','random')
            if (!(select %in% allowedCriteria))
                stop ("selection criterion for signal should be one of ",
                      paste(sapply(allowedCriteria, dQuote),collapse=','))
            index <- switch (select, first = 1, last = nrow(df),
                random = sample.int (nrow(df),1) )
            df <- df[index,,drop=FALSE]
        }
        df
    }
    #----------------------------------------------------------------------------

    if (ms(object)) {
        if (is.null(sessions)) sessions <- 1:length(object)
        temp <- lapply (object[sessions], reduce,
            columns = columns,
            outputdetector = outputdetector,
            select = select,
            dropunused = dropunused,
            verify = verify,
            by = by,
            ...)
        class(temp) <- c('list', 'capthist')
        if (length(temp) == 1) temp <- temp[[1]]
        return(temp)
    }
    else {
        if (tolower(by) == 'all')
            by <- ncol(object)
        polygons <- c('polygon','polygonX')
        transects <- c('transect','transectX')
        inputdetector <- detector(traps(object))
        ntrap <- ndetector(traps(object))  ## npoly if 'polygon' or 'transect'
        nrw <- nrow(object)
        cutval <- attr(object, 'cutval')
        if (is.null(columns)) {
##          columns <- as.list(1:ncol(object))
            columns <- split.by (1:ncol(object), by)
            if ((ncol(object) %% by) > 0)
                warning ("number of occasions is not a multiple of 'by'")
        }

        if (is.null(outputdetector)) outputdetector <- inputdetector
        if (!(outputdetector %in% .localstuff$validdetectors))
            stop ("'outputdetector' should be one of ",
                  paste(sapply(.localstuff$validdetectors, dQuote),collapse=','))
        if ((inputdetector != 'signal') && (outputdetector == 'signal'))
                stop ("cannot convert non-signal data to signal data")
        if ((!(inputdetector %in% polygons)) && (outputdetector %in% polygons))
                stop ("cannot convert non-polygon data to 'polygon' data")
        if ((!(inputdetector %in% transects)) && (outputdetector %in% transects))
                stop ("cannot convert non-transect data to 'transect' data")

        ####################################
        ## check columns
        for (i in length(columns):1) {
            occ <- columns[[i]]
            occ <- occ[occ %in% (1:ncol(object))]  ## discard nonexistent occ
            if (length(occ)==0)
                columns[[i]] <- NULL
            else
                columns[[i]] <- occ
        }
        cumocc <- numeric(0)
        for (i in length(columns):1) {
            if (any (columns[[i]] %in% cumocc))
                warning ("new columns overlap")
            cumocc <- c(cumocc, columns[[i]])
        }
        nnew <- length(columns)

        ################################

        df <- data.frame(
            trap = trap(object, names = F),
            occ = occasion(object),
            ID = animalID(object, names = F),
            alive = alive(object))

        if (outputdetector %in% c(polygons, transects)) {
            df$x <- xy(object)[,1]
            df$y <- xy(object)[,2]
        }
        if (outputdetector %in% c('signal'))
            df$signal <- signal(object)
        ################################

        validcols <- unlist(columns)
        newcols <- rep(1:nnew, sapply(columns,length))
        newcols <- factor(newcols)  ## added 2011-02-16
        df$newocc <- newcols[match(df$occ, validcols)]
        if (dropunused) {
            df$newocc <- factor(df$newocc)
            nnew <- length(levels(df$newocc))
        }
        df <- df[!is.na(df$newocc),]                   ## drop null obs
        df$newID <- factor(df$ID)                      ## assign newID
        if (outputdetector %in% .localstuff$exclusivedetectors) {
            ID.occ <- interaction(df$ID, df$newocc)
            dflist <- split(df, ID.occ)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
        }
        if (outputdetector %in% c('single')) {
            occ.trap <- interaction(df$newocc,df$trap)
            dflist <- split(df, occ.trap)
            dflist <- lapply(dflist, collapse)
            df <- do.call(rbind, dflist)
        }
        df$newID <- factor(df$ID)                     ## re-assign newID
        validrows <- (1:nrow(object)) %in% df$ID

        if (outputdetector %in% .localstuff$exclusivedetectors) {
            alivesign <- df$alive*2 - 1
            tempnew <- matrix(0, nrow = sum(validrows), ncol = nnew)
            tempnew[cbind(df$newID, df$newocc)] <- df$trap * alivesign
        }
        else {
            df$trap <- factor(df$trap, levels=1:ntrap)
            tempnew <- table(df$newID, df$newocc, df$trap)
            alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$trap),all)
            alivesign[is.na(alivesign)] <- TRUE
            alivesign <- alivesign * 2 - 1
            if (! (outputdetector %in% .localstuff$countdetectors)
                && (length(tempnew)>0)) {
                ## convert 'proximity' and 'signal' to binary
                tempnew[tempnew>0] <- 1
            }
            tempnew <- tempnew * alivesign
        }

        class(tempnew) <- 'capthist'
        session(tempnew) <- session(object)

        ################################
        ## traps

        if ((inputdetector %in% polygons) && !(outputdetector %in% polygons))
            traps(tempnew) <- poly2point(traps(object))
        else
            if ((inputdetector %in% transects) && !(outputdetector %in% transects))
                traps(tempnew) <- transect2point(traps(object))
        else
            traps(tempnew) <- traps(object)  # drop traps not used on occasions
        detector(traps(tempnew)) <- outputdetector

        ################################
        ## covariates and ancillary data

        if (!is.null(covariates(object)))
             covariates(tempnew) <- covariates(object)[validrows,,drop=F]

        detectorder <- order(df$trap, df$newocc,df$ID)  ## CHECK!
        if (outputdetector %in% c(polygons, transects))
            xy(tempnew) <- df[detectorder,c('x','y'),drop=FALSE]
        if (outputdetector %in% c('signal')) {
            signal(tempnew) <- df$signal[detectorder]
            attr(tempnew, 'cutval') <- cutval
        }

        ################################
        ## usage
        if (nrow(tempnew) > 0)
            dimnames(tempnew)[[1]] <- 1:nrow(tempnew)  ## temporary, for animalID in subset
        if (!is.null(usage(traps(tempnew)))) {
            usagematrix <- unlist(sapply (columns, fnused, max))
            usagematrix <- matrix(usagematrix, nrow = nrow(traps(tempnew)))
            usage(traps(tempnew)) <- usagematrix
            if (dropunused) {
                OK <- apply(usagematrix, 1, sum) > 0
                tempnew <- subset(tempnew, traps = OK)
            }
        }
        tempnew[is.na(tempnew)] <- 0

        ################################
        ## dimnames
        if (nrow(tempnew) > 0) {
            indices <- (1:length(validrows))[validrows]
            rowname <- rownames(object)[indices]
        }
        else
            rowname <- NULL
        if (length(dim(tempnew)) == 3)
            dimnames(tempnew) <- list(rowname,1:nnew,NULL)   # renew numbering
        else
            dimnames(tempnew) <- list(rowname,1:nnew)

        if (verify) verify(tempnew, report=1)

        tempnew
    }
}
############################################################################################
