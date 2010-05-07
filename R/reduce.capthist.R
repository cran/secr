############################################################################################
## package 'secr'
## reduce.capthist.R
## last changed 2009 12 02
#####################################################################################################

#----------------------------------------------------------------------------------------------------
# Dimensions   2       2      3	         3      3              3             3       3+      3+
#----------------------------------------------------------------------------------------------------
#              single  multi  proximity  count  quadratbinary  quadratcount  signal  polygon transect
# single	&#	&	*	   *        *           *            NA      NA      NA
# multi		&#	&	*	   *        *           *            NA      NA      NA
# proximity	&#	&	*	   *        *           *            NA      NA      NA
# count		&#@	&@	@	   *        @           *            NA      NA      NA
# quadratbinary	&#	&	*	   *        *           *            NA      NA      NA
# quadratcount	&#@	&@	@	   *        @           *            NA      NA      NA
# signal        &#	&	@	   @        @           @            $       NA      NA
# polygon     	&#~	&~	@~	   *~       @~          *~           NA       *	     NA
# transect     	&#~	&~	@~	   *~       @~          *~           NA      NA	     *

#----------------------------------------------------------------------------------------------------
#  * no loss of data
#  # must choose among animals (more than one animal)
#  & must choose among traps (animal caught more than once)
#  @ reduce to binary presence
#  $ apply rule for combining signals (first, last, random, min, max, mean)
#  ~ form new point detectors from mean of vertices (assumes symmetry)
#  NA not feasible
############################################################################################

    poly2point <- function (object, detector = 'count') {
        if (detector(object)!='polygon') stop ('need polygon input')
        if (detector=='polygon') stop ('need non-polygon output')
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
        if (detector(object)!='transect') stop ('need transect input')
        if (detector=='transect') stop ('need non-transect output')
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

reduce.capthist <- function (object, columns = NULL, outputdetector = 
    detector(traps(object)), select='last', dropunused = TRUE, verify = TRUE, sessions = 
    NULL, ...) {

# columns - list, each element gives occasions to include in new capthist

    seltrap <- function (y) {
        y <- t(y)                  # allow for occason x trap matrix in proximity, count data
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

    selsignal <- function (y) {
        if (length(y) <= 1) y
        else switch (select, 
            first = head(y,1),    # first non-null
            last = tail(y,1),     # last non-null
            min = min(y),
            max = max(y),
            mean = mean(y), 
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
    fnp2multi <- function (occ, fn) {
        occ <- occ[occ <= ncol(temp3D)]
        temp <- temp3D[,occ,,drop=F]
        temp[,,] <- aperm(apply(temp, 1:2, function(x) (1:length(x)) * x), c(2,3,1))  
        apply (temp, 1, fn)
    }
    fnp2single <- function (occ, fn) {
        occ <- occ[occ <= ncol(temp3D)]
        temp <- temp3D[,occ,,drop=F]
        temp[,,] <- apply(temp, 2:3, function(x) (1:length(x)) * x)  ## animal number
        randdropdupl <- function (x) {
            z <- (1:length(x))[abs(x)>0]
            if (length(z)>1) z <- sample(size=1, z)
            x[-z] <- 0
            x
        }
        
        temp <- apply (temp, 1:2, randdropdupl)  ## resolve multiple records of animal on one occasion
        temp <- aperm (temp, c(2,3,1))
        temp <- apply (temp, 2:3, randdropdupl)  ## resolve multiple records at one trap on one occasion
        temp <- apply (temp, 3, fn)   ## for each trap, the ID of the animal
        id <- 1:nrow(temp3D) 
        intrap <- match(id,temp)
        ifelse(is.na(intrap), 0, intrap) 
    }
    fnx <- function (x) (1-2*any(x<0)) * sum(abs(x))
    fncount <- function (occ) {
        occ <- occ[occ <= ncol(temp3D)]
        apply (temp3D[,occ,,drop=F], c(1,3), fnx)
    }
    #----------------------------------------------------------------------------

    if (inherits(object,'list')) {
        if (is.null(sessions)) sessions <- 1:length(object)
        temp <- lapply (object[sessions], reduce, 
            columns = columns, 
            outputdetector = outputdetector, 
            select = select, 
            dropunused = dropunused,
            verify = verify, ...)
        class(temp) <- c('list', 'capthist')
        if (length(temp) == 1) temp <- temp[[1]] 
        return(temp)
    }
    else {
        inputdetector <- detector(traps(object))
        if (is.null(columns)) columns <- as.list(1:ncol(object))
        if (is.null(outputdetector)) outputdetector <- inputdetector
        if (!(outputdetector %in% .localstuff$validdetectors))   
            stop ('Unrecognised output detector type')
        if ((inputdetector != 'signal') && (outputdetector == 'signal'))
                stop ('cannot convert non-signal data to signal data')
        if ((inputdetector != 'polygon') && (outputdetector == 'polygon'))
                stop ('cannot convert non-polygon data to polygon data')
        if ((inputdetector != 'transect') && (outputdetector == 'transect'))
                stop ('cannot convert non-transect data to transect data')

        nrw <- nrow(object)

        for (i in length(columns):1) {
            occ <- columns[[i]]
            occ <- occ[occ <= ncol(object)]
            if (length(occ)==0) columns[[i]] <- NULL
        } 
        nnew <- length(columns)
        ntrap <- ndetector(traps(object))  ## npoly if 'polygon' or 'transect'
        cutval <- attr(object, 'cutval')

        if (outputdetector == 'signal') {
            fnsignal <- function (occ) {
                OK <- occasion(object) %in% occ
                sig <- tapply(signal(object)[OK], list(animalID(object)[OK], trap(object)[OK]), selsignal)
                sig[!is.na(sig)]
            }
            if (!(select %in% c('first','last','min','max','mean','random'))) 
                stop ('unrecognised selection criterion for signal') 
            tempsignal <- as.numeric(sapply (columns, fnsignal))
        }

        if (inputdetector %in% c('single', 'multi')) {
            temp3D <- array(0, dim=c(dim(object), ntrap))
            trp <- as.numeric(object)
            OK <- abs(trp) > 0
            indices <- cbind(as.numeric(row(object)), as.numeric(col(object)), trp)[OK,]
            temp3D[abs(indices)] <- sign(trp[OK])
        }
        else temp3D <- as.array(object)
        if (outputdetector == 'multi') {
            if (!(select %in% c('first','last','random'))) 
                stop ('unrecognised selection criterion for trap') 
            tempnew <- sapply (columns, fnp2multi, seltrap)
        }
        else if (outputdetector == 'single') {
            if (!(select %in% c('first','last','random'))) 
                stop ('unrecognised selection criterion for single-catch trap') 
            tempnew <- sapply (columns, fnp2single, seltrap)
        }
        else { 
            tempnew <- sapply (columns, fncount)
            tempnew <- array (tempnew, dim=c(nrw, ntrap, nnew))
            tempnew <- aperm(tempnew, c(1,3,2))
            if (! (outputdetector %in% c('count', 'quadratcount','polygon','transect')) && (length(tempnew)>0)) {
                ## convert to binary
                tempnew <- sign(tempnew) * (abs(tempnew) > 0)
            }
        }
    
        ## drop empty rows
        if (nrow(tempnew)>0) {
            tempnew <- tempnew * 1    ## numeric not logical
            OK <- apply(abs(tempnew), 1, sum) > 0       ## does not apply cutval 2009 12 01
            if (outputdetector %in% c('single', 'multi')) 
                tempnew <- tempnew[OK,, drop = F]
            else 
                tempnew <- tempnew[OK,,, drop = F]
        }
        else OK <- logical(0)
        class(tempnew) <- 'capthist'
        session(tempnew) <- session(object)
        covariates(tempnew) <- covariates(object)[OK,,drop=F]
        if ((inputdetector == 'polygon') && !(outputdetector == 'polygon')) 
            traps(tempnew) <- poly2point(traps(object))
        else if ((inputdetector == 'transect') && !(outputdetector == 'transect')) 
            traps(tempnew) <- transect2point(traps(object))
        else
            traps(tempnew) <- traps(object)  # drop traps not used on occasions 
        detector(traps(tempnew)) <- outputdetector

        if (outputdetector %in% c('polygon','transect')) {
            OK <- occasion(object) %in% unlist(columns)
            attr(tempnew, 'detectedXY') <- xy(object)[OK,]
        }
        if (outputdetector %in% c('signal')) {
            attr(tempnew, 'signal') <- tempsignal 
            attr(tempnew, 'cutval') <- cutval
        }
        dimnames(tempnew)[[1]] <- 1:nrow(tempnew)    ## temporary, for animalID in subset
        if (!is.null(usage(traps(tempnew)))) {
            usagematrix <- unlist(sapply (columns, fnused, max))
            usagematrix <- matrix(usagematrix, nr=nrow(traps(tempnew)))
            usage(traps(tempnew)) <- usagematrix
            if (dropunused) {
                OK <- apply(usagematrix, 1, sum) > 0    
                tempnew <- subset(tempnew, traps=OK)

            }
        }
        tempnew[is.na(tempnew)] <- 0
    
        if (nrow(tempnew)>0)     
            rowname <- 1:nrow(tempnew)
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
