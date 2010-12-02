############################################################################################
## package 'secr'
## verify.R
## 2009 09 18, 2009 09 19, 2009 09 20 2009 10 02 2009 11 05
## 2009 11 13
## 2010 05 02 removed erroneous ref to 'areabinary' detector
############################################################################################

verify <- function (object, report, ...) UseMethod("verify")

verify.default <- function (object, report, ...) {
  cat ('no verify method for objects of class', class(object), '\n')
}
############################################################################################

overlapcells <- function (xy) {
    vertexinside <- function (a,b) {
        OK <- FALSE
        for (k in 1:4) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (a[k,]),
                as.integer (0),
                as.integer (3),
                as.integer (5),
                as.double (b),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (b[k,]),
                as.integer (0),
                as.integer (3),
                as.integer (5),
                as.double (a),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        OK
    }
    spacex <- attr(xy, 'spacex')
    spacey <- attr(xy, 'spacey')
    fuzz <- 1e-10
    spx2 <- spacex/2 - fuzz
    spy2 <- spacey/2 - fuzz
    xy <- as.matrix(xy)
    nr <- nrow(xy)
    if (nr<2)
        FALSE
    else {
        pixel <- matrix(nc=2, c(-spx2,-spx2,spx2,spx2,-spx2,-spy2,spy2,spy2,-spy2,-spy2))
        overlap <- matrix(FALSE, nrow=nr, ncol=nr)
        for (i in 1:(nr-1))
            for (j in (i+1):nr)
                {
                    verti <- t(apply(pixel, 1, function (x) xy[i,] + x))
                    vertj <- t(apply(pixel, 1, function (x) xy[j,] + x))
                    if ((length(verti)>0) && (length(vertj)>0))
                    overlap[i,j] <- vertexinside(verti,vertj)
                }
        any (overlap, na.rm=T)
    }
}
############################################################################################

overlappoly <- function (xy, polyID) {
    vertexinside <- function (a,b) {
        OK <- FALSE
        n.a <- nrow(a)
        n.b <- nrow(b)
        a <- as.matrix(a)
        b <- as.matrix(b)
        for (k in 1:n.a) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (a[k,]),
                as.integer (0),
                as.integer (n.b-1),
                as.integer (n.b),
                as.double (b),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        for (k in 1:n.b) {
            temp <- .C('inside',  PACKAGE = 'secr',
                as.double (b[k,]),
                as.integer (0),
                as.integer (n.a-1),
                as.integer (n.a),
                as.double (a),
                result = integer(1))
            if (any(as.logical(temp$result))) OK <- TRUE
        }
        OK
    }

    lxy <- split (xy, polyID)
    nr <- length(lxy)
    if (nr<2)
        FALSE
    else {
        overlap <- matrix(FALSE, nrow=nr, ncol=nr)
        for (i in 1:(nr-1))
            for (j in (i+1):nr)
                {
                    overlap[i,j] <- vertexinside(lxy[[i]], lxy[[j]])
                }
        any (overlap, na.rm=T)
    }
}
############################################################################################

xyinpoly <- function (xy, trps) {
    ptinside <- function (i,k) {
        ## is point i inside poly k?
        polyxy <- as.matrix(lxy[[k]])
        nr <- nrow(polyxy)
        temp <- .C('inside',  PACKAGE = 'secr',
            as.double (xy[i,]),
            as.integer (0),
            as.integer (nr-1),
            as.integer (nr),
            as.double (polyxy),
            result = integer(1))
        as.logical(temp$result)
    }
    lxy <- split (trps, polyID(trps))
    firstinside <- function (i) {
        for (k in 1:length(lxy))
            if (ptinside(i,k)) return(k)
        0
    }
    sapply(1:nrow(xy), firstinside)
}
############################################################################################

xyontransect <- function (xy, trps, tol=0.01) {
    ptontransect <- function (i,k) {
        ## is point i on transect k?
        transectxy <- as.matrix(lxy[[k]])
        nr <- nrow(transectxy)
        temp <- .C('ontransect',  PACKAGE = 'secr',
            as.double (xy[i,]),
            as.integer (0),
            as.integer (nr-1),
            as.integer (nr),
            as.double (transectxy),
            as.double (tol),
            result = integer(1))
        as.logical(temp$result)
    }
    lxy <- split (trps, transectID(trps))
    firsttransect <- function (i) {
        for (k in 1:length(lxy))
            if (ptontransect(i,k)) return(k)
        0
    }
    sapply(1:nrow(xy), firsttransect)
}
############################################################################################

verify.traps <- function (object, report = 2, ...) {

## Check internal consistency of 'traps' object
##
## -- Number of rows in dataframe of detector covariates differs expected
## -- Number of detectors in usage matrix differs from expected
## -- Occasions with no used detectors

    if (!inherits(object, 'traps')) {
         stop ("object must be of class 'traps'")
    }

    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report,1))
        anyerrors <- any(sapply(temp, function(x) x$errors))
        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {

        single <- detector(object) %in% c('single')
##        area <- detector(object) %in% c('quadratcount','quadratbinary')
        area <- FALSE
        poly <- detector(object) %in% c('polygon')

        usagedetectorsOK <- TRUE
        usagenonzeroOK <- TRUE
        areaOK <- TRUE
        polyIDOK <- TRUE

        if (!is.null(covariates(object)))
            if ((ncol(covariates(object)) == 0 ) |
                (nrow(covariates(object)) == 0 )) covariates(object) <- NULL

        ## 1
        trapNAOK <- !any(is.na(object))

        ## 2
        trapcovariatesOK <- ifelse (is.null(covariates(object)),
            TRUE, nrow(covariates(object)) == nrow(object))

        ## 'usage' of traps
        if (!is.null(usage(object))) {
            ## 3
            usagedetectorsOK <- nrow(usage(object)) == nrow(object)

            ## 4
            usagecount <- apply(usage(object),2,sum)
            usagenonzeroOK <- !any(usagecount == 0)
        }
        else usagecount <- rep(NA, ncol(object))

        ## 5
        if (area) {
            ## must have searchcell
            areaOK <- !is.na(searcharea(object))
            areaOK <- areaOK & !overlapcells(object)
        }
        else
        if (poly) {
            areaOK <- !overlappoly (object, polyID(object))
        }

        ## 6
        if (poly) {
            polyIDOK <- (length(polyID(object)) == nrow(object)) &
                is.factor(polyID(object))
        }

        errors <- !all(c(trapNAOK, trapcovariatesOK, usagedetectorsOK, usagenonzeroOK, areaOK, polyIDOK))

        if (report > 0) {
            if (errors) {
                if (!trapNAOK) {
                    cat ('Missing detector coordinates not allowed\n')
                }
                if (!trapcovariatesOK) {
                    cat ('Wrong number of rows in dataframe of detector covariates\n')
                    cat ('traps(capthist) :', nrow(traps(object)), 'detectors\n')
                    cat ('covariates(traps(capthist)) :', nrow(covariates(traps(object))), 'detectors\n')
                }
                if (!usagedetectorsOK) {
                    cat ('Conflicting number of detectors in usage matrix\n')
                    cat ('traps(capthist) :', nrow(traps(object)), 'detectors\n')
                    cat ('usage(traps(capthist)) :', nrow(usage(traps(object))), 'detectors\n')
                }
                if (!usagenonzeroOK) {
                    cat ("Occasions when no detectors 'used'\n")
                    cat ((1:length(usagecount))[usagecount==0], '\n')
                }
                if (!areaOK) {
                    cat ("Search areas overlap, or no search area specified \n")
                }
                if (!polyIDOK) {
                    cat ("Invalid polyID \n")
                }
            }
        }

        if ((report == 2) && !errors) cat('No errors found :-)\n')

        out <- list(errors = errors,
            trapNAOK = trapNAOK,
            trapcovariatesOK = trapcovariatesOK,
            usagedetectorsOK = usagedetectorsOK,
            usagenonzeroOK = usagenonzeroOK,
            areaOK = areaOK,
            usagecount = usagecount
        )

        invisible(out)

    }
}
############################################################################################

verify.capthist <- function (object, report = 2, tol = 0.01, ...) {

## Check internal consistency of 'capthist' object
##
## -- 'traps' component present
## -- verify(traps)
## -- No live releases
## -- Live detection(s) after reported dead
## -- More than one capture in single-catch trap(s)

## -- Number of rows in 'traps' object not compatible with reported detections
## -- Number of rows in dataframe of individual covariates differs from capthist
## -- Number of occasions in usage matrix differs from capthist
## -- Detections at unused detectors


    if (!inherits(object, 'capthist'))
        stop ("object must be of class 'capthist'")
    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report, 1))
        anyerrors <- any(sapply(temp, function(x) x$errors))
        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {

        ## preliminaries
##        dim3 <- detector(traps(object)) %in% c('proximity', 'count', 'quadratbinary', 'quadratcount', 'signal', 'polygon','transect')
        dim3 <- length(dim(object)) == 3
        count <- detector(traps(object)) %in% c('count', 'polygon','transect')
##        area <- detector(traps(object)) %in% c('quadratbinary', 'quadratcount')
        area <- FALSE
##        binary <- detector(traps(object)) %in% c('proximity', 'quadratbinary')
        binary <- detector(traps(object)) %in% c('proximity')
        single <- detector(traps(object)) %in% c('single')
        signal <- detector(traps(object)) %in% c('signal')
        poly <- detector(traps(object)) %in% c('polygon')
        transect <- detector(traps(object)) %in% c('transect')

        NAOK <- TRUE
        deadOK <- TRUE
        usageOK <- TRUE
        usageoccasionsOK <- TRUE
        usagedetectorsOK <- TRUE
        usagenonzeroOK <- TRUE
        detectorconflcts <- NULL
        singleOK <- TRUE
        binaryOK <- TRUE
        countOK <- TRUE
        cutvalOK <- TRUE
        signalOK <- TRUE
        xyOK <- TRUE
        xyinpolyOK <- TRUE
        xyontransectOK <- TRUE
        IDOK <- TRUE

        if (!is.null(covariates(object)))
            if ((ncol(covariates(object)) == 0 ) |
                (nrow(covariates(object)) == 0 ))
                covariates(object) <- NULL

        ## 1
        trapspresentOK <- !is.null(traps(object))

        ## standalone check of detectors
        if (trapspresentOK)
            trapcheck <- verify(traps(object), report = 0)  ## delay reporting
        else
            trapcheck <- list(errors=TRUE)

        ## 2
        trapsOK <- !trapcheck$errors

        ## 3
        if (length(object)==0)
            detectionsOK <- FALSE
        else  {

            detectionsOK <- sum(object[object>0]) > 0

            ## 4
            NAOK <- !any(is.na(object))

            ## 5
            if (signal) {
                ## must have cutval; deads not allowed
                if (length(attr(object,'cutval')) != 1) cutvalOK <- FALSE
                else if (any(signal(object) < attr(object,'cutval'))) cutvalOK <- FALSE
                if (length(signal(object)) != sum(abs(object)))
                    signalOK <- FALSE
            }
            ## 6
            else {
                fn <- function(x) {
                    if (dim3) x <- apply(x,1,min)
                    (min(x)<0) && (tail(x[x!=0],1)>0)
                }
                undead <- apply(object, 1, fn)
                deadOK <- !any(undead)
                if (!deadOK) {
                    if (dim3)
                        reincarnated <- object[undead,,, drop=F]
                    else
                        reincarnated <- object[undead,, drop=F]
                }
            }

            ## 7
            if (single) {
                fn <- function (x) duplicated(abs(x)[x!=0])
                multiple <- apply(object, 2, fn)
                singleOK <- !any(unlist(multiple))
            }

            ## 8
            if (binary) {
                ## must be binary
                multiples <- sum(abs(object)>1)
                binaryOK <- multiples == 0
            }

            ## 9
## blocked 2010-12-01 - no problem with 'dead' count
##            if (count) {
##                countOK <- all (object>=0)
##            }
        }

        ## 10
        if (poly | transect)
            detectornumberOK <- length(table(polyID(traps(object)))) == dim(object)[3]
        else
            detectornumberOK <- ifelse (dim3,
              dim(object)[3] == nrow(traps(object)),
              max(abs(object)) <= nrow(traps(object)))

        ## 11
        covariatesOK <- ifelse(is.null(covariates(object)),
            TRUE,
            nrow(covariates(object)) == nrow(object))

        ## is 'usage' of traps consistent with reported detections?
        if (!is.null(usage(traps(object)))) {
            conflcts <- 0

            ## 12
            usageoccasionsOK <- ncol(usage(traps(object))) == ncol(object)

            if (detectionsOK) {

                notused <- !usage(traps(object))  ## traps x occasions
                if (dim3) {
                    if (usagedetectorsOK && usageoccasionsOK) {
                        tempobj <- aperm(object, c(2,3,1))   ## occasion, traps, animal sKn
                        tempuse <- array(t(usage(traps(object))), dim=dim(tempobj))  ## replicated to fill...
                        conflcts <- (abs(tempobj)>0) && (tempuse==0)
                        tempobjmat <- array(tempobj[,,1], dim= dim(tempobj)[1:2])
                        occasion <- rep(row(tempobjmat), dim(tempobj)[3])
                        detector <- rep(col(tempobjmat), dim(tempobj)[3])
                        ID <- rep(rownames(object), rep(prod(dim(tempobj)[1:2]), nrow(object)))
                        detectorconflcts <- as.data.frame(cbind(ID,detector,occasion)[conflcts,])
                    }
                }
                else {
                    if (usagedetectorsOK && usageoccasionsOK) {
                        OK <- as.numeric(object)>0
                        occasion <- as.numeric(col(object))[OK]
                        ID <- row.names(object)[as.numeric(row(object))[OK]]
                        detector <- as.numeric(object)[OK]
                        conflcts <- notused[cbind(detector, occasion)] > 0
                        detectorconflcts <- as.data.frame(cbind(ID,detector,occasion)[conflcts,])
                    }
                }
            }

            ## 13
            usageOK <- sum(conflcts)==0

        }

        ## 14
        if (poly) {
            xy <- xy(object)
            xyOK <- nrow(xy) == sum(abs(object))
            inpoly <- xyinpoly(xy(object), traps(object))
            inpoly <- inpoly == trap(object, name = F)
            xyinpolyOK <- all(inpoly)
## check dropped 2010-11-17
##            ID <- as.numeric(animalID(object))   ## does this allow for alpha names?
##            IDOK <- all(table(ID) == apply(object,1,sum))
        }
        if (transect) {
            xy <- xy(object)
            ID <- as.numeric(animalID(object))   ## does this allow for alpha names?
            xyOK <- nrow(xy) == sum(abs(object))
            ontransect <- xyontransect(xy(object), traps(object), tol = tol)
            ontransect <- ontransect == trap(object, name = F)
            xyontransectOK <- all(ontransect)
            IDOK <- all(table(ID) == apply(object,1,sum))
        }

        errors <- !all(c(trapspresentOK, trapsOK, detectionsOK, NAOK, deadOK, singleOK, binaryOK,
            countOK, cutvalOK, signalOK, detectornumberOK, covariatesOK, usageoccasionsOK, usageOK, xyOK,
            xyinpolyOK, xyontransectOK, IDOK))

        if (report > 0) {
            if (errors) {
                cat ('Session', session(object), '\n')

                if (!trapspresentOK) {
                    cat ('No valid detectors\n')
                }
                if (!trapsOK) {
                    cat ('Errors in traps\n')
                    if (!trapcheck$trapNAOK) {
                        cat ('Missing detector coordinates not allowed\n')
                    }
                    if (!trapcheck$trapcovariatesOK) {
                        cat ('Wrong number of rows in dataframe of detector covariates\n')
                        cat ('traps(capthist) :', nrow(traps(object)), 'detectors\n')
                        cat ('covariates(traps(capthist)) :', nrow(covariates(traps(object))), 'detectors\n')
                    }
                    if (!trapcheck$usagedetectorsOK) {
                        cat ('Conflicting number of detectors in usage matrix\n')
                        cat ('traps(capthist) :', nrow(traps(object)), 'detectors\n')
                        cat ('usage(traps(capthist)) :', nrow(usage(traps(object))), 'detectors\n')
                    }
                    if (!trapcheck$usagenonzeroOK) {
                        cat ("Occasions when no detectors 'used'\n")
                        cat ((1:length(trapcheck$usagecount))[trapcheck$usagecount==0], '\n')
                    }
                }

                if (!detectionsOK) {
                    cat ('No live releases\n')
                }

                if (!NAOK) {
                    cat ('Missing values not allowed in capthist\n')
                }

                if (!deadOK) {
                    cat ('Recorded alive after dead\n')
                    print(reincarnated)
                }

                if (!singleOK) {
                    cat ('More than one capture in single-catch trap(s)\n')
                }

                if (!binaryOK) {
                    cat ('More than one detection per detector per occasion at proximity detector(s)\n')
                }

                if (!countOK) {
                    cat ('Count(s) less than zero\n')
                }

                if (!cutvalOK) {
                    cat ('Signal less than cutval or invalid cutval\n')
                }

                if (!signalOK) {
                    cat ('Signal attribute does not match detections\n')
                }

                if (!detectornumberOK) {
                    cat ('traps object incompatible with reported detections\n')
                    cat ('traps(capthist) :', nrow(traps(object)), 'detectors\n')
                    if (dim3)
                        cat ('capthist :', dim(object)[3], 'detectors\n')
                    else
                        cat ('capthist :', max(abs(object)), 'max(detector)\n')
                }

                if (!covariatesOK) {
                    cat ('Wrong number of rows in dataframe of individual covariates\n')
                    cat ('capthist :', nrow(object), 'individuals\n')
                    cat ('covariates(capthist) :', nrow(covariates(object)), 'individuals\n')
                }
                if (!usageoccasionsOK) {
                    cat ('Conflicting number of occasions in usage matrix\n')
                    cat ('capthist :', ncol(object), 'occasions\n')
                    cat ('usage(traps(capthist)) :', ncol(usage(traps(object))), 'occasions\n')
                }
                if (!usageOK) {
                    cat ("Detections at 'unused' detectors\n")
                    print(detectorconflcts)
                }
                if (!xyOK) {
                    cat ("Polygon detector xy coordinates of detections do not match counts\n")
                }
                if (!xyinpolyOK) {
                    cat ("XY coordinates not in polygon\n")
                    print (xy(object)[!inpoly,])
                }
                if (!xyontransectOK) {
                    cat ("XY coordinates not on transect\n")
                    print (xy(object)[!ontransect,])
                }
                if (!IDOK) {
                    cat ("Polygon detector mismatch between ID attribute and counts\n")
                }
            }

            if ((report == 2) && !errors) cat('No errors found :-)\n')

        }

        out <- list(errors = errors, trapcheck = trapcheck)
        if (!is.null(detectorconflcts)) out$detections.at.unused.detectors <- detectorconflcts
        invisible(out)
    }
}
############################################################################################

verify.mask <- function (object, report = 2, ...) {

## Check internal consistency of 'mask' object
##
## valid x and y coordinates
## nrow(covariates) = nrow(object)
## ...also look at attributes?

    if (!inherits(object, 'mask'))
        stop ("object must be of class 'mask'")

    if (inherits(object, 'list')) {
        temp <- lapply (object, verify, report = min(report, 1))
        anyerrors <- any(sapply(temp, function(x) x$errors))
        if ((report == 2) && !anyerrors)
            cat('No errors found :-)\n')
        invisible(list(errors = anyerrors, bysession = temp))
    }
    else {

        ## 1
        xyOK <- !(is.null(object$x) | is.null(object$y) | any(is.na(object)))
        xyOK <- xyOK && is.numeric(unlist(object))

        ## 2

        if (!is.null(covariates(object)))
            covariatesOK <- ifelse (nrow(covariates(object))>0,
            nrow(object) == nrow(covariates(object)), TRUE)
        else
            covariatesOK <- TRUE

        errors <- !all(c(xyOK, covariatesOK))

        if (report > 0) {
            if (errors) {
                ## cat ('Session', session(object), '\n')

                if (!xyOK) {
                    cat ('Invalid x or y coordinates in mask\n')
                }

                if (!covariatesOK) {
                    cat ('Number of rows in covariates(mask) differs from expected\n')
                }
            }

            if ((report == 2) && !errors) cat('No errors found :-)\n')
        }

        out <- list(errors = errors)
        invisible(out)
    }
}
############################################################################################

