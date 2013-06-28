############################################################################################
## package 'secr'
## make.capthist.R
## 2010 05 02 (transferred from methods.R) 2010 05 03, 2010-11-21, 2011-01-21
## 2012-02-08 signalnoise
## 2012-02-09 revamped sorting
## 2012-02-12 finished tidy up related to signalframe
## 2012-10-19 telemetry detector type
############################################################################################

make.capthist <- function (captures, traps, fmt = 'trapID', noccasions = NULL,
    covnames = NULL, bysession = TRUE, sortrows = TRUE, cutval = NULL, tol = 0.01,
    noncapt = 'NONE', signalcovariates = NULL)

# captures is a dataframe with the structure:
# fmt = 'trapID'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	TrapID
#   column 5    Signal   (optional)
#   column 6    Noise    (optional)

# fmt = 'XY'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	x
#   column 5	y
#   column 6    Signal    (optional)
#   column 7    Noise     (optional)

{
    session <- captures[,1]
    sessionlevels <- unique(session)  ## retains order

    ## session <- factor(session) ## 2010 04 01 automatically sorts levels
    ## use numeric sort if appropriate 2010 05 02
    if (suppressWarnings( all(!is.na(as.numeric(sessionlevels)))))
        sessionlevels <- sessionlevels[order(as.numeric(sessionlevels))]
    else
        sessionlevels <- sort(sessionlevels)
    session <- factor(session, levels=sessionlevels)
    MS <- bysession & ( length(sessionlevels) > 1)
    if (MS) {  # recursive call of make.capthist
        capturelist <- split (captures, session)
        nsession <- length(capturelist)
        traplist <- inherits(traps, 'list')
        occvector <- length(noccasions)>1
        if (traplist & (length(traps) != nsession))
            stop ("multi-session 'traps' list does not match 'captures'")
        if (occvector & (length(noccasions) != nsession))
            stop ("requires one element in 'noccasions' for each session")
        capthist <- vector('list', nsession)
        for (i in 1:nsession) {
            if (traplist)  trps <- traps[[i]] else trps <- traps
            if (occvector) nocc <- noccasions[i] else nocc <- noccasions
            capthist[[i]]  <- make.capthist (
                captures = capturelist[[i]],
                traps = trps,
                fmt = fmt,
                noccasions = nocc,
                covnames = covnames,
                bysession = FALSE,         ## 2010 04 01
                sortrows = sortrows,
                cutval = cutval,
                tol = tol)

        }
        names(capthist) <- levels(session)
        class(capthist) <- c('list','capthist')
        capthist
    }

    else ## single-session call
    {
        if (!(fmt %in% c('trapID','XY')))
            stop ("capture format not recognised")
        if (fmt!='trapID') {
          if (ncol(captures)<5)
              stop ("too few columns in capture matrix")
          if (detector(traps) %in% c('polygon','polygonX','telemetry')) {
              captTrap <- xyinpoly(captures[,4:5], traps)
              if (any(captTrap==0)) {
                  captures <- captures[captTrap>0,]  ## first! 2010-11-17
                  captTrap <- captTrap[captTrap>0]
                  warning ("detections with coordinates outside ",
                           "polygon(s) were dropped")
              }
          }
          else if (detector(traps) %in% c('transect','transectX')) {
              captTrap <- xyontransect(captures[,4:5], traps, tol)
              if (any(captTrap==0)) {
                  captTrap <- captTrap[captTrap>0]
                  captures <- captures[captTrap>0,]
                  warning ("detections with coordinates not on ",
                           "any transect were dropped")
              }
          }
          else {
              trapID    <- interaction(traps$x, traps$y)
              captTrap  <- match(interaction(captures[,4], captures[,5]), trapID)
              if (any(is.na(captTrap)))
                  stop ("failed to match some capture locations ",
                        "to detector sites")
          }
        }
        else {
            if (detector(traps) %in% .localstuff$polydetectors)
                stop ("use fmt XY to input detections from polygons or transects")
            captTrap <- match(captures[,4], row.names(traps))
            if (any(is.na(captTrap)))
                stop ("failed to match some capture locations ",
                      "to detector sites")
        }


        #  if (bysession & ( length(levels(session)) > 1)) {
        #    captures[,2] <- interaction(session, captures[,2], drop = TRUE)
        #  }

        nocc      <- max(abs(captures[,3]))
        nocc      <- ifelse (is.null(noccasions), nocc, noccasions)
        if (is.null(detector(traps)))
            stop ("'traps' must have a detector type e.g. 'multi'")
        if (is.null(cutval) && detector(traps)  %in% c('cue','signal','signalnoise'))
            stop ("missing 'cutval' (signal threshold) for signal data")

        wout <- NULL
        ID   <- NULL

        uniqueID <- unique(captures[,2])

        ## optional row sort 2009 09 26, tweaked 2010 05 01
        if (sortrows) {
            if (suppressWarnings( all(!is.na(as.numeric(uniqueID)))))
                rowOrder <- order(as.numeric(uniqueID))
            else
                rowOrder <- order (uniqueID)
            uniqueID <- uniqueID[rowOrder]
#                if (length(dim(wout))==3)
#                    wout[,,] <- wout[rowOrder,,]
#                else
#                    wout[,] <- wout[rowOrder,]
#                dimnames(wout)[[1]] <- dimnames(wout)[[1]][rowOrder]
        }

        captID <- as.numeric(factor(captures[,2], levels=uniqueID))
        nID    <- length(uniqueID)
        detectionOrder <- order(captTrap, abs(captures[,3]), captID)

        dim3 <- detector(traps) %in% .localstuff$detectors3D
        if (dim3) {
            w <- array (0, dim=c(nID, nocc, ndetector(traps)))
            ## 2011-03-19, 2011-03-27
            ## drop rows if dummy input row indicates no captures
            if (any(uniqueID == noncapt)) {
                if (any(uniqueID != noncapt))
                    stop ("cannot combine data and noncapt")
                w <- w[FALSE, , , drop = FALSE]
                dimnames(w) <- list(NULL, 1:nocc, 1:ndetector(traps))
            }
            else {
                dimnames(w) <- list(1:nID, 1:nocc, 1:ndetector(traps))
                temp <- table (captID, abs(captures[,3]), captTrap)
                d <- dimnames(temp)
                d <- lapply(d, as.numeric)
                w[d[[1]], d[[2]], d[[3]]] <- temp

                ## fix to retain deads 2010-08-06
                dead <- captures[,3]<0
                deadindices <- cbind(captID[dead], abs(captures[dead,3]), captTrap[dead])
                w[deadindices] <- w[deadindices] * -1
                #################################
            }

        }
        else {
            w     <- matrix(0, nrow = nID, ncol = nocc)
            ## 2011-03-19, 2011-03-27
            ## drop rows if dummy input row indicates no captures
            if (any(uniqueID == noncapt)) {
                if (any(uniqueID != noncapt))
                    stop ("cannot combine data and noncapt")
                w <- w[FALSE, , drop = FALSE]
            }
            else {
                ## adjusted 2009 08 13 to ensure first occurrence selected when
                ## more than one per occasion
                indices <- cbind(captID, abs(captures[,3]))
                values <- captTrap * sign(captures[,3])
                ord <- order (captID, 1:length(captID), decreasing=TRUE)
                ## drop=F is critical in next line to ensure retains dim2
                w[indices[ord,, drop=F]] <- values[ord]
            }
        }

        wout <- abind(wout, w, along=1)
        dimnames(wout)[[2]] <- 1:nocc

        if (nrow(wout) > 0) {

            dimnames(wout)[[1]] <- uniqueID

            ## check added 2010 05 01
            if (any(is.na(wout))) {
                rw <- row(wout)[is.na(wout)]
                if (dim3)
                    print (wout[rw,,,drop=F])
                else
                    print (wout[rw,,drop=F])
                stop ("missing values not allowed")
            }

            ## code to input permanent individual covariates if these are present
            zi <- NULL
            startcol <- ifelse (fmt=='trapID', 5, 6)
            if (detector(traps) %in% c('signal')) startcol <- startcol+1
            if (detector(traps) %in% c('signalnoise')) startcol <- startcol+1
            if (ncol(captures) >= startcol)
                zi <- as.data.frame(captures[,startcol:ncol(captures), drop=F])
            if (!is.null(zi)) {
                # find first match of ID with positive value for each covar
                temp <- zi[1:length(uniqueID),,drop=FALSE]
                temp[,] <- NA

                for (j in 1:ncol(zi)) {
                    nonmissing <- function(x) x[!is.na(x)][1]
                    ## levels=uniqueID to retain sort order
                    tempj <- split(zi[,j], factor(captures[,2], levels=uniqueID))
                    tempj2 <- sapply(tempj, nonmissing)
                    ## 2011-01-20 only convert to factor if character
                    if (is.character(tempj2))
                       tempj2 <- factor(tempj2)
                    temp[,j] <- tempj2
                    rownames(temp) <- names(tempj)          ## 2010 02 26
                }
                if (is.null(covnames)) names(temp) <- names(zi)
                else {
                    if (ncol(temp) != length(covnames))
                        stop ("number of covariate names does not match")
                    names(temp) <- covnames
                }
                ## added 'FALSE' 2010 02 24 -
## probably redundant given new use of uniqueID 2012-02-09
##                if (sortrows) temp <- temp[rowOrder,,drop = FALSE]
                attr(wout,'covariates') <- temp
            }
            else attr(wout,'covariates') <- data.frame()

        }

        class (wout) <- 'capthist'
        traps(wout) <- traps
        session(wout)  <- as.character(captures[1,1])

        if (nrow(wout) > 0) {

            if (detector(traps) %in% .localstuff$polydetectors) {
                ## 2011-01-21
##                xy <- captures[order(captTrap, captures[,3],captures[,2]),4:5]
                xy <- captures[detectionOrder,4:5]
                names(xy) <- c('x','y')
                attr(wout,'detectedXY') <- xy
            }
            if (detector(traps) %in% c('signal','signalnoise')) {
                if (is.null(cutval))
                    stop ("missing value for signal threshold")
                if (fmt=='XY')
                    signl <- captures[,6]
                else
                    signl <- captures[,5]
                signl <- signl[detectionOrder]
                signal(wout) <- signl
                if (detector(traps) %in% 'signalnoise') {
                    if (fmt=='XY')
                        nois <- captures[,7]
                    else
                        nois <- captures[,6]
                    nois <- nois[detectionOrder]
                    noise(wout) <- nois
                }
                if (!is.null(signalcovariates)) {
                    if (!all(signalcovariates %in% names(captures)))
                        stop ("missing signal covariate(s)")

                    attr(wout, 'signalframe') <- cbind(attr(wout, 'signalframe'),
# captures[,signalcovariates])
# need to maintain order 2012-09-15
                        captures[detectionOrder,signalcovariates])
                }
                attr(wout, 'cutval')   <- cutval
                ## dropunused = FALSE? 2012-01-11

                ## apply cutval
                wout <- subset(wout, cutval = cutval)
            }
        }

        wout
    }   ## end of single-session call
}
############################################################################################

## repeat the detections in x a number of times
## 'times' must be sorted by trap, occasion and animalID
## not sure this is useful! 2010-11-21
rep.capthist <- function (x, times) {
    if (!inherits(x, 'capthist')) {
        stop ("requires 'capthist' object")
    }
    if (ms(x)) {
        if (!is.list(times))
            stop ("multi-session capthist; 'times' should be a list")
        for (i in 1:length(x))
            x[[i]] <- rep.capthist(x, times[[i]])
        x
    }
    else {
        if (detector(traps(x)) %in% c('single','multi')) {
            stop ("cannot replicate data from this detector type")
        }
        if (sum(abs(x)) != length (times)) {
            stop ("'times' must be equal in length to the number of detections")
        }
        if (any(times < 1)) {
            stop ("'times' must be a positive integer")
        }
        ID <- match(animalID(x), dimnames(x)[[1]])
        occ <- occasion(x)
        trp <- as.numeric(as.character(trap(x)))

        x[cbind(ID, occ, trp)] <- times
        if (detector(traps(x))=='proximity') {
            detector(traps(x)) <- 'count'
            warning ("detector type changed to 'count'")
        }
        # individual covariates unchanged
        x
    }
}
