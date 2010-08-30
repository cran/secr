############################################################################################
## package 'secr'
## make.capthist.R
## last changed 2010 05 02 (transferred from methods.R) 2010 05 03
############################################################################################

make.capthist <- function (captures, traps, fmt = 'trapID', noccasions = NULL,
    covnames = NULL, bysession = TRUE, sortrows = TRUE, cutval = NULL, tol = 0.01)

# captures is a dataframe with the structure:
# fmt = 'trapID'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	TrapID
#   column 5    Signal    (optional)

# fmt = 'XY'
#   column 1	Session
#   column 2	AnimalID
#   column 3	Occasion
#   column 4	x
#   column 5	y
#   column 6    Signal    (optional)

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
            stop ('traps list does not match capture sessions')
        if (occvector & (length(noccasions) != nsession))
            stop ('noccasions does not match capture sessions')
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
                cutval = cutval)

        }
        names(capthist) <- levels(session)
        class(capthist) <- c('list','capthist')
        capthist
    }

    else ## single-session call
    {
        if (!(fmt %in% c('trapID','XY'))) stop ('Capture format not recognised')
        if (fmt!='trapID') {
          if (ncol(captures)<5) stop ('Too few columns in capture matrix')
          if (detector(traps)=='polygon') {
              captTrap <- xyinpoly(captures[,4:5], traps)
              if (any(captTrap==0)) {
                  captTrap <- captTrap[captTrap>0]
                  captures <- captures[captTrap>0,]
                  warning ('detections with coordinates outside polygon(s) were dropped')
              }
          }
          else if (detector(traps)=='transect') {
              captTrap <- xyontransect(captures[,4:5], traps, tol)
              if (any(captTrap==0)) {
                  captTrap <- captTrap[captTrap>0]
                  captures <- captures[captTrap>0,]
                  warning ('detections with coordinates not on any transect were dropped')
              }
          }
          else {
              trapID    <- interaction(traps$x, traps$y)
              captTrap  <- match(interaction(captures[,4], captures[,5]), trapID)
              if (any(is.na(captTrap))) stop ('Failed to match some capture locations to detector sites')
          }
        }
        else {
          if (detector(traps) %in% c('polygon', 'transect'))
              stop ('use fmt XY to input detections from polygons or transects')
          captTrap <- match(captures[,4], row.names(traps))
        }

        #  if (bysession & ( length(levels(session)) > 1)) {
        #    captures[,2] <- interaction(session, captures[,2], drop = TRUE)
        #  }

        nocc      <- max(abs(captures[,3]))
        nocc      <- ifelse (is.null(noccasions), nocc, noccasions)

        if (is.null(detector(traps)))
            stop ("require a detector type e.g. detector(traps) <- 'multi'")
        if (is.null(cutval) && detector(traps)=='signal')
            stop ("Missing cutval (signal threshold) for signal data")

        wout <- NULL
        ID   <- NULL

        uniqueID <- unique(captures[,2])
        captID <- as.numeric(factor(captures[,2], levels=uniqueID))
        nID    <- length(uniqueID)

        dim3 <- detector(traps) %in% c('proximity', 'quadratbinary','signal', 'count',
            'quadratcount','polygon','transect')
        if (dim3) {
            w <- array (0, dim=c(nID, nocc, ndetector(traps)))
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
        else {
            w     <- matrix(nr=nID, nc=nocc)
            w[,]  <- 0

            ## adjusted 2009 08 13 to ensure first occurrence selected when more than one per occasion
            indices <- cbind(captID, abs(captures[,3]))
            values <- captTrap * sign(captures[,3])
            ord <- order (captID, 1:length(captID), decreasing=TRUE)
            w[indices[ord,, drop=F]] <- values[ord]    ## drop=F is critical to ensure retains dim2
        }

        wout <- abind(wout, w, along=1)
        dimnames(wout)[[1]] <- uniqueID
        dimnames(wout)[[2]] <- 1:nocc

        ## check added 2010 05 01
        if (any(is.na(wout))) {
            rw <- row(wout)[is.na(wout)]
            if (dim3)
                print (wout[rw,,,drop=F])
            else
                print (wout[rw,,drop=F])
            stop ('missing values not allowed')
        }

        ## optional row sort 2009 09 26, tweaked 2010 05 01
        if (sortrows) {
            if (suppressWarnings( all(!is.na(as.numeric(uniqueID)))))
                roworder <- order(as.numeric(uniqueID))
            else
                roworder <- order (uniqueID)
            if (length(dim(wout))==3)
                wout[,,] <- wout[roworder,,]
            else
                wout[,] <- wout[roworder,]
            dimnames(wout)[[1]] <- dimnames(wout)[[1]][roworder]
        }

        ## code to input permanent individual covariates if these are present
        zi <- NULL
        startcol <- ifelse (fmt=='trapID', 5, 6)
        if (detector(traps) %in% c('signal')) startcol <- startcol+1
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
                temp[,j] <- factor(tempj2)
                rownames(temp) <- names(tempj)          ## 2010 02 26
            }
            if (is.null(covnames)) names(temp) <- names(zi)
            else {
                if (ncol(temp) != length(covnames)) stop('Number of covariate names does not match')
                names(temp) <- covnames
            }
            if (sortrows) temp <- temp[roworder,,drop = FALSE]   ## added 'FALSE' 2010 02 24
            attr(wout,'covariates') <- temp
        }
        else attr(wout,'covariates') <- data.frame()

        class (wout) <- 'capthist'
        traps(wout) <- traps
        if (detector(traps) %in% c('polygon','transect')) {
            xy <- captures[order(captures[,2],captures[,3],captTrap),4:5]
            names(xy) <- c('x','y')
            attr(wout,'detectedXY') <- xy
        }
        if (detector(traps) == 'signal') {
            if (is.null(cutval)) stop ('Missing value for signal threshold')
            if (fmt=='XY')
                attr(wout, 'signal') <- captures[,6]
            else
                attr(wout, 'signal') <- captures[,5]
            attr(wout, 'cutval')   <- cutval
        }

        session(wout)  <- as.character(captures[1,1])
        wout
    }   ## end of single-session call
}
############################################################################################
