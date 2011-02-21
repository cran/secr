############################################################################################
## package 'secr'
## read.traps.R
## last changed 2011-02-20
## Read detector locations from text file in DENSITY format
############################################################################################

renamepolyrows <- function (tr) {
    # update vertex labels for polygon and transect detectors
    if (detector(tr) %in% c('polygon','polygonX','transect','transectX')) {
        tr.freq <- table(polyID(tr))
        vertex <- unlist(sapply(tr.freq, function(x) 1:x))
        rownames(tr) <-  paste(polyID(tr),vertex,sep='.')
    }
    tr
}

read.traps <- function (file = NULL, data = NULL, detector = 'multi', covnames = NULL, ...)
## possibly add sortrows argument, but risk breaking covariates, polygon etc.
{

    # count.fields(file, sep = "", quote = "\"'", skip = 0, blank.lines.skip = TRUE,
    #     comment.char = "#")
    closepoly <- function (xyi) {    ## close polygon
        if ((tail(xyi$x,1) != xyi$x[1]) | (tail(xyi$y,1) != xyi$y[1])) {
            newrow <- xyi[1,,drop=F]
            xyi <- rbind(xyi,newrow)
        }
        xyi
    }

    if (!( detector %in% .localstuff$validdetectors ))
        stop ('invalid detector type')
    if (is.null(file) & is.null(data))
        stop ("requires 'file' or 'data'")
    ## file input
    if (!is.null(file)) {
        nfld <- count.fields (file, ...)
        if (min(nfld) < 3)
            stop ("requires 3 fields (detectorID, x, y)")
        if (min(nfld) == 3) colcl <- c('character',NA,NA)
        else colcl <- c('character',NA,NA,'character')
        if (detector %in% c('polygon', 'polygonX')) {
            temp <- read.table (file, row.names=NULL, as.is=T, colClasses=colcl, ...)
            names(temp)[1:3] <- c('polyID','x','y')

            tempID <- levels(factor(temp$polyID))
            tempindex <- match (tempID, as.character(temp$polyID))
            temp1 <- split (temp[, 1:3, drop=FALSE], temp$polyID)
            temp1 <- abind(lapply(temp1, closepoly), along=1, force.array=F)

            traps <- temp1[, c('x','y'),drop=FALSE]
            class(traps)    <- c('traps', 'data.frame')
            polyID(traps) <- factor(temp1[,1])
            temp <- temp[,-1,drop=FALSE]   ## discard polyID so usage >= col3
        }
        else if (detector %in% c('transect', 'transectX')) {
            temp <- read.table (file, row.names=NULL, as.is=T, colClasses=colcl, ...)
            names(temp)[1:3] <- c('transectID','x','y')
            tempID <- levels(factor(temp$transectID))
            tempindex <- match (tempID, as.character(temp$transectID))
            traps <- temp[, c('x','y'),drop=FALSE]
            class(traps)    <- c('traps', 'data.frame')
            transectID(traps) <- factor(temp[,1])
            temp <- temp[,-1,drop=FALSE]   ## discard transectID so usage >= col3
        }
        else {
            temp <- read.table (file, row.names=1, as.is=T, colClasses=colcl, ...)
            traps <- temp[,1:2,drop=FALSE]
            class(traps)    <- c('traps', 'data.frame')
        }
        dimnames(traps)[[2]] <- c('x','y')
    }
    ## dataframe input
    else {
        ## close polygons; tempindex used later to match covariates & usage
        if ('polyID' %in% names(data)) {
            temp <- split (data[,c('x','y','polyID'),drop=FALSE], data$polyID)
            data <- as.data.frame(abind(lapply(temp, closepoly), along=1), force.array=F)
        }
        traps <- data[,c('x','y'), drop=FALSE]
        class(traps)    <- c('traps', 'data.frame')
        if ('polyID' %in% names(data))
            polyID(traps) <- factor(data$polyID)
        if ('transectID' %in% names(data))
            transectID(traps) <- factor(data$transectID)
    }

    detector(traps) <- detector
    traps <- renamepolyrows(traps)
    usage(traps)      <- NULL
    covariates(traps) <- NULL
    if (!is.null(file)) {
        if (ncol(temp)>2) {
            if (ncol(temp)>3)
                temp2 <- apply(temp[,3:ncol(temp),drop=FALSE], 1, paste, collapse='')
            else
                temp2 <- temp[,3]
            splitfield <- matrix(unlist(strsplit(as.character(temp2),'/')),
                byrow = T, nrow = length(temp2))                            # before '/'
            used <- gsub(' ','', splitfield[,1,drop=FALSE])                 # remove blanks
            used <- gsub('//t','', used)                                    # remove tabs
            used <- gsub(',','', used)                                      # remove commas
            nocc <- max(nchar(used))

            if (nocc>0) {
                if (detector %in% c('polygon','polygonX','transect','transectX')) {
                    used <- used[tempindex]
                    nocc <- max(nchar(used))   ## recompute
                    if (any(nchar(used) != nocc))
                        stop ("'usage' fields suggest varying number of occasions")
                    usage(traps) <- matrix(unlist(strsplit(used,split='')),
                                           byrow=T, ncol = nocc)>0
                    dimnames(attr(traps, 'usage')) <- list( tempID, 1:nocc)
                }
                else {
                    if (any(nchar(used) != nocc))
                        stop ("'usage' fields suggest varying number of occasions")
                    usage(traps) <- matrix(unlist(strsplit(used,split='')),
                                           byrow=T, ncol = nocc)>0
                    dimnames(attr(traps, 'usage')) <- list(dimnames(traps)[[1]],
                                           1:nocc)
                }
            }

            if (ncol(splitfield)>1) {
                if (detector %in% c('polygon','polygonX','transect','transectX')) {
                    splitfield <- splitfield[tempindex,, drop = FALSE]
                    rownam <- tempID
                }
                else {
                    rownam <- rownames(traps)
                }

                tempcon <- textConnection(splitfield[,2, drop = FALSE])
                covariates(traps) <- data.frame(read.table(tempcon))
                rownames(covariates(traps)) <- rownam
                close(tempcon)

                if (!is.null(covnames)) {
                     ncc <- ncol(covariates(traps))
                     if (length(covnames) != ncc)
                         stop("number of covariate names does",
                              "not match number of columns", ncc )
                     names (covariates(traps)) <- covnames
                 }
            }
        }
    }
    attr(traps,'spacex') <- min(dist(unique(traps$x)))
    attr(traps,'spacey') <- min(dist(unique(traps$y)))
    spacing(traps) <- spacing(traps)   ## !!
    traps
}
############################################################################################
