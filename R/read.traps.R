############################################################################################
## package 'secr'
## read.traps.R
## 2012-12-18 binary.usage
## 2015-05-04 multi-file input
## 2015-10-12 markocc
## 2020-08-26 allow usage and covariates with data and xls input

## Read detector locations from text file in DENSITY format
############################################################################################

renamepolyrows <- function (tr) {
    # update vertex labels for polygon and transect detectors
    if (all(detector(tr) %in% .localstuff$polydetectors)) {
        tr.freq <- table(polyID(tr))
        vertex <- unlist(sapply(tr.freq, function(x) 1:x))
        rownames(tr) <-  paste(polyID(tr),vertex,sep='.')
    }
    tr
}

read.traps <- function (file = NULL, data = NULL, detector = 'multi', covnames = NULL,
    binary.usage = TRUE, markocc = NULL, trapID = NULL, ...)
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
    
    if (!all( detector %in% .localstuff$validdetectors ))
        stop ("invalid detector type")
    if (is.null(file) & is.null(data))
        stop ("requires 'file' or 'data'")
    
    if (length(file) > 1) {
        if (is.null(markocc))
            out <- mapply(read.traps, file=file, MoreArgs = list(detector = detector,
                covnames = covnames, binary.usage = binary.usage, ...), SIMPLIFY = FALSE)
        else {
            if (!inherits(markocc, 'list'))
                markocc <- list(markocc)
            out <- mapply(read.traps, file=file, markocc = markocc, MoreArgs = list(detector = detector,
                covnames = covnames, binary.usage = binary.usage, ...), SIMPLIFY = FALSE)
        }
        class (out) <- c('traps', 'list')
        out
    }
    else {
        
        ## file input
        if (!is.null(file)) {
            if (tolower(tools::file_ext(file)) %in% c("xls", "xlsx")) {
                if (!requireNamespace("readxl", quietly = TRUE))
                    stop("package readxl is required for input from Excel spreadsheets")
                data <- readxl::read_excel(file, ...)
                data <- data.frame(data)  # not a tibble
                nam <- names(data)
                if (!all(c('x','y') %in% nam)) {
                    names(data)[1:3] <- c('trapID','x', 'y')
                    warning ("assuming columns 1-3 are trapID, x, y")
                    trapID <- "trapID"
                }
                else {    
                    nam <- nam[!(nam %in% c('x','y'))]
                    if (length(nam)<1) warning("no trapID column")
                    if (is.null(trapID)) {
                        trapID <- nam[1]
                    }
                    else {
                        if (!trapID %in% nam)
                            stop("trapID is not a column name")
                    }
                }
                rownames(data) <- data[,trapID]
                out <- read.traps(data = data, markocc = markocc, detector = detector,
                    covnames = covnames, binary.usage = binary.usage)
                return(out)
            }
            else {
                nfld <- count.fields (file, ...)
                if (min(nfld) < 3)
                    stop ("requires 3 fields (detectorID, x, y)")
                if (min(nfld) == 3) colcl <- c('character',NA,NA)
                else colcl <- c('character',NA,NA,'character')
                if (all(detector %in% c('polygon', 'polygonX'))) {
                    temp <- read.table (file, row.names=NULL, as.is=T,
                        colClasses = colcl, ...)
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
                else if (all(detector %in% c('transect', 'transectX'))) {
                    temp <- read.table (file, row.names=NULL, as.is=T,
                        colClasses = colcl, ...)
                    names(temp)[1:3] <- c('transectID','x','y')
                    tempID <- levels(factor(temp$transectID))
                    tempindex <- match (tempID, as.character(temp$transectID))
                    traps <- temp[, c('x','y'),drop=FALSE]
                    class(traps)    <- c('traps', 'data.frame')
                    transectID(traps) <- factor(temp[,1])
                    temp <- temp[,-1,drop=FALSE]   ## discard transectID so usage >= col3
                }
                else {
                    temp <- read.table (file, row.names = 1, as.is = T,
                        colClasses = colcl, ...)
                    traps <- temp[,1:2,drop=FALSE]
                    class(traps)    <- c('traps', 'data.frame')
                }
                dimnames(traps)[[2]] <- c('x','y')
            }
        }
        ## dataframe input
        else {
            ## close polygons; tempindex used later to match covariates & usage
            if ('polyID' %in% names(data)) {
                temp <- split (data[,c('x','y','polyID'),drop=FALSE], data$polyID)
                data <- as.data.frame(abind(lapply(temp, closepoly), along = 1),
                    force.array = FALSE)
            }
            else {
                if (!is.null(trapID))
                    rownames(data) <- data[,trapID]
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
        
        ## 2020-08-26
        temp2 <- NULL
        if (!is.null(file) && ncol(temp)>2) {
            if (ncol(temp)>3)
                temp2 <- apply(temp[,3:ncol(temp),drop=FALSE], 1, paste, collapse=' ')
            else
                temp2 <- temp[,3]
        }
        else if (!is.null(data) && ncol(data)>3) {
            temp2 <- apply(data[,-(1:3), drop = FALSE],1,paste, collapse = ' ')
        }
        if  (!is.null(temp2)) {               
            
            splitfield <- matrix(unlist(strsplit(as.character(temp2),'/')),
                byrow = T, nrow = length(temp2))                            # before '/'
            
            if (binary.usage) {
                used <- gsub(' ','', splitfield[,1,drop=FALSE])                 # remove blanks
                used <- gsub('//t','', used)                                    # remove tabs
                used <- gsub(',','', used)                                      # remove commas
                nocc <- max(nchar(used))
                if (nocc>0) {
                    if (all(detector %in% .localstuff$polydetectors)) {
                        
                        stop("usage input not available for polygon detectors in 4.3.1")
                        used <- used[tempindex]
                        nocc <- max(nchar(used))   ## recompute
                        
                        if (any(nchar(used) != nocc))
                            stop ("'usage' fields suggest varying number ",
                                "of occasions")
                        usge <- (matrix(unlist(strsplit(used,split='')),
                            byrow = TRUE, ncol = nocc)>0) * 1
                        
                        dimnames(usge) <- list( tempID, 1:nocc)
                        usage(traps) <- usge
                    }
                    else {
                        if (any(nchar(used) != nocc))
                            stop ("'usage' fields suggest varying number ",
                                "of occasions")
                        usge <- (matrix(unlist(strsplit(used,split='')),
                            byrow = TRUE, ncol = nocc)>0) * 1
                        dimnames(usge) <- list(dimnames(traps)[[1]], 1:nocc)
                        usage(traps) <- usge
                    }
                }
            }
            else {
                ## tentative 2012-12-18
                tempcon <- textConnection(splitfield[,1, drop = FALSE])
                usge <- as.matrix(read.table(tempcon))
                close(tempcon)
                nocc <- ncol(usge)
                if (all(detector %in% .localstuff$polydetectors)) {
                    ## 2014-08-23 bug fix
                    ## usge <- usge[tempindex,]
                    usge <- usge[tempindex,,drop = FALSE]
                    dimnames(usge) <- list( tempID, 1:nocc)
                }
                else {
                    dimnames(usge) <- list(dimnames(traps)[[1]], 1:nocc)
                }
                usage(traps) <- usge
            }
            if (length(detector)>1 & length(detector) != nocc)
                stop ("differing number of occasions in detector and usage")
            if (ncol(splitfield)>1) {
                if (all(detector %in% .localstuff$polydetectors)) {
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
                        stop ("number of covariate names does ",
                            "not match number of columns ", ncc )
                    names (covariates(traps)) <- covnames
                }
                
            }
        }
        ux <- unique(traps$x)
        uy <- unique(traps$y)
        if (!is.numeric(ux) | !is.numeric(uy))
            stop("non-numeric coordinates detected")
        attr(traps,'spacex') <- ifelse (length(ux)>1, min(dist(ux)), NA)
        attr(traps,'spacey') <- ifelse (length(uy)>1, min(dist(uy)), NA)
        spacing(traps) <- spacing(traps)   ## !!
        markocc(traps) <- markocc
        traps
    }
}
############################################################################################
