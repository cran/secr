## 2011-12-01, 2011-12-07, 2012-06-08, 2012-11-04
## join returns single-session object from list of inputs

# should new trap ID be numeric , character or factor???
# 2015-01-31 join() bug fix could fail with exclusivedetector types

join <- function (object, remove.dupl.sites = TRUE, tol = 0.001) {
    onesession <- function (sess) {
        ## previous is cumulative number of preceding occasions
        CH <- object[[sess]]
        previous <- before[sess]
        newID <- animalID(CH)
        occ <- occasion(CH)
        newocc <- occ + previous
        newtrap <- trap(CH)
        df <- data.frame(newID = newID, newocc = newocc, newtrap = newtrap,
                         alive = alive(CH), sess = rep(sess, length(occ)),
                         stringsAsFactors = FALSE)
        if (!is.null(xy(CH))) {
            df$x <- xy(CH)$x
            df$y <- xy(CH)$y
        }
        if (!is.null(signal(CH)))
            df$signal <- signal(CH)
        df
    }
    #    if (!ms(object) | !inherits(object, 'capthist'))   
    if (!ms(object) | any(sapply(object, class) != 'capthist'))
        stop("requires multi-session capthist object or list of single-session capthist")
    
    outputdetector <- detector(traps(object)[[1]])
    nsession <- length(object)

    nocc <- sapply(object, ncol)
    names(nocc) <- NULL
    before <- c(0, cumsum(nocc)[-nsession])
    df <- vector(mode= 'list', length = nsession)
    for (sess in 1:nsession) {
        df[[sess]] <- onesession (sess)
    }
    df <- do.call(rbind, df)

    n <- length(unique(df$newID))
    condition.usage <- function (trp, i, nocc) {
        us <- matrix(0, nrow=nrow(trp), ncol=sum(nocc))
        s1 <- c(1, cumsum(nocc)+1)[i]
        s2 <- cumsum(nocc)[i]
        if (is.null(usage(trp)))
            us[,s1:s2] <- 1
        else
            us[,s1:s2] <- usage(trp)
        usage(trp) <- us
        trp
    }

    ## resolve traps
    temptrp <- traps(object)
    for (i in 1:length(temptrp))
        temptrp[[i]] <- condition.usage(temptrp[[i]], i, nocc)

    alltemptrp <- do.call(rbind, c(temptrp, renumber = FALSE))
    newtraps <- alltemptrp
    df$newtrap <- paste(df$newtrap,df$sess, sep=".")

    # ensure retain all occasions
    df$newocc <- factor(df$newocc, levels = 1:sum(nocc))

    if (outputdetector %in% .localstuff$exclusivedetectors) {
        alivesign <- df$alive*2 - 1
        tempnew <- matrix(0, nrow = n, ncol = sum(nocc))
        dimnames(tempnew) <- list(unique(df$newID), 1:sum(nocc))
        
        ## bug fix 2015-01-31
        ## tempnew[cbind(df$newID, df$newocc)] <- as.numeric(df$newtrap) * alivesign
        trapindex <- match(df$newtrap, unique(df$newtrap))        
        tempnew[cbind(df$newID, df$newocc)] <- trapindex * alivesign
    }
    else {
        df$newtrap <- factor(df$newtrap, levels=rownames(newtraps))
        tempnew <- table(df$newID, df$newocc, df$newtrap)
        alivesign <- tapply(df$alive, list(df$newID,df$newocc,df$newtrap),all)
        alivesign[is.na(alivesign)] <- TRUE
        alivesign <- alivesign * 2 - 1
        tempnew <- tempnew * alivesign
    }
    class(tempnew) <- 'capthist'
    traps(tempnew) <- newtraps
    session(tempnew) <- 1
    neworder <- order (df$newocc, df$newID, df$newtrap)

    if (!is.null(df$x))
        xy(tempnew) <- df[neworder,c('x','y'), drop = FALSE]

    if (!is.null(df$signal))
        signal(tempnew) <- df[neworder,'signal']

    if (!is.null(covariates(object))) {
        tempcov <- do.call(rbind, covariates(object))
        IDcov <- unlist(lapply(object,rownames))
        ## use first match
        tempcov <- tempcov[match(rownames(tempnew), IDcov),,drop = FALSE]
        rownames(tempcov) <- rownames(tempnew)
        covariates(tempnew) <- tempcov
    }

    if (remove.dupl.sites)
        tempnew <- reduce(tempnew, span=tol, dropunused = FALSE, verify = FALSE)
    attr(tempnew, 'interval') <- unlist(sapply(nocc, function(x) c(1,rep(0,x-1))))[-1]

    tempnew

}

RMarkInput <- function (object, grouped = FALSE, covariates = TRUE) {
    if (!inherits(object, "capthist"))
        stop ("requires single-session capthist object")
    if (ms(object))
        stop ("requires single-session capthist object - use 'join'")
    if (length(dim(object)) != 2)
        CH0 <- reduce(object, outputdetector = 'multi')
    else
        CH0 <- object
    ntimes <- ncol(object)
    alive <- apply(object,1,function(x) all(x>=0))
    CH <- pmin(abs(CH0),1)

    if (is.logical(covariates)) {
        if (covariates) {
            if (is.null(covariates(CH0)))
                stop("no covariates in object")
            covnames <- names(covariates(CH0))
        }
        else
            covnames <- ""
    }
    else {
        covnames <- covariates
        if (is.character(covariates)) {
            if (is.null(covariates(CH0)))
                stop("no covariates in object")
        }
        found <- covnames %in% names(covariates(CH0))
        if (any(!found)) {
            stop(paste(covnames[!found], collapse=','), " not in covariates(object)")
        }
    }

    if (any(covnames != "")) {
        if (grouped) {
            warning("'grouped' is incompatible with individual covariates and will be ignored")
            grouped <- FALSE
        }
    }

    if (grouped)   ## bug fix 2012-07-04
        CH <- cbind(CH, alive) ## add single-digit code as last column
    CH <- data.frame(ch = apply(CH, 1, paste, collapse=''),
        stringsAsFactors = FALSE)

    if (grouped) {
        temp <- table(CH$ch)
        alive <- as.numeric(substring(names(temp),ntimes+1,ntimes+1))
        CH <- data.frame(ch = substring(names(temp),1,ntimes),
                         freq = as.numeric(temp),
                         stringsAsFactors = FALSE)
        CH$freq <- ifelse(alive, CH$freq, -CH$freq)
        CH <- CH[order(CH$ch, decreasing = TRUE),]
        row.names(CH) <- 1:nrow(CH)
    }
    else {
        CH$freq <- ifelse(alive,1,-1)
        if (any(covnames != "")) {
            CH[,covnames] <- covariates(CH0)[,covnames]
        }
    }
    attr(CH, "interval") <- attr(object, "interval")
    if (is.null(attr(CH,"interval")))
        attr(CH, "interval") <- rep(0,ntimes-1)
    CH
}

# temp <- join(ovenCH)
# ovenint <- c(rep(0,8),rep(c(1,rep(0,9)),4))
# secr.fit(temp, details=list(interval = ovenint))

# plot(join(ovenCH), tracks=T)

unjoin <- function (object, interval, ...) {
    if (missing(interval) & is.null(attr(object,'interval')))
        stop ("requires 'interval' to define sessions")
    if (missing(interval) )
        interval <- attr(object,"interval")
    session <- c(0,cumsum(interval>0))+1
    nsess <- max(session)
    if (nsess<2) {
        warning ("interval define only one session")
        return(object)
    }
    newobj <- vector(mode='list', length=nsess)
    for (sess in 1:nsess) {
        newobj[[sess]] <- subset(object, occasions = (session==sess), ...)
    }
    class (newobj) <- c('list','capthist')
    session(newobj) <- 1:nsess
    return(newobj)
}

unRMarkInput <- function(df, covariates = TRUE) {
    if (!is.data.frame(df))
        stop("requires dataframe input")
    if (!all(c('ch','freq') %in% names(df)))
        stop ("ch and freq are required fields")
    nocc <- nchar(df$ch)
    if (length(unique(nocc))>1)
        stop ("ch must be a constant-length string of 0s and 1s")
    nocc <- nocc[1]
    freq <- df$freq
    alive <- sign(freq)
    freq <- abs(freq)
    freq <- rep(1:nrow(df), freq)
    alive <- alive[freq]
    ch <- df$ch[freq]
    CH <- matrix(as.numeric(unlist(sapply(ch, strsplit, ''))), byrow = TRUE, ncol = nocc)
    class(CH) <- 'capthist'
    # allow deads
    last <- function(x) which.max(cumsum(x))
    CH[cbind(1:nrow(CH), apply(CH,1,last))] <- alive
    # transfer covariates if present
    if (ncol(df)>2) {
        if (is.logical(covariates)) {
            if (covariates)
                covnames <- names(df)[-match(c('ch','freq'),names(df))]
            else
                covnames <- ""
        }
        else {
            covnames <- covariates[covariates %in% names(df)]
        }

        if (any(covnames != ""))
        covariates(CH) <- df[freq, covnames, drop = FALSE]
    }
    CH
}
