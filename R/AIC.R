###############################################################################
## package 'secr'
## AIC.R
## 2018-01-10
###############################################################################

oneline.secr <- function (secr, k) {
    
    if (ms(secr$capthist)) {
        n <- sum(sapply (secr$capthist, nrow))
        ncapt <- sum(sapply (secr$capthist, function (x) sum(abs(x>0))))
    }
    else  {
        n  <- nrow(secr$capthist)     # number caught
        ncapt <- sum( abs( secr$capthist)>0)
    }
    
    ## 2015-03-31, 2017-02-10
    Npar <- nparameters(secr)   ## see utility.R
    AICval <- 2*secr$fit$value + k *Npar
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

AIC.secr <- function (object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c('AICc','AIC')) {
    allargs <- list(...)
    modelnames <- (c ( as.character(match.call(expand.dots=FALSE)$object),
                       as.character(match.call(expand.dots=FALSE)$...) ))
    allargs <- secrlist(object, allargs)
    names(allargs) <- modelnames
    AIC(allargs, sort=sort, k=k, dmax=dmax, criterion=criterion)
}
############################################################################################
############################################################################################

AIC.secrlist <- function (object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c('AICc','AIC')) {
    
    if (k != 2)
        warning ("k != 2 and AIC.secr output may be mis-labelled")
    
    if (length(list(...)) > 0)
        warning ("... argument ignored in 'AIC.secrlist'")
    
    if (length(object) > 1) {
        ## check added 2013-10-14
        hcovs <- sapply(object, function(x) if (is.null(x$hcov)) '' else x$hcov)
        if (length(unique(hcovs)) > 1)
            stop ("AIC invalid when models use different hcov")
    }
    
    criterion <- match.arg(criterion)
    modelnames <- names(object)
    allargs <- object
    if (any(sapply(allargs,class) != 'secr'))
        stop ("components of 'object' must be 'secr' objects")
    
    output <- data.frame(t(sapply(allargs, oneline.secr, k = k)), stringsAsFactors=F)
    for (i in 3:6)
        output[,i] <- as.numeric(output[,i])
    
    output$delta <- output[,criterion] - min(output[,criterion])
    OK <- abs(output$delta) < abs(dmax)
    sumdelta <- sum(exp(-output$delta[OK]/2))
    output$wt <- ifelse ( OK, round(exp(-output$delta/2) / sumdelta,4), 0)
    row.names(output) <- modelnames
    if (sort) output <- output [order(output[,criterion]),]
    names(output)[7] <- paste('d',criterion,sep='')
    names(output)[8] <- paste(criterion,'wt',sep='')
    if (nrow(output)==1) { output[,8] <- NULL; output[,7] <- NULL}
    
    output
}
############################################################################################
