############################################################################################
## package 'secr'
## scoretest.R
## score test for secr models
## last changed 2009 01 19, 2009 07 10, 2009 08 15, 2009 10 07 (q)
## 2009 08 25 fixed AIC small sample adjustment
## 2009 08 25 block re-zeroing of AICc vs constrained model
## requires nlme
############################################################################################
## source ('d:\\density secr 1.1\\secr\\R\\scoretest.R')
############################################################################################

prepare <- function (secr, newmodel) {

    ## Prepare detection design matrices and lookup
    ## secr provides raw data etc.

    capthist <- secr$capthist
    mask     <- secr$mask
    timecov  <- secr$timecov
    sessioncov <- secr$sessioncov
    groups   <- secr$groups
    dframe   <- secr$dframe  ## never used?
    CL       <- secr$CL
    detectfn <- secr$detectfn
    link     <- secr$link
    details  <- secr$details

    sessionlevels <- session(capthist)
    if (is.null(sessionlevels)) sessionlevels <- '1'
    grps  <- group.levels(capthist, groups)

    design <- secr.design.MS (capthist, newmodel, timecov, sessioncov, groups, dframe)
    design0 <- secr.design.MS (capthist, newmodel, timecov, sessioncov, groups, dframe, naive=T, bygroup=!CL)

    if (CL) D.designmatrix <- NULL
    else {
        temp <- D.designdata ( mask, newmodel$D, grps, sessionlevels, sessioncov)
        D.designmatrix <- model.matrix(newmodel$D, temp)
        attr(D.designmatrix, 'dimD') <- attr(temp, 'dimD')
    }

    #############################
    # Parameter mapping (general)
    #############################

    D.modelled <- !CL                                           ## & is.null(fixed$D)
    np <- sapply(design$designMatrices, ncol)
    if (D.modelled) np <-  c(D = ncol(D.designmatrix), np)

    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)
    if (!D.modelled) parindx$D <- NULL

## DOES THIS DEAL WITH SESSION COVARIATES OF DENSITY? 2009 08 16

    ################
    # Variable names
    ################

    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    if (!CL) betanames <- c(paste('D', colnames(D.designmatrix), sep='.'), betanames)
    betanames <- sub('..(Intercept))','',betanames)

    if (detector(traps(capthist)) %in% c('polygon', 'polygonX', 'transect', 'transectX',
                                         'signal','cue'))
        savedlogmultinomial <- 0
    else
        savedlogmultinomial <- logmultinom(capthist, group.factor(capthist, groups))

    ################
    list (
       link      = link,
       parindx   = parindx,
       capthist  = capthist,
       mask      = mask,
       CL        = CL,
       detectfn  = detectfn,

       D.design  = D.designmatrix,
       design    = design,
       design0   = design0,

       groups    = groups,
       betanames = betanames,
       details   = details,
       logmult   = savedlogmultinomial
     )

 }   # end of 'prepare'
############################################################################################

mapbeta <- function (parindx0, parindx1, beta0, betaindex)

## Extend beta vector from simple model (beta0) to a more complex (i.e. general)
## model, inserting neutral values (zero) as required.
## For each real parameter, a 1:1 match is assumed between
## beta values until all beta values from the simpler model are
## used up. THIS ASSUMPTION MAY NOT BE JUSTIFIED.
## betaindex is a user-controlled alternative.

{
    ## list of zeroed vectors, one per real parameter
    beta1 <- lapply(parindx1, function (x) {x[]<-0; x})

    if (!is.null(betaindex)) {
        beta1 <- unlist(beta1)
        if (sum(betaindex>0) != length(beta0))
            stop ("invalid 'betaindex'")
        beta1[betaindex] <- beta0
        beta1
    }
    else {
        indx <- lapply(parindx0, function(x) x-x[1]+1)
        for (j in 1:length(beta1))
          beta1[[j]][indx[[j]]] <- beta0[parindx0[[j]]]
        unlist(beta1)
    }
}
############################################################################################

score.test <- function (secr, ..., betaindex = NULL, trace = FALSE) {

    if (!inherits(secr, 'secr'))
        stop ("only for 'secr' objects")

    models <- list(...)

    if (length(models)>1) { ##  ) {
        # apply to each component of 'model' list
        score.list <- lapply (models, score.test, secr=secr)
        class(score.list) <- c('score.list', class(score.list))
        score.list
    }
    else {
        if (is.list(models) & !('formula' %in% class(models[[1]])) )
            model <- models[[1]]
        else
            model <- models

        if (inherits (model, 'secr')) model <- model$model
        # model may be incompletely specified
        model <- stdform (model)

        model <- replace (secr$model, names(model), model)
        if (secr$CL) model$D <- NULL
        if (secr$detectfn!=1) model$z <- NULL

        #########################################################
        # check components on which model differs from secr$model

        if (length(secr$model) != length(model))
            stop ("models differ in number of real parameters")
        if (any(names(secr$model) != names(model)))
            stop ("real parameters do not match")
        differ <- logical(length(model))
        for (j in 1:length(model)) differ[j] <- secr$model[[j]] != model[[j]]
        if (sum(differ)==0)
            stop ("models are identical")

        #########################################################
        # construct essential parts of new secr model

        newsecr <- prepare (secr, model)
        newsecr$details <- replace (secr$details, 'trace', trace)  ## override

        beta0 <- secr$fit$par
        beta1 <- mapbeta(secr$parindx, newsecr$parindx, beta0, betaindex)
        names(beta1) <- newsecr$betanames

        if (trace) {
            ## .localstuff is environment for temp variables in namespace
            .localstuff$iter <- 0
            cat('Eval     logLik', formatC(newsecr$betanames, format='f', width=11), '\n')
        }

        loglikfn <- function (beta, design) {
           -secr.loglikfn(
               beta     = beta,
               parindx  = design$parindx,
               link     = design$link,
               fixed    = list(),
               designD  = design$D.design,
               design   = design$design,
               design0  = design$design0,
               capthist = design$capthist,
               mask     = design$mask,
               detectfn = design$detectfn,
               CL       = design$CL,
               groups   = design$groups,
               details  = design$details,
               logmult  = design$logmult,
               dig      = 4,
               betaw    = 11)
        }

        # cat('logLik',loglikfn(beta1, design=newsecr), '\n')   ## testing
        # cat('logmult', newsecr$logmult, '\n')

        ## fdHess is from package nlme
        require(nlme)
        grad.Hess <- fdHess(beta1, fun = loglikfn, design = newsecr,
                         .relStep = 0.001,
                         minAbsPar=0.1)

        u.star <- grad.Hess$gradient
        i.star <- -grad.Hess$Hessian

## debug 2009 08 25
## library(numDeriv)
## u.star <- grad(func = loglikfn, x = beta1, method="Richardson", list(eps=1e-3, d=0.1,
##      zero.tol=0.001, r=4, v=2), design = newsecr)
## i.star <- -hessian(func = loglikfn, x = beta1, method="Richardson", list(eps=1e-3, d=0.1,
##      zero.tol=0.001, r=4, v=2), design = newsecr)
## print(grad.Hess$mean)
## print (u.star)
## print(i.star)

        score <- try (solve(i.star), silent = TRUE)
        if (inherits(score, "try-error")) {
            warning ("could not invert information matrix")
            score <- NA
        }
        else  score <- u.star %*% score %*% u.star

        ## score statistic cf -2(logLik1 - logLik0)
        statistic <- c('X-squared' = as.numeric(score))

        ## number caught, for AICc
        n <- ifelse (is.list(secr$capthist),
            sum(sapply (secr$capthist, nrow)),
            nrow(secr$capthist) )

        np0 <- length(beta0)
        np1 <- length(beta1)
        np  <- c(np0, np1)
        parameter <- c(df = np1-np0)

        H0    <- model.string(secr$model)
        H1    <- model.string(model)
        call0 <-  paste('[', format(secr$call), ']', sep = '')
        AIC   <- c(0, -statistic + 2 * parameter)
        AICc  <- ifelse ((n-np-1)>0,
            AIC + 2 * np * (np+1) / (n-np-1),
            NA)
##        AICc  <- AICc - AICc[1]  # relative to simpler (constrained) model
        names(AIC) <- c('AIC.0','AIC.1')
        names(AICc) <- c('AICc.0','AICc.1')

        temp <- list(
             statistic = statistic,
             parameter = parameter,
             p.value = 1 - pchisq(score, parameter),
             method = 'Score test for two SECR models',
             data.name = paste(H0, 'vs', H1, call0),
             np0 = np0,
             H0 = H0,
             H1 = H1,
             H1.beta = beta1,
             AIC = AIC,
             AICc = AICc
        )
        class(temp) <- c('score.test', 'htest')
        temp

    }
}
############################################################################################

score.table <- function (object, ..., sort = TRUE, dmax = 10) {

    if ('score.test' %in% class(object)) object <- list(object)
    if (length(match.call(expand.dots=F)$...)==0) score.list <- object
    else score.list <- c(object, list(...))
    if ( any( sapply(score.list, function(x) !('score.test' %in% class(x)) )) )
        stop ("arguments must be 'score.test' objects")

    output <- data.frame(model= sapply (score.list, function(x) x$H1))
    output$X.squared <- sapply(score.list, function(x) x$statistic)
    output$df <- sapply(score.list, function(x) x$parameter)
    output$p.value <- round(sapply(score.list, function(x) x$p.value),6)
    output$npar <- score.list[[1]]$np0 + sapply(score.list, function(x) x$parameter)
    output$AIC <- sapply(score.list, function(x) x$AIC[2])
    output$AICc <- sapply(score.list, function(x) x$AICc[2])

    output <- rbind(data.frame(model=score.list[[1]]$H0,
        X.squared = NA,
        df = NA,
        p.value = NA,
        npar = score.list[[1]]$np0,
        AIC = score.list[[1]]$AIC[1],
        AICc = score.list[[1]]$AICc[1]), output)

    if (nrow(output)>1) {
        output$dAICc <- output$AICc - min(output$AICc)
        OK <- abs(output$dAICc) < abs(dmax)
        sumdAICc <- sum(exp(-output$dAICc[OK]/2))
        output$AICwt <- ifelse ( OK, exp(-output$dAICc/2) / sumdAICc, 0)
    }

    output[,c('AIC','AICc','dAICc')] <- round(output[,c('AIC','AICc','dAICc')],3)
    if (!is.null(output$AICwt)) output[,'AICwt'] <- round(output[,'AICwt'],3)

    if (sort) output <- output [order(output$AICc),]
    else      output <- output [,]
    row.names(output) <- NULL

    output
}
############################################################################################

# print(score.test (secr0CL, list(g0=~b, sigma=~1)))
# m1 <- list(D = ~1, g0 = ~b, sigma = ~1)
# m2 <- list(D = ~1, g0 = ~1, sigma = ~b)
############################################################################################

