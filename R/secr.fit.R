################################################################################
## package 'secr'
## secr.fit.R
## moved from methods.R 2011-01-30
## 2011-10-20 generalized designD
## 2011-10-20 renamed 'grps' 'grouplevels'
## 2011-12-16 streamlined preparation of models
################################################################################

secr.fit <- function (capthist, model = list(D~1, g0~1, sigma~1), mask = NULL,
    buffer = NULL, CL = FALSE, detectfn = NULL, binomN = NULL, start = NULL,
    link = list(), fixed = list(), timecov = NULL, sessioncov = NULL,
    groups = NULL, dframe = NULL, details = list(), method = 'Newton-Raphson',
    verify = TRUE, trace = NULL, ...)

{
# Fit spatially explicit capture recapture model
#
# Arguments:
#
#  capthist   -  capture history object (includes traps object as an attribute)
#  model      -  formulae for real parameters in terms of effects and covariates
#  mask       -  habitat mask object
#  buffer     -  default buffer width should mask not be provided
#  CL         -  logical switch : conditional likelihood (T) or full likelihood (F)
#  detectfn   -  code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential etc.
#  start      -  start values for maximization (numeric vector link scale);
#                if NULL then 'autoini' function is used
#  link       -  list of parameter-specific link function names 'log', 'logit', 'identity',
#                'sin', 'neglog'
#  fixed      -  list of fixed values for named parameters
#  timecov    -  data for time covariates if these are used in 'model'
#  sessioncov -  dataframe of session-level covariates
#  groups     -  vector of names to group fields in attr(capthist,'covariates') dataframe
#  dframe     -  optional data frame of design data for detection model (tricky & untested)
#  details    -  list with several additional settings, mostly of special interest
#  method     -  optimization method (indirectly chooses
#  verify     -  logical switch for pre-check of capthist and mask with verify()
#  trace      -  logical; if TRUE output each likelihood as it is calculated
#  ...        -  other arguments passed to nlm() or optim()


    if (!inherits(capthist, 'capthist'))
        stop ("requires 'capthist' object")
    if ((detector(traps(capthist))=='cue') & (is.null(groups) | !CL))
        stop ("cue detector requires CL = TRUE and groups")

    #################################################
    ## Remember start time and call

    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")

    cl   <- match.call(expand.dots = TRUE)

    #################################################
    ## Default detection function

    if (is.null(detectfn)) {
        if (detector(traps(capthist)) %in% c('cue', 'signal')) {
            detectfn <- 10
            warning ("detectfn not specified; using signal strength (10)")
        }
        else {
            detectfn <- 0
        }
    }

    else {
        if (detector(traps(capthist)) == 'presence')
            detectfn <- valid.detectfn(detectfn, 0:8)
        else
            detectfn <- valid.detectfn(detectfn)
    }

    #################################################
    ## Use input 'details' to override various defaults

    defaultdetails <- list(distribution = 'poisson', scalesigma = FALSE,
        scaleg0 = FALSE, hessian = 'auto', trace = TRUE, LLonly = FALSE,
        centred = FALSE, binomN = 1, cutval = 0,
        minprob = 1e-50, tx = 'identity', param = 0)
    if (detector(traps(capthist)) %in% .localstuff$countdetectors)
        defaultdetails$binomN <- 0   ## Poisson
    if (!is.null(attr(capthist,'cutval')))
        defaultdetails$cutval <- attr(capthist,'cutval')
    if (is.logical(details$hessian))
        details$hessian <- ifelse(details$hessian, 'auto', 'none')
    details <- replace (defaultdetails, names(details), details)
    if (!is.null(trace)) details$trace <- trace
    if (!is.null(binomN)) details$binomN <- binomN   ## 2011 01 28
    if (details$LLonly)  details$trace <- FALSE
    if (detector(traps(capthist)) != 'multi') details$param <- 0

    ## 2011-02-06 to be quite clear -
    if (detector(traps(capthist)) %in% c(.localstuff$exclusivedetectors, 'proximity','signal'))
        details$binomN <- 1;


    #################################################
    ## NEW multi-session capthist object 12/2/09
    ## MS - indicator TRUE if multi-session (logical)
    ## sessionlevels - names of sessions (character)

    MS <- ms(capthist)
    sessionlevels <- session(capthist)

    if (is.null(sessionlevels)) sessionlevels <- '1'
    anycount <- any(detector(traps(capthist)) %in% .localstuff$countdetectors) &
        (details$binomN != 1)    ## 2011-01-29
    anypoly  <- any(detector(traps(capthist)) %in% c('polygon',  'polygonX'))
    anytrans <- any(detector(traps(capthist)) %in% c('transect', 'transectX'))

    if (MS) {
       if (any (sapply(traps(capthist), detector) == 'single'))
        warning ("multi-catch likelihood used for single-catch traps")
    }
    else {
       if (detector(traps(capthist)) == 'single')
        warning ("multi-catch likelihood used for single-catch traps")
    }

    #################################################
    ## Optional data check added 2009 09 19

    if (verify) {
        memo ('Checking data', details$trace)
        test <- verify(capthist, report = 1)
        if (test$errors)
            stop ("'verify' found errors in 'capthist' argument")

        if (!is.null(mask)) {
            if (MS & ms(mask)) {
                ## list of masks
                test <- lapply(mask, verify, report = 1)
                notOK <- any(unlist(test))
            }
            else notOK <- verify(mask, report = 1)$errors
            if (notOK)
                stop ("'verify' found errors in 'mask' argument")
        }
    }

    #################################################
    ## Ensure valid mask
    ## assume traps(capthist) will extract a list of trap layouts
    ## if multi-session (MS == TRUE)

    usebuffer <- is.null(mask)    ## flag for later check
    if (usebuffer) {
        if (is.null(buffer)) {
            buffer <- 100
            if (!(detector(traps(capthist))=='presence'))
                warning ("using default buffer width 100 m")
        }
        if (MS) mask <- lapply (traps(capthist), make.mask, buffer = buffer)
        else    mask <- make.mask(traps(capthist), buffer = buffer)
    }
    else {
      if (MS & !ms(mask)) {
          ## inefficiently replicate mask for each session!
          mask <- lapply(sessionlevels, function(x) mask)
          class (mask) <- c('list', 'mask')
          names(mask) <- sessionlevels
      }
    }

    nc <- ifelse (MS, sum(sapply(capthist, nrow)), nrow(capthist))
    if (nc < 1)
        warning (nc, " detection histories")

    if (MS) {
        q  <- attr(capthist[[1]],'q')
        Tm <- attr(capthist[[1]],'Tm')
    }
    else {
        q <- attr(capthist,'q')
        Tm <- attr(capthist,'Tm')
    }
    nonID <- !is.null(Tm)   ## were marked animals recorded if unidentified?
    if (!is.null(q) & CL)
        stop ("mark-resight incompatible with CL")

    #################################################
    ## optional centring of traps and mask 2010 04 27
    if (details$centred) {
        centre <- function (xy, dxy) {
            xy[,] <- sweep(xy, MARGIN = 2, FUN='-', STATS = dxy)
            xy
        }
        if (MS) {
            nsess <- length(traps(capthist))
            offsetxy <- lapply(traps(capthist), function(xy) apply(xy, 2, mean))
            for (i in 1:nsess) {
                temptraps <- centre(traps(capthist[[i]]), offsetxy[[i]])
                traps(capthist[[i]]) <- temptraps
                mask[[i]] <- centre(mask[[i]], offsetxy[[i]])
                attr(mask[[i]], 'meanSD')[1,1:2] <- attr(mask[[i]], 'meanSD')[1,1:2] -
                    offsetxy[[i]]
                attr(mask[[i]], 'boundingbox') <- centre(attr(mask[[i]], 'boundingbox'),
                    offsetxy[[i]])
            }
        }
        else {
            offsetxy <- apply(traps(capthist), 2, mean)
            traps(capthist) <- shift(traps(capthist), -offsetxy)
            mask <- shift.traps(mask, -offsetxy)
            attr(mask, 'meanSD')[1,1:2] <- attr(mask, 'meanSD')[1,1:2] - offsetxy
            attr(mask, 'boundingbox') <- centre(attr(mask, 'boundingbox'), offsetxy)
        }
    }

    #################################################
    ## Use input formula to override defaults

    if ('formula' %in% class(model)) model <- list(model)
    model <- stdform (model)  ## named, no LHS
    defaultmodel <- list(D=~1, g0=~1, sigma=~1, z=~1, w=~1, pID=~1,
        beta0=~1, beta1=~1, sdS=~1, b0=~1, b1=~1, phi=~1)
    model <- replace (defaultmodel, names(model), model)

#    if (!(detectfn %in% c(1,3,7,8))) model$z <- NULL
#    if (!(detectfn %in% c(5,6))) model$w <- NULL
#    if (!(detectfn %in% c(9))) { model$b0 <- NULL; model$b1 <- NULL }
#    if (!(detectfn %in% c(10,11))) { model$beta0 <- NULL; model$beta1 <- NULL; model$sdS <- NULL }
#    if (!(detectfn %in% 0:8)) { model$g0 <- NULL; model$sigma <- NULL }
#    if (is.null(details$intervals)) { model$phi <- NULL }
#    if (is.null(q) | !nonID)
#        model$pID <- NULL    ## no use for this parameter
#    else
#        if (model$pID != ~1)
#            stop ("'pID' must be constant in this implementation")
#     if (CL) model$D <- NULL

    pnames <- switch (detectfn+1,
        c('g0','sigma'),           # 0 halfnormal
        c('g0','sigma','z'),       # 1 hazard rate
        c('g0','sigma'),           # 2 exponential
        c('g0','sigma','z'),       # 3
        c('g0','sigma'),           # 4
        c('g0','sigma','w'),       # 5
        c('g0','sigma','w'),       # 6
        c('g0','sigma','z'),       # 7
        c('g0','sigma','z'),       # 8
        c('b0','b1'),              # 9
        c('beta0','beta1','sdS'),  # 10
        c('beta0','beta1','sdS'))  # 11

    if (!CL) pnames <- c('D', pnames)
    if (!is.null(details$intervals)) pnames <- c(pnames, 'phi')
    if (!is.null(q) & nonID) {
	pnames <- c(pnames, 'pID')
        if (model$pID != ~1)
            stop ("'pID' must be constant in this implementation")
    }
    pnames <- pnames[!(pnames %in% names(fixed))]  ## drop fixed real parameters
    model[!(names(model) %in% pnames)] <- NULL     ## select real parameters
    vars <-  unlist(lapply(model, all.vars))

    ############################################
    # Finite mixtures - 2009 12 10, 2011 12 16
    ############################################
    nmix <- get.nmix(model)
    if ((nmix>1) & (nmix<4)) {
        model$pmix <- as.formula(paste('~h', nmix, sep=''))
        if (!all(all.vars(model$pmix) %in% c('session','g','h2','h3')))
            stop ("formula for pmix may include only 'session', 'g' or '1'")
        pnames <- c(pnames, 'pmix')
    }
    details$nmix <- nmix

    #################################################
    ## Specialisations

    if (CL & !(is.null(groups) | (detector(traps(capthist))=='cue'))) {
        groups <- NULL
        warning ("groups not valid with CL; groups ignored")
    }
    if (CL && var.in.model('g', model))
        stop ("'g' is not a valid effect when 'CL = TRUE'")
    if ((length(model) == 0) & !is.null(fixed))
        stop ("all parameters fixed")     ## assume want only LL
    if (details$scalesigma) {
        if (CL)
            stop ("cannot use 'scalesigma' with 'CL'")
        if (!is.null(fixed$D))
            stop ("cannot use 'scalesigma' with fixed density")
        if (!(model$D == formula(~1) |
              model$D == formula(~session)))
            stop ("cannot use 'scalesigma' with inhomogenous density or groups")
        if (!is.null(groups))
            stop ("cannot use 'scalesigma' with groups")
    }
    if (details$scaleg0) {
        if (!is.null(groups))
            stop ('Cannot use scaleg0 with groups')
    }

    #################################
    # Link functions (model-specific)
    #################################
    defaultlink <- list(D='log', g0='logit', sigma='log', z='log', w='log', pID='logit',
        beta0='identity', beta1='neglog', sdS='log', b0='log', b1='neglog', pmix='logit',
        cuerate='log', phi = 'logit')
    if (anycount) defaultlink$g0 <- 'log'
    link <- replace (defaultlink, names(link), link)
    link[!(names(link) %in% pnames)] <- NULL

#    if (CL) link$D <- NULL
#    if (!(detectfn %in% c(1,3,7,8))) link$z <- NULL
#    if (!(detectfn %in% c(5,6))) link$w <- NULL
#    if (!(detectfn %in% c(9))) { link$b0 <- NULL; link$b1 <- NULL }
#    if (!(detectfn %in% c(10,11))) { link$beta0 <- NULL; link$beta1 <- NULL; link$sdS <- NULL }
#    if (!(detectfn %in% 0:8)) { link$g0 <- NULL; link$sigma <- NULL }
#    if (is.null(q) | !nonID) link$pID <- NULL
#    if (details$nmix==1) link$pmix <- NULL

    if (details$scaleg0) link$g0 <- 'log'  ## Force log link in this case as no longer 0-1
    if (!(detector(traps(capthist))=='cue')) link$cuerate <- NULL

    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################

    memo ('Preparing detection design matrices', details$trace)
    design <- secr.design.MS (capthist, model, timecov, sessioncov, groups, dframe)
    design0 <- secr.design.MS (capthist, model, timecov, sessioncov, groups, dframe,
        naive = T, bygroup = !CL)

    ##############################################
    # Prepare turnover design matrices and lookup
    # experimental addition 2011-11-30
    ##############################################

    if (!is.null(details$intervals)) {
        memo ('Preparing turnover design matrices', details$trace)
        intervalcov <- NULL
        nTurnoverParameters <- integer(0)
    }
    designphi <- phi.designdata (capthist, model, intervalcov, sessioncov,
        groups, dframe, intervals = details$intervals)

    ###############################
    # Prepare density design matrix
    ###############################
    D.modelled <- !CL & is.null(fixed$D)
    if (!D.modelled) {
       designD <- matrix(nrow = 0, ncol = 0)
       grouplevels <- 1    ## was NULL
       attr(designD, 'dimD') <- NA
       nDensityParameters <- integer(0)
    }
    else {
        grouplevels  <- group.levels(capthist,groups)
        if (!is.null(details$userDfn)) {
            ## may provide a function used by getD in functions.R
            ## userDfn(mask, beta[parindx$D], ngrp, nsession)
            designD <- details$userDfn
            if (!is.function(designD))
                stop ("details$userDfn should be a function")
            ## this form of call returns only coefficient names
            Dnames <- designD('parameters', mask)
        }
        else {
            memo ('Preparing density design matrix', details$trace)
            temp <- D.designdata( mask, model$D, grouplevels, sessionlevels, sessioncov)
            D.designmatrix <- model.matrix(model$D, temp)
            attr(D.designmatrix, 'dimD') <- attr(temp, 'dimD')
            Dnames <- colnames(D.designmatrix)
            designD <- D.designmatrix
        }
        nDensityParameters <- length(Dnames)
    }

    #############################
    # Parameter mapping (general)
    #############################

    np <- sapply(design$designMatrices, ncol)
    np <- c(D = nDensityParameters, np)
    if (!is.null(details$intervals))
        np <- c(np,  sapply(designphi$designMatrices, ncol))
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)
    if (!D.modelled) parindx$D <- NULL

    ###################################################
    # Option to generate start values from previous fit
    ###################################################

    if (inherits(start, 'secr')) {
        ## use 'mapbeta' from score.test.R
        start <- mapbeta(start$parindx, parindx, coef(start)$beta, NULL)
    }

    #############################
    # Single evaluation option
    #############################
    .localstuff$iter <- 0
    if (details$LLonly) {
      if (is.null(start))
          stop ("must provide transformed parameter values in 'start'")
      if (!is.null(q))
          stop ("not for mark-resight")

      LL <- - secr.loglikfn (beta = start,
                       parindx    = parindx,
                       link       = link,
                       fixed      = fixed,
                       designD    = designD,
                       design     = design,
                       design0    = design0,
                       designphi  = designphi,
                       capthist   = capthist,
                       mask       = mask,
                       detectfn   = detectfn,
                       CL         = CL,
                       groups     = groups,
                       details    = details,
                       logmult    = TRUE,     ## add if possible
                       )

      return(c(logLik=LL))
    }

    ###############################
    # Start values (model-specific)
    ###############################
    ## 'start' is vector of beta values (i.e. transformed)

    if (!is.null(start)) {
        if (detector(traps(capthist))=='cue')
            stopifnot (length(start) == (NP+1))
        else
            stopifnot (length(start) == NP)
    }
    else {

        start3 <- list(D=NA, g0=NA, sigma=NA)

          if (!(detectfn %in% c(9,10,11)) && !anypoly && !anytrans) {  ## not for signal attenuation
            memo('Finding initial parameter values...', details$trace)
            ## autoini uses default buffer dbar * 4
            if (MS)
                start3 <- autoini (capthist[[1]], mask[[1]]) ## Use session 1 - can be risky
            else
                start3 <- autoini (capthist, mask)
            if (any(is.na(unlist(start3)))) {
                warning ("'secr.fit' failed because initial values not found",
                    " (data sparse?); specify transformed values in 'start'")
                return (list(call=cl, fit=NULL))
            }

            if (details$scaleg0 & anycount)
                stop ("'scaleg0' not compatible with count detectors")
            ## next two stmts must be this order (g0 then sigma)
            if (details$scaleg0) start3$g0 <- start3$g0 * start3$sigma^2
            if (details$scalesigma) start3$sigma <- start3$sigma * start3$D^0.5

            memo(paste('Initial values ', paste(paste(c('D', 'g0', 'sigma'),
                '=', round(unlist(start3),5)),collapse=', ')),
                details$trace)
        }
        else warning ("using default starting values")

        ## assemble start vector
        default <- list(
            D     = ifelse (is.na(start3$D), 1, start3$D),
            g0    = ifelse (is.na(start3$g0), 0.1, start3$g0),
            sigma = ifelse (is.na(start3$sigma), unlist(RPSV(capthist))[1], start3$sigma),
            z     = 5,
            w     = 10,
            pID   = 0.7,
            beta0 = details$cutval + 30,
            beta1 = -0.2,
            sdS   = 2,
            b0    = 2,        ## changed from 15 2010-11-01
            b1    = -0.1,
            pmix  = 0.25,
            phi   = 0.7
        )
        if (detectfn %in% c(6)) {
            default$w <- default$sigma
            default$sigma <- default$sigma/2
        }
        if (detectfn %in% c(7)) {
            default$z <- default$sigma/5
        }
        if (detectfn %in% c(8)) {
            default$z <- 1    ## cumulative gamma
        }
        if (anypoly | anytrans) {
            if (MS) {
                tempcapthist <- capthist[[1]]
                tempmask <- mask[[1]]
            }
            else {
                tempcapthist <- capthist
                tempmask <- mask
            }
            default$D <- 2 * nrow(tempcapthist) / (nrow(tempmask)*attr(tempmask,'area'))
            default$g0 <- sum(tempcapthist) / nrow(tempcapthist) / ncol(tempcapthist)
            default$sigma <- RPSV(tempcapthist)
        }
        if (is.na(default$sigma)) default$sigma <- 20
        getdefault <- function (par) transform (default[[par]], link[[par]])

        start <- rep(0, NP)
        for ( i in 1:length(parindx) )
            start[parindx[[i]][1]] <- getdefault (names(model)[i])
        if(details$nmix>1)
             ## scaled by mlogit.untransform
            start[parindx[['pmix']]] <- (2:details$nmix)/(details$nmix+1)
        if (detector(traps(capthist))=='cue')
            start <- c(start, log(3))    ## cuerate
        ## could be nrow(capthist)/length(levels(covariates(capthist)$animal

        # D/ngrp when figure out where to calculate this

        ## if (!(is.null(q) | !nonID) & is.null(fixed$pID))
        ##     start[parindx$pID[1]] <- getdefault('pID')
    }

    ##########################
    # Fixed beta parameters
    ##########################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        start <- start[is.na(fb)]  ## drop unwanted betas; remember later to adjust parameter count
    }
    ##########################
    # Variable names (general)
    ##########################
    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    realnames <- names(model)
    if (D.modelled) betanames <- c(paste('D', Dnames, sep='.'), betanames)
    if (!is.null(details$intervals))   ## allow turnover 2011-11-30
        betanames <- c(betanames, unlist(sapply(designphi$designMatrices, colnames)))
    betanames <- sub('..(Intercept))','',betanames)
    ## allow for fixed beta parameters 2009 10 19
    if (!is.null(details$fixedbeta))
        betanames <- betanames[is.na(details$fixedbeta)]

    #################################
    # Variable names (model-specific)
    #################################

    if (detector(traps(capthist))=='cue') {
        betanames <- c(betanames, 'cuerate')
        realnames <- c(realnames, 'cuerate')
    }
    betaw <- max(max(nchar(betanames)),8)   # for 'trace' formatting

    #####################
    # Maximize likelihood
    #####################

    memo('Maximizing likelihood...', details$trace)
    if (details$trace) cat('Eval     Loglik', formatC(betanames, format='f', width=betaw), '\n')

    loglikefn <- secr.loglikfn
    if (!is.null(q))
        stop ("mark-resight option not operative")
        ## loglikefn <- MRsecr.loglikfn

    if (tolower(method) %in% c('newton-raphson', 'nr')) {
        args <- list (p         = start,
                        f          = loglikefn,
                        link       = link,
                        fixed      = fixed,
                        parindx    = parindx,
                        capthist   = capthist,
                        mask       = mask,
                        CL         = CL,
                        detectfn   = detectfn,
                        designD    = designD,
                        design     = design,
                        design0    = design0,
                        designphi  = designphi,
                        groups     = groups,
                        details    = details,
                        logmult    = TRUE,     ## add if possible
                        betaw      = betaw,   # for trace format
                        hessian    = tolower(details$hessian)=='auto',
                        stepmax    = 10)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (nlm, args)
        this.fit$par <- this.fit$estimate     # copy for uniformity
        this.fit$value <- this.fit$minimum    # copy for uniformity
        if (this.fit$code > 2)
            warning ("possible maximization error: nlm returned code ",
                this.fit$code, ". See ?nlm")
    }
    else {
        args <- list(par     = start,
                        fn         = loglikefn,
                        link       = link,
                        fixed      = fixed,
                        parindx    = parindx,
                        capthist   = capthist,
                        mask       = mask,
                        CL         = CL,
                        detectfn   = detectfn,
                        designD    = designD,
                        design     = design,
                        design0    = design0,
                        designphi  = designphi,
                        groups     = groups,
                        details    = details,
                        logmult    = TRUE,     ## add if possible
                        betaw      = betaw,   # for trace format
                        hessian    = tolower(details$hessian)=='auto',
                        method     = method)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (optim, args)
        # default method = 'BFGS', control=list(parscale=c(1,0.1,5))
        if (this.fit$convergence != 0)
            warning ("probable maximization error: optim returned convergence ",
                this.fit$convergence, ". See ?optim")
    }

    this.fit$method <- method         ## remember what method we used...
    covar <- NULL
    N <- NULL
    if (this.fit$value > 1e9) {     ## failed
        this.fit$beta[] <- NA
    }
    else {


        ############################
        # Variance-covariance matrix
        ############################

        if (tolower(details$hessian)=='fdhess') {
            require (nlme)
            memo ('Computing Hessian with fdHess in nlme', details$trace)
            loglikfn <- function (beta) {
               -secr.loglikfn(
                            beta       = beta,
                            parindx    = parindx,
                            link       = link,
                            fixed      = fixed,
                            designD    = designD,
                            design     = design,
                            design0    = design0,
                            designphi  = designphi,
                            capthist   = capthist,
                            mask       = mask,
                            detectfn   = detectfn,
                            CL         = CL,
                            groups     = groups,
                            details    = details,
                            logmult    = TRUE,     ## add if possible
                            betaw      = betaw,    ## for trace format
                     )
            }
            grad.Hess <- fdHess(this.fit$par, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
            this.fit$hessian <- -grad.Hess$Hessian
        }

        if (!is.null(this.fit$hessian)) {
            covar <- try(solve(this.fit$hessian))
            if (inherits(covar, "try-error")) {
                warning ("could not invert Hessian to compute ",
                         "variance-covariance matrix")
                covar <- matrix(nrow = NP, ncol = NP)
            }
            else if (any(diag(covar)<0)) {
                if (method == "BFGS")
                    suffix <- ""
                else
                    suffix <- " - try method = 'BFGS'"
                warning ("variance calculation failed ", suffix)
                covar <- matrix(nrow = NP, ncol = NP)
            }
            dimnames(covar) <- list(betanames, betanames)
        }

        ## predicted D across mask
        if (!CL) {
            D <- getD (designD, this.fit$par, mask, parindx, link, fixed,
                       grouplevels, sessionlevels)
            N <- t(apply(D, 2:3, sum, drop = FALSE))
            cellarea <- if (ms(mask)) sapply(mask, attr, 'area')
                        else cellarea <- attr(mask,'area')
            N <- sweep(N, FUN = '*', MARGIN = 1, STATS = cellarea)
        }
    }

    desc <- packageDescription("secr")  ## for version number

    temp <- list (call = cl,
                  capthist = capthist,
                  mask = mask,
                  detectfn = detectfn,
                  CL = CL,
                  timecov = timecov,
                  sessioncov = sessioncov,
                  groups = groups,
                  dframe = dframe,
                  design = design,      ## added 2009 09 05
                  design0 = design0,    ## added 2009 06 25
                  designphi = designphi, ## added 2011 11 30

                  start = start,        ## added 2009 09 09
                  link = link,
                  fixed = fixed,
                  parindx = parindx,
                  model = model,
                  details = details,

                  vars = vars,
                  betanames = betanames,
                  realnames = realnames,

                  fit = this.fit,
                  beta.vcv = covar,
#                 D = D,                   ## dropped 2011-11-10
                  N = N,                   ## added 2011-11-10
                  version = desc$Version,  ## added 2009 09 21
                  starttime = starttime,   ## added 2009 09 21
                  proctime = (proc.time() - ptm)[1]
             )

    attr (temp, 'class') <- 'secr'

    ## check introduced 2010-12-01 & adjusted 2011-09-28
    ## bias.D not for polygon & transect detectors
    if (verify & usebuffer & (this.fit$value < 1e9) &
        (detector(traps(capthist)) %in% .localstuff$pointdetectors) &
        !(detector(traps(capthist)) %in% c('cue','unmarked','presence'))) {
        if (MS) {
            nsess <- length(capthist)
            bias <- numeric(nsess)
            for (i in 1:nsess) {
                bias[i] <- bias.D(buffer, traps(capthist)[[i]], detectfn = temp$detectfn,
                    detectpar = detectpar(temp)[[i]], noccasions = ncol(capthist[[i]]))$RB.D
            }
            if (any(bias > 0.01))
                warning ("predicted relative bias exceeds 0.01 with ",
                         " buffer = ", buffer)
        }
        else {
            bias <- bias.D(buffer, traps(capthist), detectfn=temp$detectfn,
                detectpar = detectpar(temp), noccasions = ncol(capthist))$RB.D
            if (bias > 0.01)
                warning ("predicted relative bias exceeds 0.01 with ",
                         "buffer = ", buffer)
        }
    }

    memo(paste('Completed in ', round(temp$proctime,2), ' seconds at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"),
        sep=''), details$trace)
    temp

}
