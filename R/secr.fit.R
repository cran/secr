################################################################################
## package 'secr'
## secr.fit.R
## moved from methods.R 2011-01-30
## 2011-10-20 generalized designD
## 2011-10-20 renamed 'grps' 'grouplevels'
## 2011-12-16 streamlined preparation of models
## 2012-01-22 purged phi/turnover
## 2012-01-31 experimental addition of parameter cut
## 2012-04-06 'fixed' bug fixed (see functions.r)
## 2012-07-24 unmash component of details
## 2013-04-11 hcov
## 2013-04-19 lambda0
## 2013-04-20 default mask type changed to trapbuffer
## 2013-07-02 esa parameterisation details$param = 2
## 2013-07-19 a0 parameterisation details$param = 3
## 2013-10-14 tweak a0 biasD
## 2013-10-28 revise pmix models
## 2013-11-09 rearrange model check code to put in warning
## 2014-03-12 sigmak parameterisation details$param = 4
## 2014-03-12 bufferbiascheck now in suggest.buffer
## 2014-03-24 allow combination esa, sigmak (param 6)
###############################################################################

secr.fit <- function (capthist, model = list(D~1, g0~1, sigma~1), mask = NULL,
    buffer = NULL, CL = FALSE, detectfn = NULL, binomN = NULL, start = NULL,
    link = list(), fixed = list(), timecov = NULL, sessioncov = NULL, hcov = NULL,
    groups = NULL, dframe = NULL, details = list(), method = 'Newton-Raphson',
    verify = TRUE, biasLimit = 0.01, trace = NULL, ncores = 1, ...)

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


    #################################################
    ## Remember start time and call

    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")

    cl   <- match.call(expand.dots = TRUE)

    ## 2014-02-13
    if (is.character(capthist)) {
        capthist <- get(capthist, pos=-1)
    }
    if (is.character(mask)) {
        mask <- get(mask, pos=-1)
    }
    if (is.character(dframe)) {
        dframe <- get(dframe, pos=-1)
    }

    if (!inherits(capthist, 'capthist'))
        stop ("requires 'capthist' object")
    if ((detector(traps(capthist))=='cue') & (is.null(groups) | !CL))
        stop ("cue detector requires CL = TRUE and groups")

    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = FALSE)
        clusterEvalQ(clust, library(secr))
    }
    else
        clust <- NULL

    #################################################
    ## Default detection function

    if (is.null(detectfn)) {
        if (detector(traps(capthist)) %in% c('cue', 'signal')) {
            detectfn <- 10
            warning ("detectfn not specified; using signal strength (10)")
        }
        else if (detector(traps(capthist)) %in% c('signalnoise')) {
            detectfn <- 12
            warning ("detectfn not specified; using signal-noise (12)")
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

    defaultdetails <- list(distribution = 'poisson',
                           scalesigma = FALSE,
                           scaleg0 = FALSE,
                           hessian = 'auto',
                           trace = TRUE,
                           LLonly = FALSE,
                           centred = FALSE,
                           binomN = 1,
                           cutval = 0,
                           minprob = 1e-50,
                           tx = 'identity',
                           param = 0,
                           unmash = FALSE,
                           telemetrytype = 'concurrent',
                           telemetrysigma = FALSE,
                           telemetrybvn = FALSE,
                           ignoreusage = FALSE,
                           debug = FALSE,
                           intwidth2 = 0.8,
                           normalize = FALSE,
                           usecov = NULL
##                           , maxcallsize = 512
                           )

    if (detector(traps(capthist)) %in% .localstuff$countdetectors)
        defaultdetails$binomN <- 0   ## Poisson
    if (!is.null(attr(capthist,'cutval')))
        defaultdetails$cutval <- attr(capthist,'cutval')
    else if (ms(capthist) & !is.null(attr(capthist[[1]],'cutval')))   ## 2012-09-04
        defaultdetails$cutval <- attr(capthist[[1]],'cutval')
    if (is.logical(details$hessian))
        details$hessian <- ifelse(details$hessian, 'auto', 'none')
    details$telemetrytype <- match.arg(details$telemetrytype, c('none', 'independent',
        'dependent', 'concurrent'))
    if (details$telemetrytype == 'independent')
        details$telemetrysigma <- TRUE
    details <- replace (defaultdetails, names(details), details)
    if (details$telemetrytype == 'none' && details$telemetrysigma == TRUE)
        details$telemetrytype <- 'independent'
    if (!is.null(trace)) details$trace <- trace
    if (!is.null(binomN)) {
        if (detector(traps(capthist)) == 'count') {
            if (tolower(binomN) == 'usage')
                binomN <- 1   ## code for 'binomial size from usage' 2012-12-22
            if (tolower(binomN) == 'poissonhazard')
                binomN <- -1  ## interpret g() as true detection function 2013-01-07
        }
        details$binomN <- binomN   ## 2011 01 28
    }
    if (details$LLonly)  details$trace <- FALSE

    if (!(detector(traps(capthist)) %in% c('single','multi')) & (details$param == 1)) {
        warning ("Gardner & Royle parameterisation not appropriate, using param = 0")
        details$param <- 0
    }

    ## 2011-02-06 to be quite clear -
    if (detector(traps(capthist)) %in% c(.localstuff$exclusivedetectors,
                                         'proximity','signal','signalnoise'))
        details$binomN <- 1;


    #################################################
    ## MS - indicator TRUE if multi-session (logical)
    ## sessionlevels - names of sessions (character)

    MS <- ms(capthist)
    sessionlevels <- session(capthist)

    if (is.null(sessionlevels)) sessionlevels <- '1'
    anycount <- any(detector(traps(capthist)) %in% .localstuff$countdetectors)
    anypoly  <- any(detector(traps(capthist)) %in% c('polygon',  'polygonX'))
    anytrans <- any(detector(traps(capthist)) %in% c('transect', 'transectX'))
    alltelem <- all(detector(traps(capthist)) %in% c('telemetry'))
    if (alltelem) CL <- TRUE

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
            if (!(detector(traps(capthist))=='presence') & !alltelem)
                warning ("using default buffer width 100 m")
        }
        if (MS) mask <- lapply (traps(capthist), make.mask, buffer = buffer, type = "trapbuffer")
        else    mask <- make.mask(traps(capthist), buffer = buffer, type = "trapbuffer")
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

    #################################################
    ## orphan mark-resight code - not currently used

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
    if (sum(c('g0','a0','lambda0') %in% names(model)) > 1)
        stop ("model should include only one of g0, a0, lambda0")

    esa <- "esa" %in% names(model)
    a0 <- "a0" %in% names(model)
    sigmak <- "sigmak" %in% names(model)
    if (esa & !sigmak) {
        details$param <- 2
    }
    if (a0 & !sigmak) {
        details$param <- 3
    }
    if (sigmak) {
        if (esa) {
            details$param <- 6
        }
        else {
            if (CL)
                stop ("sigmak parameterization requires full model, not CL, unless also 'esa'")
            details$param <- ifelse(a0, 5, 4)
        }
    }
    if (details$param > 0)
        warning ("Using parameterization details$param = ", details$param)

    ## modelled parameters
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
        c('beta0','beta1','sdS'),  # 11
        c('beta0','beta1','sdS'),  # 12  cf parnames() in utility.R: muN, sdN?
        c('beta0','beta1','sdS'),  # 13  cf parnames() in utility.R: muN, sdN?
        c('lambda0','sigma'),      # 14 hazard halfnormal
        c('lambda0','sigma','z'),  # 15 hazard hazard rate
        c('lambda0','sigma'),      # 16 hazard exponential
        c('lambda0','sigma','w'),  # 17
        c('lambda0','sigma','z'))  # 18

    if (details$param %in% c(2,6))
        pnames[1] <- "esa"
    if (details$param %in% c(3,5))
        pnames[1] <- "a0"
    if (details$param %in% 4:6) {
        pnames[2] <- "sigmak"
        pnames <- c(pnames, "c")
    }

    if (!CL) pnames <- c('D', pnames)
    if (!is.null(q) & nonID) {
	pnames <- c(pnames, 'pID')
        if (model$pID != ~1)
            stop ("'pID' must be constant in this implementation")
    }
    if (alltelem) {
        rnum <- match(c('D','lambda0','a0','esa','g0'), pnames)
        rnum[is.na(rnum)] <- 0
        pnames <- pnames[-rnum]
    }
    if (details$param %in% 4:6) {
        if (! ("c" %in% names(model))) {
            ## default to fixed c = 0
            if (!("c" %in% names(fixed)))
                fixed$c <- 0
        }
    }
    fnames <- names(fixed)
    pnames <- pnames[!(pnames %in% fnames)]        ## drop fixed real parameters

    ## this warning may be superfluous, and distracting
    ## OK <- names(model) %in% pnames
    ## if (any(!OK))
    ##    warning ("parameters in model not consistent with detectfn : ",
    ##             paste(names(model)[!OK], collapse = ', '))

    ## 2013-11-09 only now fillout model
    ## pmix, lambda0 added to defaultmodel 2013-04
    defaultmodel <- list(D=~1, g0=~1, lambda0=~1,  esa=~1, a0=~1, sigma=~1, sigmak=~1, z=~1,
                         w=~1, c=~1, pID=~1, beta0=~1, beta1=~1, sdS=~1, b0=~1, b1=~1)
                         ## pmix=~1) ## 2013-10-27
    model <- replace (defaultmodel, names(model), model)

    ######################################################
    # Finite mixtures - 2009 12 10, 2011 12 16, 2013 10 27
    ######################################################

    nmix <- get.nmix(model, capthist, hcov)
    if (nmix > 3)
        stop ("number of latent classes exceeds 3")
    if ((nmix>1) & !is.null(hcov) & !is.null(groups))
        stop ("hcov mixture model incompatible with groups")
    if ((nmix == 1) & ('pmix' %in% c(fnames,names(model))))
        stop ("pmix specified for invariant detection model")

    if ((nmix>1) & !('pmix' %in% fnames)) {
        if (is.null(model$pmix))
            model$pmix <- ~1
        pmixvars <- all.vars(model$pmix)
        if (!any (pmixvars %in% c('h2','h3'))) {
            model$pmix <- if (nmix == 2)
                update(model$pmix, ~. + h2)
            else
                update(model$pmix, ~. + h3)
        }
        else {
            ## badvar <- pmixvars %in% c('t','T','b','B','bk')
            badvar <- !(pmixvars %in% c('session','Session',sessioncov,'h2','h3'))
            if (any(badvar))
                stop ("formula for pmix may not include ", pmixvars[badvar])
        }
        pnames <- c(pnames, 'pmix')
    }
    else
        model$pmix <- NULL
    details$nmix <- nmix

    ############################################
    ## Tidy up
    ############################################

    model[!(names(model) %in% pnames)] <- NULL     ## select real parameters
    vars <-  unlist(lapply(model, all.vars))

    ############################################
    ## Specialisations
    ############################################
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
        warning ("scaleg0 is deprecated and will soon be dropped from secr; use a0 parameterization instead (param == 3)")
        if (!is.null(groups))
            stop ('Cannot use scaleg0 with groups')
    }
    if (details$scaleg0) {
        warning ("scalesigma is deprecated and will soon be dropped from secr; use sigmak parameterization instead (param == 4")
    }
    ############################################
    # Link functions (model-specific)
    ############################################
    defaultlink <- list(D='log', g0='logit', lambda0='log', esa='log', a0='log',
                        sigma='log', sigmak='log', z='log', w='log', c='identity', pID='logit',
                        beta0='identity', beta1='neglog', sdS='log', b0='log', b1='neglog',
                        pmix='logit', cuerate='log', cut='identity')
    if (anycount) defaultlink$g0 <- 'log'
    link <- replace (defaultlink, names(link), link)
    if (details$scaleg0) link$g0 <- 'log'  ## Force log link in this case as no longer 0-1
    if (!(detector(traps(capthist))=='cue')) link$cuerate <- NULL
    link[!(names(link) %in% c(fnames,pnames))] <- NULL

    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################

    memo ('Preparing detection design matrices', details$trace)
    design <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                              dframe)
    design0 <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                               dframe, naive = T, bygroup = !CL)

    ############################################
    # Prepare density design matrix
    ############################################
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

    ############################################
    # Parameter mapping (general)
    ############################################
    if (!is.null(design$designMatrices))
        np <- sapply(design$designMatrices, ncol)
    else {
        np <- c(detectpar = 0)
    }
    np <- c(D = nDensityParameters, np)
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)[np>0]
    if (!D.modelled) parindx$D <- NULL

    ############################################
    # Optionally start from previous fit
    ############################################

    if (inherits(start, 'secr')) {
        ## use 'mapbeta' from score.test.R
        start <- mapbeta(start$parindx, parindx, coef(start)$beta, NULL)
    }

    ############################################
    # send data to worker processes
    # do it once, not each evaluation
    ############################################
    if (ncores > 1) {
        clusterExport(clust, c("capthist", "mask", "groups",
            "design","design0", "detectfn"), envir = environment())
    }

    ############################################
    # Single evaluation option
    ############################################
    .localstuff$iter <- 0
    if (details$LLonly) {
      if (is.null(start))
          stop ("must provide transformed parameter values in 'start'")
      if (!is.null(q))
          stop ("not for mark-resight")

      LL <- - secr.loglikfn (beta = start,
                       parindx    = parindx,
                       link       = link,
                       fixedpar   = fixed,
                       designD    = designD,
                       design     = design,
                       design0    = design0,
                       capthist   = capthist,
                       mask       = mask,
                       detectfn   = detectfn,
                       CL         = CL,
                       hcov       = hcov,
                       groups     = groups,
                       details    = details,
                       logmult    = TRUE,     ## add if possible
                       ncores     = ncores,
                       clust      = clust
                       )

      return(c(logLik=LL))
    }
    ############################################
    # Start values (model-specific)
    # 'start' is vector of beta values (i.e. transformed)
    ############################################
    if (is.null(start)) {
        start3 <- list(D = NA, g0 = NA, sigma = NA)
        ## not for signal attenuation
        if (!(detectfn %in% c(9,10,11,12,13)) && !anypoly && !anytrans) {
            memo('Finding initial parameter values...', details$trace)
            ## autoini uses default buffer dbar * 4
            if (MS)
                ## Using session 1, but this can be risky
                start3 <- autoini (capthist[[1]], mask[[1]],
                                   binomN = details$binomN,
                                   ignoreusage = details$ignoreusage)
            else
                start3 <- autoini (capthist, mask,
                                   binomN = details$binomN,
                                   ignoreusage = details$ignoreusage)

            if (any(is.na(unlist(start3)))) {
                warning ("'secr.fit' failed because initial values not found",
                         " (data sparse?); specify transformed values in 'start'")
                return (list(call = cl, fit = NULL))
            }
            if (details$unmash & !CL) {
                nmash <- attr(capthist[[1]], 'n.mash')
                if (!is.null(nmash)) {
                    n.clust <- length(nmash)
                    start3$D <- start3$D / n.clust
                }
            }
            nms <- c('D', 'g0', 'sigma')
            nms <- paste(nms, '=', round(unlist(start3),5))
            memo(paste('Initial values ', paste(nms, collapse=', ')),
                details$trace)
        }
        else warning ("using default starting values")

        #--------------------------------------------------------------
        # assemble start vector
        rpsv <- unlist(RPSV(capthist))[1]
        n <- ifelse (ms(capthist), nrow(capthist[[1]]), nrow(capthist))
        default <- list(
            D     = ifelse (is.na(start3$D), 1, start3$D),
            g0    = ifelse (is.na(start3$g0), 0.1, start3$g0),
            lambda0 = -log(1-ifelse (is.na(start3$g0), 0.1, start3$g0)),
            sigma = ifelse (is.na(start3$sigma), rpsv, start3$sigma),
            z     = 5,
            w     = 10,
            pID   = 0.7,
            beta0 = details$cutval + 30,
            beta1 = -0.2,
            sdS   = 2,
            b0    = 2,      ## changed from 15 2010-11-01
            b1    = -0.1,
            pmix  = 0,      ## superceded below

# esa   = ifelse (is.na(start3$g0), 10, start3$g0 * start3$sigma^2 *
#     ndetector(traps(capthist))[1] / 2500),
# replaced 2014-03-24
            esa   = ifelse (is.na(start3$D), n / 5, n / start3$D),
            a0    = ifelse (is.na(start3$g0), 0.1 * rpsv^2, start3$g0 *
                start3$sigma^2) / 10000 * 2 * pi
        )

        ## moved from inside autoini block 2013-07-19
        if (details$scaleg0 & anycount)
            stop ("'scaleg0' not compatible with count detectors")
        ## next two stmts must be this order (g0 then sigma)
        if (details$scaleg0) default$g0 <- default$g0 * default$sigma^2
        if (details$scaleg0) default$lambda0 <- default$lambda0 * default$sigma^2
        if (details$scalesigma) default$sigma <- default$sigma * default$D^0.5
        if (details$param %in% 4:6) {
            default$sigmak <- default$sigma * default$D^0.5
            default$c <- 0 ## but problems if take log(c)
        }

        if (detectfn %in% c(6)) {
            default$w <- default$sigma
            default$sigma <- default$sigma/2
        }
        if (detectfn %in% c(7)) {
            default$z <- default$sigma/5
        }
        if (detectfn %in% c(8, 18)) {
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
            default$lambda0 <- -log(1-default$g0)
            if (details$binomN > 1)
                default$g0 <- default$g0 / details$binomN
            if ((details$binomN == 1) & (detector(traps(capthist)) %in%
                 c('polygon','transect'))) {
                ## assume using usage for binomN
                usge <- usage(traps(capthist))
                default$g0 <- default$g0 / mean(usge[usge>0])
            }
            default$sigma <- RPSV(tempcapthist)
        }
        if (is.na(default$sigma)) default$sigma <- 20
        getdefault <- function (par) {
            transform (default[[par]], link[[par]])
        }

        start <- rep(0, NP)
        for ( i in 1:length(parindx) ) {
            # print (getdefault (names(model)[i]))
            start[parindx[[i]][1]] <- getdefault (names(model)[i])
        }
        if ((details$nmix>1) & !('pmix' %in% fnames))
            ## new starting values 2013-04-20
#            start[parindx[['pmix']]] <- clean.mlogit((1:nmix)-0.5)[-1]
        {
            start[parindx[['pmix']][1]] <- clean.mlogit((1:nmix)-0.5)[2]
        }
        if (detector(traps(capthist))=='cue')
            start <- c(start, log(3))    ## cuerate
        if (detectfn %in% c(12,13))
            start <- c(start, 46,3)    ## muN, sdN

        # D/ngrp when figure out where to calculate this

        ## if (!(is.null(q) | !nonID) & is.null(fixed$pID))
        ##     start[parindx$pID[1]] <- getdefault('pID')

        # start vector completed
        #--------------------------------------------------------------
    }

    ############################################
    ## ad hoc fix for experimental parameters
    ############################################
    nmiscparm <- 0
    if (detector(traps(capthist)) %in% c('cue'))
        nmiscparm <- 1
    if (detector(traps(capthist)) %in% c('signalnoise'))
        nmiscparm <- 2

    NP <- NP + nmiscparm
    stopifnot (length(start) == NP)

    ############################################
    # Fixed beta parameters
    ############################################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        ## drop unwanted betas; remember later to adjust parameter count
        start <- start[is.na(fb)]
        NP <- length(start)
    }
    ############################################
    # Variable names (general)
    ############################################
    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    realnames <- names(model)
    if (D.modelled) betanames <- c(paste('D', Dnames, sep='.'), betanames)
    betanames <- sub('..(Intercept))','',betanames)

    ############################################
    # Variable names (model-specific)
    ############################################

    if (detector(traps(capthist))=='cue') {
        betanames <- c(betanames, 'cuerate')
        realnames <- c(realnames, 'cuerate')
    }
    if (detectfn %in% c(12,13)) {
        betanames <- c(betanames, 'muN', 'sdN')
        realnames <- c(realnames, 'muN', 'sdN')
    }
    ## allow for fixed beta parameters
    if (!is.null(details$fixedbeta))
        betanames <- betanames[is.na(details$fixedbeta)]
    betaw <- max(max(nchar(betanames)),8)   # for 'trace' formatting

    ############################################
    # Maximize likelihood
    ############################################

    memo('Maximizing likelihood...', details$trace)
    if (details$trace)
        cat('Eval     Loglik', formatC(betanames, format='f', width=betaw), '\n')

    loglikefn <- secr.loglikfn
    if (!is.null(q))
        stop ("mark-resight option not operative")
        ## loglikefn <- MRsecr.loglikfn

    ## arguments always passed to loglikefn
    secrargs <- list(
                     link       = link,
                     fixedpar   = fixed,
                     parindx    = parindx,
                     capthist   = capthist,
                     mask       = mask,
                     CL         = CL,
                     detectfn   = detectfn,
                     designD    = designD,
                     design     = design,
                     design0    = design0,
                     hcov       = hcov,
                     groups     = groups,
                     details    = details,
                     logmult    = TRUE,     # add if possible
                     ncores     = ncores,
                     clust      = clust,
                     betaw      = betaw)    # for trace format

    ############################################
    ## calls for specific maximisation methods
    lcmethod <- tolower(method)
    ## 2013-04-21
    if (NP == 1) {
        lcmethod <- "optimise"
        signs <- c(-1,1) * sign(start)
        args <- list (f         = loglikefn,
            interval  = start * (1 + details$intwidth2 * signs))
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- try(do.call (optimise, args))
        if (inherits(this.fit, 'try-error'))
            warning ("univariate search for minimum failed")
        this.fit$par <- this.fit$minimum
        this.fit$value <- this.fit$objective
        if (details$hessian != "none")
            details$hessian <- "fdHess"
    }
    else if (lcmethod %in% c('newton-raphson')) {
        args <- list (p         = start,
                      f         = loglikefn,
                      hessian   = tolower(details$hessian)=='auto',
                      stepmax   = 10)
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (nlm, args)
        this.fit$par <- this.fit$estimate     # copy for uniformity
        this.fit$value <- this.fit$minimum    # copy for uniformity
        if (this.fit$code > 2)
            warning ("possible maximization error: nlm returned code ",
                this.fit$code, ". See ?nlm")
    }
    #-----------------------------------------------------------------
    else if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
                             "SANN", "Brent")) {
        args <- list(par     = start,
                     fn      = loglikefn,
                     hessian = tolower(details$hessian)=='auto',
                     method  = method)
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (optim, args)
        if (this.fit$convergence != 0)
            warning ("probable maximization error: optim returned convergence ",
                this.fit$convergence, ". See ?optim")
    }
    #-----------------------------------------------------------------
    # Hessian-only 2013-02-23
    else if (lcmethod %in% 'none') {
        memo ('Computing Hessian with fdHess in nlme', details$trace)
        loglikfn <- function (beta) {
            do.call(secr.loglikfn, c(list(beta=beta), secrargs))
        }
        grad.Hess <- nlme::fdHess(start, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
#        this.fit <- list (value = 0, par = start,
        # let's keep the loglik 2014-01-25
        this.fit <- list (value = loglikfn(start), par = start,
                          gradient = grad.Hess$gradient,
                          hessian = grad.Hess$Hessian)
        biasLimit <- NA   ## no bias assessment
    }
    else stop ("maximization method", method, "not recognised")
    ############################################################################

    this.fit$method <- method         ## remember which method we used...
    covar <- NULL
    N <- NULL
    if (this.fit$value > 1e9) {     ## failed
        # this.fit$beta[] <- NA
        this.fit$par[] <- NA
    }
    else {

    ############################################
    ## Variance-covariance matrix
    ############################################

        if (tolower(details$hessian)=='fdhess') {
                memo ('Computing Hessian with fdHess in nlme', details$trace)
                loglikfn <- function (beta) {
                    do.call(secr.loglikfn, c(list(beta=beta), secrargs))
                }
                grad.Hess <- nlme::fdHess(this.fit$par, fun = loglikfn,
                                          .relStep = 0.001, minAbsPar = 0.1)
                this.fit$hessian <- grad.Hess$Hessian
        }
        if (!is.null(this.fit$hessian)) {
            covar <- try(MASS::ginv(this.fit$hessian))
            if (inherits(covar, "try-error")) {
                warning ("could not invert Hessian to compute ",
                         "variance-covariance matrix")
                covar <- matrix(nrow = NP, ncol = NP)
            }
            else if (any(diag(covar)<0)) {
                warning ("at least one variance calculation failed ")
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

    ############################################
    ## form output list
    ############################################

    desc <- packageDescription("secr")  ## for version number

    ## ad hoc censoring of evaluated calls from do.call 2014-02-10
    ## wrong, and ugly...
    ## if (object.size(cl) > details$maxcallsize)
    ##    cl <- paste('secr.fit(',substring(as.character(cl)[2],1,details$maxcallsize),'...',sep='')

    output <- list (call = cl,
                  capthist = capthist,
                  mask = mask,
                  detectfn = detectfn,
                  CL = CL,
                  timecov = timecov,
                  sessioncov = sessioncov,
                  hcov = hcov,
                  groups = groups,
                  dframe = dframe,
                  design = design,
                  design0 = design0,
                  start = start,
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
                  N = N,
                  version = desc$Version,
                  starttime = starttime,
                  proctime = (proc.time() - ptm)[3]
             )

    class(output) <- 'secr'

    if (usebuffer & !is.na(biasLimit)) {
        test <- try(bufferbiascheck(output, buffer, biasLimit))
    }

    memo(paste('Completed in ', round(output$proctime,2), ' seconds at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"),
        sep=''), details$trace)

    if (ncores > 1) {
        stopCluster(clust)
    }

    output

}
################################################################################
