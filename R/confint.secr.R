############################################################################################
## package 'secr'
## confint.secr.R
## last changed 2009 06 11, 2009 07 16 2009 10 20
## could be sped up by adopting Venzon & Moolgavkar algorithm e.g. in Bhat package 

############################################################################################

confint.secr <- function (object, parm, level = 0.95, newdata = NULL, tracelevel = 1,
    tol = 0.0001, ...) {

    ## profile likelihood interval for estimated secr parameters
    ## cf confint.glm in MASS

    ## interpret character 'parm' as real parameter
    ## interpret numeric 'parm' as beta parameter

    ## case 1 - real parameter not 1:1 beta so require lagrange
    ## case 2 - real parameter but model ~1 so 1:1 beta
    ## case 3 - beta parameter

    #---------------------------------------------------------------------------------------

    profileInterval <- function (parm, ...) {

        predicted <- function (beta) {
            temp <- secr.lpredictor (newdata = newdata, model = object$model[[parm]], 
                indx = object$parindx[[parm]], beta = beta)[1,'estimate']
            untransform(temp, object$link[[parm]])
        }

        #######################
        ## case 1 - Lagrange

        profileLL.lagrange <- function (gamma, parm) {
            lagrange <- function (beta2, gamma) {
                ## return quantity to be maximized
                templl <- secr.loglikfn (
                    beta       = beta2,
                    link       = object$link,
                    fixed      = object$fixed,
                    parindx    = object$parindx,
                    capthist   = object$capthist,
                    mask       = object$mask,
                    CL         = object$CL,
                    detectfn   = object$detectfn,
                    designD    = D.designmatrix,  
                    design     = object$design,
                    design0    = object$design0,
                    groups     = object$groups,
                    details    = details,
                    logmult    = logmult,
                    betaw      = max(max(nchar(object$betanames)),8))
                templl - gamma * predicted (beta2)   
            }

            ## maximize for fixed gamma (equivalent to fixed 'parm')
            lagrange.fit <- nlm (p = object$fit$par, f = lagrange, gamma = gamma, hessian = FALSE)  
            .localstuff$beta <- lagrange.fit$estimate
            lp <- - secr.loglikfn (
                beta       = lagrange.fit$estimate,
                link       = object$link,
                fixed      = object$fixed,
                parindx    = object$parindx,
                capthist   = object$capthist,
                mask       = object$mask,
                CL         = object$CL,
                detectfn   = object$detectfn,
                designD    = D.designmatrix,  
                design     = object$design,
                design0    = object$design0,
                groups     = object$groups,
                details    = details,
                logmult    = logmult,
                betaw      = max(max(nchar(object$betanames)),8))
            ## cat ('gamma ', gamma, '  coef ', lagrange.fit$estimate, '  lp ', lp, 
            ##     '  lp - targetLL ', lp-targetLL, '\n')
            lp - targetLL
        }

        #######################
        ## cases 2,3 - fix beta

        profileLL3 <- function (x, parm) {
            ## fix required beta parameter to x and evaluate discrepancy
            fb <- object$details$fixedbeta
            if (is.null(fb)) fb <- rep(NA, np)
            fb[parm] <- x
            details$fixedbeta <- fb
            fit <- secr.fit (capthist = object$capthist, mask = object$mask, 
                buffer = object$buffer, CL = object$CL, detectfn = object$detectfn,   
                start = object$start, 
                link = object$link, fixed = object$fixed,  model = object$model, 
                timecov = object$timecov, sessioncov = object$sessioncov, 
                groups = object$groups, dframe = object$dframe, 
                details = details, verify = FALSE, ...)$fit 
            - fit$value - targetLL
        }

        #######################
        getlimit <- function (start, step, trialvals) {
            OK <- FALSE
            bound.0 <- start
            f.0 <- profileLL (start, parm)
            for (i in trialvals) { 
               bound.1 <- start + step * i 
               f.1 <- profileLL (bound.1, parm)
               if (prod(sign(c(f.0, f.1)))<0) {  ## different
                   OK <- TRUE
                   break
               }
               else {
                   bound.0 <- bound.1
                   f.0 <- f.1
               }                
            }
            if (!OK) {
                if (is.character(parm))
                    stop(paste('did not find root within', max(trialvals), 'units of MLE'))
                else
                    stop(paste('did not find root within', max(trialvals), 'SE of MLE'))
            }
            if (step < 0)
                c(lower = bound.1, upper = bound.0, f.lower = f.1, f.upper = f.0)
            else
                c(lower = bound.0, upper = bound.1, f.lower = f.0, f.upper = f.1)
        }

        #######################
        
        if (is.numeric(parm)) 
            profileLL <- profileLL3
        else 
            profileLL <- profileLL.lagrange

        estimate <- coef(object, alpha=1-level)[parm,]

        if (is.numeric(parm))  {
            startlow <- getlimit(estimate$beta, -2 * estimate$SE.beta, c(1,2,4,8))
            startupp <- getlimit(estimate$beta, +2 * estimate$SE.beta, c(1,2,4,8))
        }
        else {
            startlow <- getlimit (0, -1, c(2,5,40,200,1000))
            startupp <- getlimit (0, +1, c(2,5,40,200,1000))
        } 

        temproot <- uniroot (profileLL, 
            parm    = parm,
            lower   = startlow['lower'],
            upper   = startlow['upper'],
            f.lower = startlow['f.lower'],
            f.upper = startlow['f.upper'],
            tol     = tol)

        if (is.numeric(parm))  
            LCL <- temproot$root
        else 
            LCL <- predicted (.localstuff$beta)

        temproot <- uniroot (profileLL, 
            parm    = parm,
            lower   = startupp['lower'],
            upper   = startupp['upper'],
            f.lower = startupp['f.lower'],
            f.upper = startupp['f.upper'],
            tol     = tol)

        if (is.numeric(parm))  
            UCL <- temproot$root
        else 
            UCL <- predicted (.localstuff$beta)

        c(LCL, UCL)
    }
    # end of profileInterval
    #---------------------------------------------------------------------------------------

    memo ('Profile likelihood interval(s)...', tracelevel > 0)

    if (!inherits(object, 'secr')) stop ('Requires secr object')
    np <- length(object$betanames)  ## number of beta parameters

    ## case 1 - real parameter not 1:1 beta so require lagrange
    ## case 2 - real parameter but model ~1 so 1:1 beta
    ## case 3 - beta parameter

    if (is.character(parm)) {
        OK <- (parm %in% object$realnames)
        if (any(!OK)) 
            stop ('requested parameter(s) not in model')
        case <- ifelse (object$model[parm] != ~1, 1, 2)
        if (any(case==1)) {
            if (is.null(newdata))
                newdata <- secr.make.newdata (object)[1,, drop = FALSE]  ## default base levels
            logmult <- logmultinom(object$capthist, group.factor(object$capthist, object$groups))

            ## reconstruct density design matrix
            D.modelled <- !object$CL & is.null(object$fixed$D)
            if (!D.modelled) {
                D.designmatrix <- matrix(nr=0, nc=0)
                attr(D.designmatrix, 'dimD') <- NA
            }
            else {
                grps  <- group.levels(object$capthist, object$groups)
                temp <- D.designdata( object$mask, object$model$D, grps, 
                    session(object$capthist), object$sessioncov)
                D.designmatrix <- model.matrix(object$model$D, temp) 
                attr(D.designmatrix, 'dimD') <- attr(temp, 'dimD')
            }
        }
    }
    else { 
        if (any ((parm<1) | (parm>np))) 
            stop ('invalid beta parameter number')
        case <- rep(3, np)
    }
    #---------------------------------------------------------------------------------------

    targetLL <- - object$fit$value -  qchisq(level,1)/2   # -1.92 for 95% interval
    details <- replace (object$details, 'hessian', FALSE)  ## no need for vcov matrix
    details$trace <- tracelevel > 1
    out <- matrix(nr=length(parm), nc=2)

    for (i in 1:length(parm)) {
        parmn <- ifelse (case[i] == 2, match(parm[i], object$betanames), parm[i])        
        out[i,] <- profileInterval(parmn)   ## character value if 'real'
        if (case[i] == 2) out[i,] <- untransform(out[i,], object$link[[parm[i]]])
    }
    if (case == 3) 
        dimnames(out) <- list(object$betanames[parm], c('beta.lcl', 'beta.ucl'))
    else
        dimnames(out) <- list(parm, c('lcl', 'ucl'))
    out
}
############################################################################################
