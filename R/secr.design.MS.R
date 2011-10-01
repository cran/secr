###############################################################################
## package 'secr'
## secr.design.MS.R
## 2009 08 20 disable bk, Bk
## 2009 10 08 sighting
## 2009 12 10 mixtures
## 2009 12 13 mixtures pmix ~ h2

################################################################################
## source ('d:\\density secr 1.5\\secr\\R\\secr.design.MS.R')
################################################################################

secr.design.MS <- function (capthist, models, timecov = NULL, sessioncov = NULL,
    groups = NULL, dframe = NULL, naive = FALSE, bygroup = FALSE, ...) {

## Generate design matrix, reduced parameter array, and parameter index array (PIA)
## for detection function parameters
## 'capthist' must be of class 'capthist' or 'list'
##
## uses pad1 to pad session-specific covar to constant length with first value,
## defined in 'functions.R'

    findvars.MS <- function (cov, vars, dimcov) {
        ## function to add covariates to a design data frame 'dframe'
        ## cov may be a dataframe or list of dataframes, one per session (R > 1),
        ## if list, then require predictors to appear in all sessions
        ## uses pad1 and insertdim from functions.R
        ## NOT to be used to add group variables
        ## Does not yet standardize numeric covariates if (!is.factor(vals)) vals <- stdfn(vals)

        if (is.null(cov) | (length(cov)==0) | (length(vars)==0)) return()
        else {
            found <- ''
            if (!is.data.frame(cov)) {   ## therefore assume is a list
                if (!is.list(cov) | (R==1))
                    stop ("irregular covariates; check multisession structure")
                covnames <- lapply(cov, names)
                varincov <- sapply(covnames, function(nam) vars %in% nam)
                if (length(vars)>1) found <- vars[apply(varincov,1,all)]
                else found <- vars[all(varincov)]

                for (variable in found) {
                    vals <- unlist(lapply(cov, function(x) pad1(x[,variable],
                        dims[dimcov])))
                    if (any(is.na(vals)))
                        stop ("covariate missing values not allowed")

                    dframe[,variable] <<- insertdim (vals, dimcov, dims)
                }
            }
            else
            {
                found <- names(cov) %in% vars
                if (is.data.frame(cov) & any(found)) {
                    found <- names(cov)[found]
                    values <- as.data.frame(cov[,found])
                    names(values) <- found
                    if (length(values)>0) {
                        for (variable in found) {
                            vals <- values[,variable]
                            dframe[,variable] <<- insertdim (vals, dimcov, dims)
                        }
                    }
                }
            }
            vars <<- vars[!(vars %in% found)]
        }
    }

    # bygroup = T results in one row per group instead of one row per individual
    # This setting is used for 'naive' table
    # groups is a vector of factor names whose intersection defines group
    # groups only defined for CL = T
    # use of 'g' requires valid groups definition
    # grouping variables are also added individually to dframe

    models$D <- NULL                          # drop density model
    npar     <- length(models)                # real parameters
    parnames <- names(models)                 # typically c('g0', 'sigma', 'z')
    vars     <- unique (unlist(sapply (models, all.vars)))
    vars     <- vars[!(vars %in% groups)]     # groups treated separately
    nmix     <- get.nmix(models)

##    if (sum(c('b','B','bk','Bk') %in% vars) > 1)
    if (sum(c('b','B') %in% vars) > 1)
    stop ("model should not use more than one type of behavioural response")
    trps    <- traps(capthist)                 # session-specific trap array
    used    <- usage(trps)                     # session-specific usage
    zcov    <- covariates(capthist)            # session-specific individual covariates
    trapcov <- covariates(trps)                # session-specific trap covariates
    if (('g' %in% vars) & is.null(groups))
        stop ("requires valid 'groups' covariate")
    grps    <- group.levels(capthist,groups)
    ngrp    <- max(1,length(grps))

    ## 'session-specific' list if MS
    MS   <- inherits(capthist, 'list') # logical for multi-session
    sessionlevels <- session(capthist)
    if (is.null(sessionlevels)) sessionlevels <- '1'
    getnk <- function (object) {
        ifelse (detector(object) %in% .localstuff$polydetectors,
            length(table(polyID(object))),
            nrow(object))
    }

    if (MS) {
        R <- length(capthist)
        n <- ifelse (bygroup, ngrp, max(sapply(capthist, nrow))) # max over sessions
        S <- max(sapply(capthist, ncol))                         # max over sessions
        K <- max(sapply(trps, getnk))                            # max over sessions
    }
    else {
        R <- 1
        n <- ifelse (bygroup, ngrp, nrow(capthist))
        S <- ncol(capthist)
        K <- getnk(trps)
    }
    ## cover unmarked case
    if (n == 0) n <- 1

    ## timecov may be a vector or a dataframe or a list of vectors or a list of data frames
    ## conform to either dataframe or list of dataframes (1 per session)
    if (!is.data.frame(timecov)) {
        if (is.list(timecov)) {
            if (length(timecov) != R)
                stop ("wrong number of sessions in 'timecov' list")
            timecov <- lapply(timecov, as.data.frame)
        }
        else timecov <- as.data.frame(timecov)
    }

    if (!is.null(sessioncov)) {
        sessioncov <- as.data.frame(sessioncov)
        if (nrow(sessioncov) != R)
            stop("number of rows in 'sessioncov' should equal ",
                 "number of sessions")
    }

    dims <- c(R,n,S,K,nmix)  # 'virtual' dimensions
    autovars <- c('session','Session','g','t','T','b','B','kcov','tcov', 'sighting', 'h2', 'h3')
##    autovars <- c('session','Session','g','t','T','b','B','bk','Bk',
##                  'kcov','tcov', 'sighting', 'h2', 'h3')
    if (is.null(dframe)) dframe <- data.frame(dummy=rep(1, R*n*S*K*nmix))

## BUG FIXED 2009 12 22
## Default action of following was to sort levels
## but these are allowed to be in arbitrary order
##    dframe$session <- factor( insertdim (sessionlevels, 1, dims))

    dframe$session <- factor( insertdim (sessionlevels, 1, dims), levels = sessionlevels)
    dframe$Session <- insertdim (0:(R-1), 1, dims)

    #--------------------------------------------------------------------------
    if ('t' %in% vars)
        dframe$t <- factor( insertdim (1:S, 3, dims) )

    #--------------------------------------------------------------------------
    ## 2009 10 07 mark-resight
    q <- attr(capthist, 'q')
    if (any((q > 0) & (q < S))) {
        ## marking  on occasions 1:q, sighting on occasions (q+1):S
        ## following line must be revised for vector q (varying by session)
        dframe$sighting <- factor( insertdim (rep(0:1, c(q, S-q)), 3, dims) )
    }

    #--------------------------------------------------------------------------
    if (!is.null(groups)) {

        ################
        # add g factor

        if (bygroup)
            gvar <- factor(grps)
        else {
            gvar <- group.factor(capthist, groups)          # list if MS
            if (MS) gvar <- lapply(gvar, pad1, n)           # constant length
        }
        # by animal within session
        dframe$g <- insertdim ( unlist(gvar), c(2,1), dims)  ## unlist works on factors, too

        #################################
        # also add separate group factors

        # if bygroup, just use group levels for each session
        # else get group membership from covariates(capthist) for each session,
        # and pad to max(n) if needed

        if (bygroup) {
            MSlevels <- function (facname) {
                ## assume all sessions the same without checking
                if (inherits(capthist, 'list')) levels(zcov[[1]][,facname])
                else levels(zcov[,facname])
            }

            ## list with levels for each grouping factor named in 'groups'
            grouplist <- lapply(groups, MSlevels)
            grouping <- expand.grid(grouplist) ## one column per group
            names(grouping) <- groups
            for (i in groups) insertdim(grouping[,i], c(2,1), dims)  ## all sessions the same
        }
        else {
            for (i in groups) {
                if (MS) {
                    grouping <- lapply(zcov, function(x) x[,i])
                    grouping <- unlist(lapply(grouping, pad1, n))
                }
                else grouping <- zcov[,i]
                insertdim(grouping, c(2,1), dims)
            }
        }
    }
    #--------------------------------------------------------------------------

    ## general for behavioural response fields

    ## assume sessions same type
    #  prox <- detector(trps)[[1]] %in% c('proximity', 'count', 'polygon','transect')
    ## 2011-09-26
    prox <- detector(trps)[[1]] %in% .localstuff$detectors3D
    makeb <- function (caphist) {      ## global response
        if (prox) {
           temp0 <- apply(caphist, 1:2, sum)
           t(apply(temp0, 1, prevcapt))
        }
        else t(apply(abs(caphist),1, prevcapt))
    }

##    makebk <- function (caphist) {     ## trap-specific response
##        if (prox) {
##            temp <- apply(caphist, c(1,3), prevcapt)
##        }
##        else {
##            temp0 <- array (dim=c(n,S,K))
##            temp0[,,] <- 0
##            indices <- cbind(
##                as.numeric(row(caphist)),
##                as.numeric(col(caphist)),
##                as.numeric(abs(caphist)))
##            temp0[indices] <- 1
##            temp <- apply(temp0, c(1,3), prevcapt)
##        }
##        aperm(temp, c(2,1,3))
##    }

    #--------------------------------------------------------------------------

    if ('b' %in% vars) {
        if (naive) dframe$b <- rep(0, R*n*S*K)
        else {
            prevcapt <- function(x) c(F, cumsum(x[-S])>0)
            if (MS) {
                temp <- lapply(capthist, makeb)
                temp <- lapply(temp, padarray, c(n,S))
                temp <- unlist(temp)
            }
            else temp <- makeb(capthist)  # one session
            dframe$b <- insertdim (as.vector(temp), c(2,3,1), dims)
        }
    }

    #------------------------------------------------

    if ('B' %in% vars) {
        if (naive) dframe$B <- rep(0, R*n*S*K)
        else {
            prevcapt <- function(x) c(F, x[-S]>0)
            if (MS) {
                temp <- lapply(capthist, makeb)
                temp <- lapply(temp, padarray, c(n,S))
            }
            else temp <- makeb(capthist)  # one session
            dframe$B <- insertdim (as.vector(unlist(temp)), c(2,3,1), dims)
       }
    }

    #------------------------------------------------

##    if ('bk' %in% vars) {
##        if (naive) dframe$bk <- rep(0, R*n*S*K)
##        else {
##            prevcapt <- function(x) c(F, cumsum(x[-S])>0)
##            if (MS) {
##                temp <- lapply(capthist, makebk)
##                temp <- lapply(temp, padarray, c(n,S,K))
##            }
##            else temp <- makebk(capthist)  # one session
##            dframe$bk <- insertdim(as.vector(unlist(temp)), c(2,3,4,1), dims)
##        }
##    }

    #------------------------------------------------
##    if ('Bk' %in% vars) {
##        if (naive) dframe$Bk <- rep(0, R*n*S*K)
##        else {
##            prevcapt <- function(x) c(F, x[-S]>0)
##            if (MS) {
##                temp <- lapply(capthist, makebk)
##                temp <- lapply(temp, padarray, c(n,S,K))
##            }
##            else temp <- makebk(capthist)  # one session
##            dframe$Bk <- insertdim(as.vector(unlist(temp)), c(2,3,4,1), dims)
##       }
##    }
    #--------------------------------------------------------------------------

    if ('tcov' %in% vars) {
        if (is.null(timecov))
            stop ("requires valid time covariate 'timecov'")
        if (length(unlist(timecov)) > S)                                             ### CHECK
            warning ("length of 'timecov' exceeds number of occasions")
        if (is.data.frame(timecov)) {
            if (nrow(timecov) < S)
                stop ("requires valid time covariate 'timecov'")
            timecov <- timecov[,1,drop=F]  ## retain only first
        }
        else {
            if (any(sapply(timecov, nrow) < S))
                stop ("requires valid time covariate 'timecov'")
            timecov <- lapply (timecov, function(x) pad1(x[,1], S))
        }
        dframe$tcov <- insertdim (unlist(timecov), c(3,1), dims)
    }

    if ('T' %in% vars) {
        dframe$T <- insertdim (0:(S-1), c(3,1), dims)
    }

    #--------------------------------------------------------------------------
    if ('kcov' %in% vars) {
        if (is.null(trapcov))
            stop ("model uses trap covariates, but valid covariate ",
                  "data not found")
        if (is.data.frame(trapcov)) trapcov <- trapcov[,1,drop=F]  ## retain only first
        else trapcov <- lapply (trapcov, function(x) pad1(x[,1], K))
        dframe$kcov <- insertdim(unlist(trapcov), c(4,1), dims)
    }

    #--------------------------------------------------------------------------

    ## h2 or h3
    if (nmix > 1) {
        mixture <- paste('h',nmix,sep='')
        dframe[,mixture] <- insertdim(factor(1:nmix), 5, dims)
    }

    #--------------------------------------------------------------------------
    ## all autovars should have now been dealt with
    vars <- vars[!vars %in% autovars]

    #--------------------------------------------------------------------------
    # add zcov, sessioncov, trapcov, timecov

    if ((length(vars)>1) & (R>1))
        warning ("implementation of user covariates with multi-session ",
                 "data is not fully tested")

    if (!bygroup) findvars.MS (zcov, vars, 2)   ## CL only
    ## possible alternative for consideration 2009 09 25
    ## findvars.MS (zcov, vars, 2)

    findvars.MS (sessioncov, vars, 1)
    findvars.MS (timecov, vars, 3)      ## session-specific list
    findvars.MS (trapcov, vars, 4)      ## session-specific list

    if (length(vars)>0) {
        if (!is.null(zcov)) {
            if (is.data.frame(zcov))
                znames <- names(zcov)
            else
                znames <- unlist(lapply(zcov, names))
            if (any (vars %in% znames))
                stop ("seems you are trying to use individual covariates ",
                      " in a full-likelihood model")
        }
        stop ("covariate(s) ", paste(vars,collapse=","), " not found")
    }

    #------------------------------------------------
    # DOES NOT YET STANDARDIZE COVARIATES (cf Dmat)

    make.designmatrix <- function (formula, prefix, ...) {
     # combine formula and dframe to generate design matrix
        if (is.null(formula)) {
            list (model = NULL, index = rep(1,R*n*S*K*nmix))
        }
        else {
            tempmat <- model.matrix(formula, dframe, ...)
            ## drop pmix beta0 column from design matrix
            if (prefix=='pmix') tempmat <- tempmat[,-1,drop=FALSE]
            temp <- make.lookup (tempmat)   # retain unique rows
            list (model=temp$lookup, index=temp$index)
        }
    }

    ## replace NA with dummy value '0' to stop model.matrix dropping rows added as padding
    ## might be more efficient to allow it to drop and then match later
    dframe[is.na(dframe)] <- 0

    # list with one component per real parameter
    # each of these is a list with components 'model' and 'index'

    designMatrices <- sapply (1:length(models), simplify=FALSE,
        function (x) make.designmatrix(models[[x]], names(models[x])))
    names(designMatrices) <- names(models)

    ## dim(indices) = c(R*n*S*K*nmix, npar)
    indices <- sapply (designMatrices, function(x) x$index)

    indices <- matrix(unlist(indices), ncol = npar)

    # retain just the 'model' components of 'designMatrices'
    designMatrices <- lapply (designMatrices, function(x)x$model )

    # prefix column names in 'designMatrices' with parameter name
    for (i in 1:npar)
        colnames(designMatrices[[i]]) <- paste (parnames[i], '.',
            colnames(designMatrices[[i]]), sep='')

    # repackage indices to define unique combinations of parameters
    indices2 <- make.lookup(indices)

    #--------------------------------------------------------------------
    # PIA = Parameter Index Array
    #       index to row of parameterTable for a given R,n,S,K,nmix
    # dim(parameterTable) = c(uniqueparcomb, npar)
    #       index to row of designMatrix for each real parameter
    #--------------------------------------------------------------------

    PIA <- array(indices2$index, dim=c(R,n,S,K,nmix))

    parameterTable <- indices2$lookup

    colnames(parameterTable) <- parnames

    #--------------------------------------------------------------------
    # Zero the index of trap+time pairs that were 'not set'
    # the external C code checks for this and sets p(detection) to zero
    #--------------------------------------------------------------------

    # 'used' is list if MS

## OOPS - TOO TIRED TO THINK THROUGH THE IMPLEMENTATION OF USE WITH MIXTURES
## UNSURE FROM HERE 2009 12 10

    if ((!is.null(used)) & (length(used)>0)) {
        allused <- unlist(used)
        if (!is.null(allused))
        if (any(!allused)) {
            if (!MS) {
                PIA[1, , , ,] <- PIA[1, , , ,] * rep(rep(t(used),rep(n,S*K)),nmix)
            }
            else for (r in 1:R) {
                use <- array (0, dim=c(S,K))
                temp <- t(used[[r]])  # dim (S',K')
                use[1:nrow(temp), 1:ncol(temp)] <- temp  # padding
                PIA[r, , , ,] <- PIA[r, , , ,] * rep(rep(use,rep(n,S*K)),nmix)
            }
        }
    }

    #--------------------------------------------------------------------

    list(designMatrices = designMatrices, parameterTable = parameterTable, PIA = PIA, R = R)
}
############################################################################################


