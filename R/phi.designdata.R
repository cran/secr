###############################################################################
## package 'secr'
## phi.designdata.R
## 2011 11 30
## Initially only 'phi', but potentially other turnover parameters

################################################################################
## source ('d:\\density secr 2.3\\secr\\R\\phi.designdata.R')
################################################################################

phi.designdata <- function (capthist, models, intervalcov = NULL, sessioncov = NULL,
    groups = NULL, dframe = NULL, bygroup = FALSE, intervals, ...) {

## Generate design matrix, reduced parameter array, and parameter index array (PIA)
## for survival parameter
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
        ## Does not yet standardize numeric covariates if (!is.factor(vals)) vals <- scale(vals)

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
    # groups only defined for CL = FALSE  [corrected 2011-11-27]
    # use of 'g' requires valid groups definition
    # grouping variables are also added individually to dframe

    if (is.null(intervals)) {
        return(list(designMatrices = list(phi = matrix(0,nrow=0,ncol=1)),
            parameterTable = matrix(0,nrow=0,ncol=1),
            PIA = array(0,dim=c(1,1,1,1))))
    }


    models   <- models['phi']
    npar     <- length(models)                # real parameters
    parnames <- names(models)                 # typically c('phi')
    vars     <- unique (unlist(sapply (models, all.vars)))
    vars     <- vars[!(vars %in% groups)]     # groups treated separately
    nmix     <- get.nmix(models)

    trps    <- traps(capthist)                 # session-specific trap array
    zcov    <- covariates(capthist)            # session-specific individual covariates
    if (('g' %in% vars) & is.null(groups))
        stop ("requires valid 'groups' covariate")
    grouplevels  <- group.levels(capthist,groups)
    ngrp    <- max(1,length(grouplevels))

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
        gaps <- lapply(intervals, function (x) x[x>0])
        j <- lapply(intervals, function(x) c(1,cumsum(x>0)+1))   # s to j
        J <- max(sapply(gaps, length))                           # max over sessions
    }
    else {
        R <- 1
        n <- ifelse (bygroup, ngrp, nrow(capthist))
        gaps <- intervals[intervals>0]
        j <- c(1,cumsum(intervals>0))                            # s to j
        J <- length(gaps)
    }
    ## cover unmarked case
    if (n == 0) n <- 1

    #--------------------------------------------------------------------------
    # intervalcov may be a vector or a dataframe or a list of vectors or a list of data frames
    # conform to either dataframe or list of dataframes (1 per session)
    if (!is.null(intervalcov)) {
        if (!is.data.frame(intervalcov)) {
            if (is.list(intervalcov)) {
                if (length(intervalcov) != R)
                    stop ("wrong number of sessions in 'intervalcov' list")
                intervalcov <- lapply(intervalcov, as.data.frame)
            }
            else intervalcov <- as.data.frame(intervalcov)
        }
    }

    #--------------------------------------------------------------------------
    # session covariates
    if (!is.null(sessioncov)) {
        sessioncov <- as.data.frame(sessioncov)
        if (nrow(sessioncov) != R)
            stop("number of rows in 'sessioncov' should equal ",
                 "number of sessions")
    }

    #--------------------------------------------------------------------------
    dims <- c(R,n,J,nmix)  # 'virtual' dimensions
    dframenrow <- prod(dims)      # number of rows
    autovars <- c('g','b','B','interval', 'h2', 'h3')
    #--------------------------------------------------------------------------
    # user-specified dframe
    if (is.null(dframe)) {
        dframe <- data.frame(dummy=rep(1, dframenrow))
        dframevars <- ""
    }
    else {
        ## special treatment for design0:
        ## dframe must be thinned by rejecting rows for repeats from each group
        tempn <- n
        if (bygroup) {
            OK <- rep(FALSE, nrow(dframe) )
            tempn <- length(OK)/(R*J*nmix)
            dim(OK) <- c(R,tempn,J,nmix)
            if (is.null(groups)) {
                firstofgroup <- 1   ## no groups
                OK[,firstofgroup,,,] <- TRUE
                dframe <- dframe[OK,,drop=FALSE]
            }
            else
                stop ("dframe specification does not work with grouping at present, sorry")
        }
        if (nrow(dframe) !=  dframenrow )
            stop ("dframe should have ", R*tempn*J*nmix, " rows ( R*n*J*nmix )")
        dframevars <- names(dframe)
    }
    #--------------------------------------------------------------------------

    dframe$interval <- factor( insertdim (1:J, 3, dims))

    #--------------------------------------------------------------------------
    if (!is.null(groups)) {

        ################
        # add g factor

        if (bygroup)
            gvar <- factor(grouplevels)
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
            for (i in groups)
                dframe[,i] <- insertdim(grouping[,i], c(2,1), dims)  ## all sessions the same
        }
        else {
            for (i in groups) {
                if (MS) {
                    grouping <- lapply(zcov, function(x) x[,i])
                    grouping <- unlist(lapply(grouping, pad1, n))
                }
                else grouping <- zcov[,i]
                dframe[,i] <- insertdim(grouping, c(2,1), dims)
            }
        }
    }
    #--------------------------------------------------------------------------

    ## behavioural response fields
    ## assume sessions are same type   ******************************
    prox <- detector(trps)[[1]] %in% .localstuff$detectors3D
    makeb <- function (caphist) {      ## global response
        ## convert to 2-D CH
        if (prox)
            temp0 <- apply(abs(caphist), 1:2, sum)
        else
            temp0 <- abs(caphist)
        ## convert to n x (J+1) CH
        temp0 <- t(apply(temp0, 1, function(x) (1:(J+1)) %in% j[x>0])) * 1
        t(apply(temp0, 1, prevcapt))
    }
    #--------------------------------------------------------------------------

    if ('b' %in% vars) {
        prevcapt <- function(x) cumsum(x[-(J+1)])>0
        if (MS) {
            temp <- lapply(capthist, makeb)
            temp <- lapply(temp, padarray, c(n,J))
            temp <- unlist(temp)
        }
        else temp <- makeb(capthist)  # one session
        dframe$b <- insertdim (as.vector(temp), c(2,3,1), dims)
    }

    #------------------------------------------------

    if ('B' %in% vars) {
        prevcapt <- function(x) x[-(J+1)]>0
        if (MS) {
            temp <- lapply(capthist, makeb)
            temp <- lapply(temp, padarray, c(n,J))
        }
        else temp <- makeb(capthist)  # one session
        dframe$B <- insertdim (as.vector(unlist(temp)), c(2,3,1), dims)
    }

    #------------------------------------------------

    ## h2 or h3
    if (nmix > 1) {
        mixture <- paste('h',nmix,sep='')
        dframe[,mixture] <- insertdim(factor(1:nmix), 5, dims)
    }

    #--------------------------------------------------------------------------
    ## all autovars should have now been dealt with
    vars <- vars[!(vars %in% c(autovars, dframevars))]

    #--------------------------------------------------------------------------
    # add zcov, sessioncov, intervalcov

    if (!bygroup) findvars.MS (zcov, vars, 2)   ## CL only
    findvars.MS (sessioncov, vars, 1)
    findvars.MS (intervalcov, vars, 3)      ## session-specific list

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

    make.designmatrix <- function (formula, prefix, ...) {
     # combine formula and dframe to generate design matrix
        if (is.null(formula)) {
            list (model = NULL, index = rep(1,dframenrow))
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
    ## DOES THIS MESS UP FACTOR LEVELS?
    dframe[is.na(dframe)] <- 0

    # list with one component per real parameter
    # each of these is a list with components 'model' and 'index'

    designMatrices <- sapply (1:length(models), simplify=FALSE,
        function (x) make.designmatrix(models[[x]], names(models[x])))
    names(designMatrices) <- names(models)

    ## dim(indices) = c(R*n*J*nmix, npar)
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
    #       index to row of parameterTable for a given R,n,J,nmix
    # dim(parameterTable) = c(uniqueparcomb, npar)
    #       index to row of designMatrix for each real parameter
    #--------------------------------------------------------------------

    PIA <- array(indices2$index, dim = dims)
    parameterTable <- indices2$lookup
    colnames(parameterTable) <- parnames

    #--------------------------------------------------------------------

    list(designMatrices = designMatrices, parameterTable = parameterTable, PIA = PIA)
}
############################################################################################


