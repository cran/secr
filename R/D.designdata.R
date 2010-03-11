############################################################################################
## package 'secr'
## D.design.MS.R
## Prepare density design matrix
## 
## force levels=sessionlevels for dframe$session 2010 02 25

############################################################################################
## NOTE does not standardize sessioncov, maskcov
############################################################################################

D.designdata <- function (mask, Dmodel, grps, sessionlevels, sessioncov = NULL) {

    ## mask may be a list of the same length as sessionlevels
    ## grps = group.levels(capthist,groups)

    stdfn <- function(x) if (is.numeric(x)) (x-mean(x))/sqrt(var(x)) else x

    #######################################
    findvars <- function (cov, vars, dimcov, std=TRUE) 
    ## function to add covariates to a design data frame
    {
        found <- names(cov) %in% vars
        if (is.data.frame(cov) & any(found)) {
            found <- names(cov)[found]
            values <- as.data.frame(cov[,found])
            names(values) <- found
            if (length(values)>0) {
                for (i in 1:ncol(values)) {    
                    vals <- values[,i]
                    if (!is.factor(vals) & std) vals <- stdfn(vals)
                    eval(parse(text=paste ('dframe$', found[i],
                          ' <<- insertdim (vals, dimcov, dims)', sep='') ))
                }
                vars <<- vars[!(vars %in% found)]
            }
        }
    }
    #######################################
    
    MS <- length(sessionlevels) > 1  ## !inherits(mask, 'mask')
    if (MS) nmask <- max(sapply(mask, nrow))
    else nmask <- nrow(mask)

    vars  <- all.vars(Dmodel)

    ##  Always conform to the global number of groups 2009 06 23
    ngrp  <- length(grps)

    R     <- length(sessionlevels)
    dims  <- c(nmask, ngrp, R)

    ## maybe return from here if no D model... 2009 03 04

    getcol <- function (msk, colnum) {
        pad1 (stdfn(msk[,colnum]), nmask)
    }

    if (MS) {
        x <- lapply(mask, getcol, 1)
        y <- lapply(mask, getcol, 2)   ## bug 2009 09 04: changed 1 -> 2
    }
    else {
        x <- getcol(mask,1)
        y <- getcol(mask,2)
    }

    dframe <- as.data.frame( list (
      x = rep(as.vector(unlist(x)), ngrp),
      y = rep(as.vector(unlist(y)), ngrp)
    ))

    if ('g' %in% vars) {   
        if (length(grps)<1) stop('No groups specified')
        dframe$g <- insertdim(factor(grps), 2, dims)
    } 

    ## force levels=sessionlevels 2010 02 25
    if ('session' %in% vars) {
       dframe$session <- insertdim(factor(sessionlevels, levels=sessionlevels), 3, dims)
    } 

    ## all autovars should have now been dealt with 
    vars <- vars[!vars %in% c('g','x','y','session')]

    if (!is.null(sessioncov)) {
       findvars (sessioncov, vars, 3, std = FALSE)   
    }

    if (!is.null(covariates(mask))) {
       findvars (covariates(mask), vars, 1, std = FALSE)   
    }

    if (length(vars)>0) stop (paste(paste(vars,collapse=','),'not found'))

    attr(dframe, 'dimD') <- c(nmask, ngrp, R)
    dframe
}
############################################################################################

