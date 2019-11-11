###############################################################################
## package 'secr'
## fastsecrloglik.R
## likelihood evaluation functions
## last changed 2019-09-01
###############################################################################

#--------------------------------------------------------------------------------
allhistfast <- function (realparval, gkhk, pi.density, PIA, 
                         nk2ch, usge, pmixn, maskusage,
                         grain, binomN) {
    nc <- dim(PIA)[2]
    nmix <- dim(PIA)[5]
    m <- length(pi.density)
    sump <- numeric(nc)
    for (x in 1:nmix) {
        temp <- fasthistoriescpp(
            as.integer(m),
            as.integer(nc),
            as.integer(nrow(realparval)),
            as.integer(grain),
            as.integer(binomN),
            matrix(nk2ch[,,1], nrow=nc),
            matrix(nk2ch[,,2], nrow=nc),
            as.double (gkhk$gk),  ## precomputed probability 
            as.double (gkhk$hk),  ## precomputed hazard 
            as.double (pi.density),
            as.integer(PIA[1,,,,x]),
            as.integer(usge),
            as.matrix (maskusage))
        sump <- sump + pmixn[x,] * temp
    }
    sump
}
#--------------------------------------------------------------------------------

integralprw1fast <- function (realparval0, gkhk, pi.density, PIA0, 
                              nk2ch0, usge, pmixn, grain, binomN) {
    nc <- dim(PIA0)[2]
    nr <- 1    ## nrow(nk2ch0)
    nmix <- dim(PIA0)[5]
    m <- length(pi.density)
    sump <- numeric(nc)
    for (x in 1:nmix) {
        temp <- fasthistoriescpp(
            as.integer(m),
            as.integer(nr),    ## 1 
            as.integer(nrow(realparval0)),
            as.integer(0),                     ## force grain = 0 as single capture history
            as.integer(binomN),
            matrix(nk2ch0[,,1], nrow = nr),
            matrix(nk2ch0[,,2], nrow = nr),
            as.double (gkhk$gk),        ## precomputed probability 
            as.double (gkhk$hk),        ## precomputed hazard
            as.double (pi.density),
            as.integer(PIA0[1,,,,x]),
            as.integer(usge),
            as.matrix (matrix(TRUE, nrow = nc, ncol = m)))
        sump <- sump + pmixn[x,1:nc] * (1 - temp)
    }
    sump
}

#######################################################################################
fastsecrloglikfn <- function (
    beta, 
    parindx, 
    link, 
    fixedpar, 
    designD, 
    designNE, 
    design, 
    design0, 
    CL, 
    detectfn,
    learnedresponse,
    sessionlevels,
    data,
    details,
    dig = 3, betaw = 10, neglik = TRUE)

# Return the negative log likelihood for spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function
#    0 = halfnormal, 1 = hazard, 2 = exponential etc.
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
# details$trace=T sends a one-line report to the screen

{
    #--------------------------------------------------------------------------------
    sessionLL <- function (data, sessnum, like = 0) {
        ## log likelihood for one session
        ## in multi-session case must get session-specific data from lists
        ###################################################################
        #---------------------------------------------------
      PIA <- design$PIA[sessnum, 1:data$nc, 1:data$s, 1:data$K, ,drop=FALSE]
        ## unmodelled beta parameters, if needed
        miscparm <- getmiscparm(details$miscparm, detectfn, beta, parindx, details$cutval)
        density <- getmaskpar(!CL, D, data$m, sessnum, details$unmash, 
                              attr(data$capthist, 'n.mash'))
        if (CL) {
            pi.density <- rep(1/data$m, data$m)  
        }
        else {
            pi.density <- density / sum(density)
        }
        #---------------------------------------------------
        ## allow for scaling of detection
        Dtemp <- if (D.modelled) mean(D[,1,sessnum]) else NA
        Xrealparval <- reparameterize (realparval, detectfn, details,
                                       data$mask, data$traps, Dtemp, s)
        ## check valid parameter values
        if (!all(is.finite(Xrealparval))) {
            cat ('beta vector :', beta, '\n')
            warning ("extreme 'beta' in 'fastsecrloglikfn' ",
                     "(try smaller stepmax in nlm Newton-Raphson?)")
            return (1e10)
        }
        
        ## DOES NOT ALLOW FOR GROUP VARIATION IN DENSITY
        ## more thoughts 2015-05-05
        ## could generalize by
        ## -- making Dtemp a vector of length equal rows in realparval
        ## -- matching either
        ##      first group (as before)
        ##      sum of all groups
        ##      own group [PROBLEM: locating group of each realparval row]
        ## in all cases density is the mean over mask points
        
        ## CHECK use of Dtemp in regionN.R, sim.secr.R
        ## PERHAPS for consistency make a function to construct Dtemp vector
        ## given mask, model, group matching rule (first, sum, own)
        
        ##--------------------------------------------------------------

        #####################################################################
        pmixn <- getpmix (data$knownclass, PIA, Xrealparval)  ## membership prob by animal
        if (is.function(details$userdist)) {
          noneuc <- getmaskpar(!is.null(NE), NE, data$m, sessnum, FALSE, NULL)
          distmat2 <- getuserdist(data$traps, data$mask, details$userdist, sessnum, 
                                  noneuc[,1], density[,1], miscparm)
        }
        else
            distmat2 <- data$distmat2
        
        ## precompute gk, hk for point detectors
        if (data$dettype[1] %in% c(0,1,2,5,8)) {
            gkhk <- makegkParallelcpp (as.integer(detectfn),
                                       as.integer(details$grain),
                                       as.matrix(Xrealparval),
                                       as.matrix(distmat2),
                                       as.double(miscparm))
            if (details$anycapped) {   ## capped adjustment
              gkhk <- cappedgkhkcpp (
                as.integer(nrow(Xrealparval)),
                as.integer(nrow(data$traps)),
                as.double(attr(data$mask, "area")),
                as.double(density[,1]),
                as.double(gkhk$gk), as.double(gkhk$hk))  
            }
        }
        prw <- allhistfast (Xrealparval, gkhk, pi.density, PIA, 
                            data$CH, data$usge, pmixn, data$maskusage, 
                            details$grain, details$binomN)
        pdot <- integralprw1fast (Xrealparval, gkhk, pi.density, PIA, 
                                  data$CH0, data$usge, pmixn, details$grain, details$binomN)
        
        comp <- matrix(0, nrow = 5, ncol = 1)
        comp[1,1] <- if (any(is.na(prw) || prw<=0)) NA else sum(log(prw))
        comp[2,1] <- if (any(is.na(pdot) || pdot<=0)) NA else -sum(log(pdot))
        if (!CL) {
            N <- sum(density[,1]) * getcellsize(data$mask)
            meanpdot <- data$nc / sum(1/pdot)
            comp[3,1] <- switch (data$n.distrib+1,
                                 dpois(data$nc, N * meanpdot, log = TRUE),
                                 lnbinomial (data$nc, N, meanpdot),
                                 NA)
        }
        ## adjustment for mixture probabilities when class known
        known <- sum(data$knownclass>1)
        if (details$nmix>1 && known>0) {
            nm <- tabulate(data$knownclass, nbins = max(data$knownclass))
            pmix <- attr(pmixn, 'pmix')
            for (x in 1:details$nmix) {
                # need group-specific pmix
                comp[4,1] <- comp[4,1] + nm[x+1] * log(pmix[x]) 
            }
        }
        
        if (details$debug>1) {
            comp <- apply(comp,1,sum)
            cat(comp[1], comp[2], comp[3], comp[4], comp[5], comp[6], '\n')
        }
        sum(comp) + data$logmult
        
    } ## end sessionLL
    
    ###############################################################################################
    ## Main line of fastsecrloglikfn

    nsession <- length(sessionlevels)
    #--------------------------------------------------------------------
    # Fixed beta
    beta <- fullbeta(beta, details$fixedbeta)
    
    #--------------------------------------------------------------------
    # Detection parameters
    detparindx <- parindx[!(names(parindx) %in% c('D', 'noneuc'))]
    detlink <- link[!(names(link) %in% c('D', 'noneuc'))]
    realparval  <- makerealparameters (design, beta, detparindx, detlink, fixedpar)
    
    #--------------------------------------------------------------------
    # Density
    D.modelled <- !CL & is.null(fixedpar$D)
    if (!CL ) {
        sessmask <- lapply(data, '[[', 'mask')
        grplevels <- unique(unlist(lapply(data, function(x) levels(x$grp))))
        D <- getD (designD, beta, sessmask, parindx, link, fixedpar,
                   grplevels, sessionlevels, parameter = 'D')
        if (!is.na(sumD <- sum(D)))
            if (sumD <= 0)
                warning ("invalid density <= 0")
    }
    
    #--------------------------------------------------------------------
    # Non-Euclidean distance parameter
    sessmask <- lapply(data, '[[', 'mask')
    NE <- getD (designNE, beta, sessmask, parindx, link, fixedpar,
                levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
    #--------------------------------------------------------------------
    # typical likelihood evaluation
    loglik <- sum(mapply (sessionLL, data, 1:nsession))
    .localstuff$iter <- .localstuff$iter + 1   ## moved outside loop 2011-09-28
    if (details$trace) {
        fixedbeta <- data[[1]]$details$fixedbeta
        if (!is.null(fixedbeta))
            beta <- beta[is.na(fixedbeta)]
        cat(format(.localstuff$iter, width=4),
            formatC(round(loglik,dig), format='f', digits=dig, width=10),
            formatC(beta, format='f', digits=dig+1, width=betaw),
            '\n')
        flush.console()
    }
    loglik <- ifelse(is.finite(loglik), loglik, -1e10)
    ifelse (neglik, -loglik, loglik)
}  ## end of fastsecrloglikfn
############################################################################################

