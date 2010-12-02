secr.cellprob <- function (beta, parindx, link, fixed, designD, design, design0 = NULL,
                           capthist, mask, detectfn = 0, CL = T, groups = NULL,
                           details, dig = 3)

# Return the cell probabilities for inhomogeneous Poisson spatial capture-recapture model

# Transformed parameter values (density, g0, sigma, z etc.) are passed in the vector 'beta'
# 'detectfn' is integer code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential
# 'CL' is logical for conditional (CL=T) vs full (CL=F) likelihood
#

{
    MS <- inherits(capthist,'list')
    if (MS) sessionlevels <- session(capthist)   ## was names(capthist) 2009 08 15
    else sessionlevels <- 1
    nsession <- length(sessionlevels)

    if (!detectfn %in% c(0:2,5))
        stop ("'detectfn' can only take values 0 (halfnormal), 1 (hazard-rate),",
              " 2 (exponential) or 5 (w-exponential) for now")

    #--------------------------------------------------------------------
    # Groups
    grp  <- group.factor (capthist, groups)
    ngrp <- max(1,length(group.levels(capthist, groups)))
    #--------------------------------------------------------------------
    # Detection parameters

    realparval  <- makerealparameters (design, beta, parindx, link, fixed)
    realparval0 <- makerealparameters (design0, beta, parindx, link, fixed)  # naive

    ## place holders for now 2009 09 27
    binomN <- 0
    cut <- 1.0
    spherical <- TRUE

    #--------------------------------------------------------------------
    # Density

    D.modelled <- !CL & is.null(fixed$D)
    if (!CL) {
         D <- getD (designD, beta, mask, parindx, link, fixed, MS, ngrp, nsession)
        if (sum(D)<=0) stop ('Invalid density')
    }

    #--------------------------------------------------------------------

if (nsession>1)
    stop ("not yet for multiple sessions")
if (CL)
    stop ("not yet for conditional likelihood")

    ####################################################
    ## start loop over sessions

    for (sessnum in 1:nsession) {

        ## in multi-session case must get session-specific data from lists
        if (MS) {
            session.capthist <- capthist[[sessnum]]
            session.traps    <- traps(capthist)[[sessnum]]
            session.mask     <- mask[[sessnum]]
            session.grp      <- grp[[sessnum]]
        }
        else {
            session.capthist <- capthist
            session.traps    <- traps(capthist)
            session.mask     <- mask
            session.grp      <- grp
        }

        if (nrow(session.capthist)==0) stop(paste('no data for session',sessnum))
        nc   <- nrow(session.capthist)
        s    <- ncol(session.capthist)
        m    <- nrow(session.mask)
        cell <- attr(session.mask,'area')
        session.mask <- as.matrix(session.mask[,1:2])

        session.detector <- detector(session.traps)
        if (session.detector %in% c('multi', 'single')) dettype <- 0
        else if (session.detector == 'proximity') dettype <- 1
        else if (session.detector == 'signal') dettype <- 2
        else if (session.detector == 'count') dettype <- 3
        else if (session.detector == 'areasearch') dettype <- 4
        else stop('Unrecognised detector type')

        if (dettype == 2) {    # signal strength
            session.capthist <- switch( details$tx,
                log = ifelse( session.capthist == 0, -1e20, log(session.capthist)),
                logit = ifelse( session.capthist == 0, -1e20, logit(session.capthist)),
                identity = ifelse( session.capthist == 0, -1e20, session.capthist)
            )
        }
        if (dettype == 4) {    # area search
            if (is.null(attr(session.traps,'spacing')))
                stop ("cell spacing required")
            searchcell <- attr(session.traps,'spacing')^2 * 0.0001
        }
        else searchcell <- 1    # for proximity

        k     <- nrow(session.traps)
        trps  <- unlist(session.traps, use.names=F)

        #--------------------------------
        ## differentiate so density & g do not both need to use sessions
        if (CL)
            density <- 0
        else
            density <- D[1:m,,min(dim(D)[3],sessnum)]

        sessg <- min (sessnum, design$R)

        #------------------------------------------
        # allow for scaling of detection parameters

##        Xrealparval  <- realparval
##        Xrealparval0 <- realparval0
##        sigmaindex <- 2
##        g0index <- 1
##        if (details$scalesigma) {   ## assuming previous check that scalesigma OK...
##            Xrealparval[,sigmaindex]  <- Xrealparval[,sigmaindex] / D[1,1,sessnum]^0.5
##            Xrealparval0[,sigmaindex] <- Xrealparval0[,sigmaindex] / D[1,1,sessnum]^0.5
##        }
##        if (details$scaleg0)    {   ## assuming previous check that scaleg0 OK...
##            Xrealparval[,g0index]  <- Xrealparval[,g0index] / Xrealparval[,sigmaindex]^2
##            Xrealparval0[,g0index] <- Xrealparval0[,g0index] / Xrealparval0[,sigmaindex]^2
##        }
##
        Dtemp <- ifelse (D.modelled, D[1,1,sessnum], NA)
        Xrealparval <- scaled.detection (realparval, details$scalesigma, details$scaleg0, Dtemp)
        Xrealparval0 <- scaled.detection (realparval0, details$scalesigma, details$scaleg0, Dtemp)

        #------------------------------------------

        ## For conditional likelihood, supply a value for each
        ##  animal, not just groups
        if (CL) tempPIA0 <- design0$PIA[sessg,1:nc,1:s,1:k]
        else tempPIA0 <- design0$PIA[sessg,1:ngrp,1:s,1:k]

cat ('dettype', dettype, '\n')
cat ('*cc', nrow(Xrealparval), '\n')

    temp <- .C('secrcellprob', PACKAGE = 'secr',
            as.integer(dettype),  # 0 = multicatch, 1 = proximity, 2 = signal, 3 = count, 4 = area search
            as.integer(session.capthist),
            as.integer(session.grp),
            as.integer(nc),
            as.integer(s),
            as.integer(k),
            as.integer(m),
            as.integer(ngrp),
            as.double(trps),
            as.double(session.mask),
            as.double(density),                                   # density at each mask point x ngrp cols
            as.double(Xrealparval),
            as.double(Xrealparval0),
            as.integer(nrow(Xrealparval)),                        # number of rows in lookup table
            as.integer(nrow(Xrealparval0)),                       # ditto, naive
            as.integer(design$PIA[sessg,1:nc,1:s,1:k]),           # index of nc,S,K to rows in Xrealparval
            as.integer(tempPIA0),                                 # index of ngrp,S,K to rows in Xrealparval0
            as.double(cell),
            as.integer(detectfn),
            as.integer(binomN),
            as.double(cut),
            as.integer(spherical),
            as.double(details$minprob),
            value=double(nc),
            resultcode=integer(1)
    )

        if (temp$resultcode != 0)
            stop ("error in external function 'secr.cellprob'")
    }
    ## end loop over sessions
    ####################################################

    pi.n <- temp$value
    names(pi.n) <- rownames(session.capthist)
    pi.n
}
############################################################################################




