##############################################################################
## package 'secr'
## fxi.R
## 2019-08-17 fxi.secr uses C++ call
###############################################################################

fxi2SPDF <- function (x, ID, levels) {
  if (missing(ID))
    ID <- 1:length(x)
  if (missing(levels))
    levels <- names(x[[1]])[names(x[[1]]) != 'mode']
  getxy1 <- function(one)
    lapply(one[levels], function (xx) Polygon(cbind(xx$x, xx$y)))
  oneanimal <- function (x,id)
    Polygons(x, id)
  xy <- lapply(x[ID], getxy1)
  modes <- t(sapply(x[ID], '[[', 'mode') )
  modes <- matrix(unlist(modes), ncol = 2)
  listSrs <- mapply(oneanimal, xy, ID)
  SpP <- SpatialPolygons(listSrs)
  df <- data.frame(modex = modes[,1], modey = modes[,2], row.names = ID)
  SpatialPolygonsDataFrame(SpP, df)
}

##    writeSpatialShape(SPDF, ...)

fxi.contour <- function (object, i = 1, sessnum = 1, border = 100, nx = 64,
                         levels = NULL, p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE, fitmode =
                           FALSE, plotmode = FALSE, fill = NULL, SPDF = FALSE,  ncores = NULL, ...) {
  if (inherits(object$mask, 'linearmask'))
    stop("contouring fxi is not appropriate for linear habitat")
  if (ms(object)) {
    session.traps <- traps(object$capthist)[[sessnum]]
  }
  else {
    session.traps <- traps(object$capthist)
  }
  tempmask <- make.mask (session.traps, border, nx = nx, type = 'traprect')
  xlevels <- unique(tempmask$x)
  ylevels <- unique(tempmask$y)
  
  fxi <- function (ni) {
    z <- allz[[ni]]
    if (is.null(levels)) {
      temp <- sort(z, decreasing = T)
      cump <- cumsum(temp) / sum(temp)
      levels <- suppressWarnings(approx (x = cump, y = temp, xout = p)$y)
      labels <- p
    }
    else
      labels <- levels
    
    templines <- contourLines(xlevels, ylevels, matrix(z, nrow = nx), levels = levels)
    ## extra effort to apply correct labels
    getlevels <- function(clines.i) sapply(clines.i, function(q) q$level)
    label.levels <- function (island) {
      which.levels <- match (getlevels(island), levels)
      names(island) <- labels[which.levels]
      island
    }
    templines <- label.levels(templines)
    
    wh <- which.max(unlist(lapply(templines, function(y) y$level)))
    if (length(templines) > 0) {  
      cc <- templines[[wh]]
      cc <- data.frame(cc[c('x','y')])
      templines$mode <- data.frame(x=mean(cc$x), y=mean(cc$y))
      if (fitmode) {
        templines$mode <- fxi.mode(object, sessnum = sessnum,
                                   start = templines$mode, i = ni)
      }
      if (plt) {
        labels <- labels[!is.na(levels)]
        levels <- levels[!is.na(levels)]
        
        contour (xlevels, ylevels, matrix(z, nrow = nx), add = add,
                 levels = levels, labels = labels, ...)
        
        ## optional fillin
        if (!is.null(fill)) {
          z[z < (0.999 * min(levels))] <- NA
          levels <- rev(c(1,levels,0))
          .filled.contour(xlevels, ylevels,  matrix(z, nrow = nx), levels= levels,
                          col = fill)
        }
        
        if (plotmode) {
          points(templines$mode, col = 'red', pch = 16)
        }
        
      }
    }
    templines
  }
  
  allz <- fxi.secr(object, i=i, sessnum = sessnum, X = tempmask, ncores = ncores)
  if (!is.list(allz))
    allz <- list(allz)
  temp <- lapply(1:length(allz), fxi)
  
  if (SPDF)
    temp <- fxi2SPDF(temp)
  
  if (plt)
    invisible(temp)
  else
    temp
}

###############################################################################

fxi.mode <- function (object, i = 1, sessnum = 1, start = NULL, ncores = NULL, ...) {
  if (length(i)>1) stop ("fxi.mode takes single i")
  if (ms(object))
    session.capthist <- object$capthist[[sessnum]]
  else
    session.capthist <- object$capthist
  start <- unlist(start)
  if (is.null(start)) {
    session.traps <- traps(session.capthist)
    animal <- animalID(session.capthist, names=F) == i
    trp <- trap(session.capthist)[animal]
    start <- sapply(traps(session.capthist)[trp,],mean)
  }
  if (is.character(i))
    i <- match(i, row.names(session.capthist))
  if (is.na(i) | (i<1) | (i>nrow(session.capthist)))
    stop ("invalid i in fxi.secr")
  fn <- function(xy,i) {
    -fxi.secr(object, i=i, sessnum = sessnum, X = matrix(xy, ncol=2), ncores = ncores)[[1]]
  }
  temp <- nlm(f = fn, p = start, i = i, typsize = start, ...)$estimate
  data.frame(x=temp[1], y=temp[2])
}

###############################################################################

## mask if specified should be for a single session
## ... passes newdata df to predict.secr

fx.total <- function (object, sessnum = 1, mask = NULL, ncores = NULL, ...)
{
  if (ms(object)) {
      n <- nrow(object$capthist[[sessnum]])
      if (is.null(mask)) mask <- object$mask[[sessnum]]
      detectpar <- detectpar(object, ...)[[sessnum]]
      CH <- object$capthist[[sessnum]]
  }
  else {
      n <- nrow(object$capthist)
      if (is.null(mask)) mask <- object$mask
      detectpar <- detectpar(object, ...)
      CH <- object$capthist
  }
  fxi <- fxi.secr(object, sessnum = sessnum, X = mask, ncores = ncores)
  fx <- do.call(cbind, fxi)
  fxt <- apply(fx, 1, sum)
  fxt <- fxt/getcellsize(mask)
  D <- predictDsurface(object, mask = mask)
  D <- covariates(D)$D.0
  pd <- pdot(X = mask, traps = traps(CH), detectfn = object$detectfn,
             detectpar = detectpar, noccasions = ncol(CH), ncores = ncores)
  nct <- D * (1 - pd)
  covariates(mask) <- data.frame(D.fx = fxt, D.nc = nct, D.sum = fxt + nct)
  class(mask) <- c("Dsurface", class(mask))
  mask
}

###############################################################################

allhistfxi <- function (m, realparval, haztemp, gkhk, pi.density, PIA, usge,
                        CH, binomN, grp, pmixn, grain) {
    nc <- nrow(CH)
    nmix <- nrow(pmixn)
    sump <- matrix(0, nrow = nc, ncol = m)
    for (x in 1:nmix) {
        hx <- if (any(binomN==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
        hi <- if (any(binomN==-2)) haztemp$hindex else -1                   ## index to hx
   
        temp <- simplehistoriesfxicpp(
            as.integer(x-1),
            as.integer(m),
            as.integer(nc),
            as.integer(nrow(realparval)),
            as.integer(grain),
            as.integer(binomN),
            as.integer(CH),   
            as.integer(grp)-1L,
            as.double (gkhk$gk),     ## precomputed probability 
            as.double (gkhk$hk),     ## precomputed hazard
            as.matrix (pi.density),
            as.integer(PIA),
            as.matrix(usge),
            as.matrix (hx),                
            as.matrix (hi))
        sump <- sump + sweep(temp, MARGIN=1, STATS = pmixn[x,], FUN = "*")
    }
    sump
}

allhistpolygonfxi <- function (detectfn, realparval, haztemp, hk, H, pi.density, PIA, 
                               CH, xy, binomNcode, grp, usge, mask, pmixn, maskusage, 
                               grain, minprob) {
    nc <- nrow(CH)
    nmix <- nrow(pmixn)
    m <- length(pi.density)
    s <- ncol(usge)
    sump <- matrix(0, nrow = nc, ncol = m)
    for (x in 1:nmix) {
        hx <- if (any(binomNcode==-2)) matrix(haztemp$h[x,,], nrow = m) else -1 ## lookup sum_k (hazard)
        hi <- if (any(binomNcode==-2)) haztemp$hindex else -1                   ## index to hx
        temp <- polygonfxicpp(
            as.integer(nc),
            as.integer(detectfn[1]),
            as.integer(grain),
            as.double(minprob),          
            as.integer(binomNcode),   # vector length s
            as.integer(CH),   
            as.matrix(xy$xy),
            as.vector(xy$start),
            as.integer(as.numeric(grp))-1L,
            as.double(hk),
            as.double(H),
            as.matrix(realparval),
            matrix(1,nrow=s, ncol=nmix),  ## pID?
            as.matrix(mask),
            as.matrix (pi.density),
            as.integer(PIA[1,,,,x]),
            as.matrix(usge),
            as.matrix (hx),                
            as.matrix (hi),      
            as.matrix(maskusage)
        )
        sump <- sump + sweep(temp, MARGIN=1, STATS = pmixn[x,], FUN = "*")
    }
    sump
}

fxi.secr <- function (object, i = NULL, sessnum = 1, X = NULL, ncores = NULL) {
    
    ## temporary fix for lack of fastproximity code
    object$details$fastproximity <- FALSE   ## 2020-08-30
    
    ## data for a single session
    data <- prepareSessionData(object$capthist, object$mask, object$details$maskusage, 
                             object$design, object$design0, object$detectfn, object$groups, 
                             object$fixed, object$hcov, object$details)
    
  sessionlevels <- session(object$capthist)
  beta <- coef(object)$beta
  beta <- fullbeta(beta, object$details$fixedbeta)
  detparindx <- object$parindx[!(names(object$parindx) %in% c('D', 'noneuc'))]
  detlink <- object$link[!(names(object$link) %in% c('D', 'noneuc'))]
  realparval  <- makerealparameters (object$design, beta, detparindx,
                                     detlink, object$fixed)
  data <- data[[sessnum]]
  reusemask <- is.null(X)
  if (reusemask) {
    X <- data$mask
  }
  else {
    X <- matrix(unlist(X), ncol = 2)
  }
  #----------------------------------------
  # restrict to selected individuals
  xy <- data$xy 
  if (is.null(i)) {
    ok <- 1:nrow(data$CH)
  }
  else {
    ok <- i
    if (!is.null(xy)) {
        xy <- getxy(data$dettype, selectCHsession(subset(object$capthist, ok), sessnum))
    }
  }
  if (length(dim(data$CH)) == 2) {
    CH <- data$CH[ok,,drop=FALSE]
  }
  else {
    CH <- data$CH[ok,,,drop=FALSE]
  }
  grp <- data$grp[ok]

  setNumThreads(ncores)
  
  #----------------------------------------
  # Density
  if (is.null(object$model$D))
    D.modelled <- FALSE
  else {
    if (!is.null(object$fixed$D))
      D.modelled <- FALSE
    else
      D.modelled <- (object$model$D != ~1)
  }
  if (D.modelled) {
    predD <- predictDsurface (object)
    if (ms(object))
      predD <- predD[[sessnum]]
    D <- covariates(predD)$D.0  ## does not apply if groups
    pimask <- D / sum(D)   ## vector of probability mass for each mask cell
  }
  else {
    mm <- nrow(data$mask)
    pimask <- rep(1, mm)  ## could be 1/mm, but as we normalise anyway...
  }
  ## fetch predicted density at each new point X
  ## covariates(session.mask) <- data.frame(pi = pimask)
  if (!is.null(covariates(data$mask)))
    covariates(data$mask) <- cbind(data.frame(pi = pimask), covariates(data$mask))
  else
    covariates(data$mask) <- data.frame(pi = pimask)
  ## does this work for linearmask?
  tmpmask <- suppressWarnings(addCovariates(X, data$mask, strict = TRUE))
  piX <- covariates(tmpmask)$pi
  piX[is.na(piX)] <- 0
  #----------------------------------------
  
  ## TO BE FIXED
   
  # NE <- getD (object$designNE, beta, object$mask, object$parindx, object$link, object$fixed,
  #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
  # NEX <- getD (object$designNE, beta, X, object$parindx, object$link, object$fixed,
  #             levels(data$grp[[1]]), sessionlevels, parameter = 'noneuc')
  # 
  NE <- NULL
  
  #---------------------------------------------------
  ## allow for scaling of detection
  Dtemp <- if (D.modelled) mean(D) else NA
  Xrealparval <- reparameterize (realparval, object$detectfn, object$details,
                                 data$mask, data$traps, Dtemp, s)
  PIA <- object$design$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
  PIA0 <- object$design0$PIA[sessnum, ok, 1:data$s, 1:data$K, ,drop=FALSE]
  pmix <- getpmix (data$knownclass[ok], PIA, Xrealparval)  ## membership prob by animal

  ## unmodelled beta parameters, if needed
  miscparm <- getmiscparm(object$details$miscparm, object$detectfn, object$beta, 
                          object$parindx, object$details$cutval)

  gkhk <- makegk (data$dettype, object$detectfn, data$traps, data$mask, object$details, sessnum, 
                  NE, D, miscparm, Xrealparval)
  haztemp <- gethazard (data$m, data$binomNcode, nrow(realparval), gkhk$hk, PIA, data$usge)
  
  ## 2020-01-26 conditional on point vs polygon detectors
  if (data$dettype[1] %in% c(0,1,2,5,8,13)) {
      
      prmat <- allhistfxi (data$m, Xrealparval, haztemp, gkhk, pimask, PIA, data$usge,
                           CH, data$binomNcode, grp, pmix, object$details$grain)
  }
  else {
      # warning ("fxi.secr experimental for polygon detector types")
      prmat <- allhistpolygonfxi (object$detectfn, Xrealparval, haztemp, gkhk$hk, gkhk$H, pimask, PIA, 
          CH, xy, data$binomNcode, grp, data$usge, data$mask,
          pmix, data$maskusage, object$details$grain, object$details$minprob)
  }
  
  pisum <- apply(prmat,1,sum)
  
  if (reusemask) {
    out <- sweep(prmat, MARGIN=1, STATS=pisum, FUN="/")
  }
  else {
    nX <- nrow(X)
    gkhkX <- makegk (data$dettype, object$detectfn, data$traps, X, object$details, sessnum, NE, D, miscparm, Xrealparval)
    haztempX <- gethazard (nX, data$binomNcode, nrow(realparval), gkhkX$hk, PIA, data$usge)
    
    if (data$dettype[1] %in% c(0,1,2,5,8,13)) {
        ## point detectors
        prmatX <- allhistfxi (nX, Xrealparval, haztempX, gkhkX, piX, PIA, data$usge,
                          CH, data$binomNcode, grp, pmix, object$details$grain)
    }
    else {
        ## polygon-like detectors
        prmatX <- allhistpolygonfxi (object$detectfn, Xrealparval, haztempX, gkhkX$hk, gkhkX$H, piX, PIA, 
                             CH, xy, data$binomNcode, grp, data$usge, X,
                             pmix, data$maskusage, object$details$grain, object$details$minprob)
    }
    out <- sweep(prmatX, MARGIN=1, STATS=pisum, FUN="/")
  }
  out <- as.list(as.data.frame(t(out)))
  names(out) <- row.names(CH)
  out
}