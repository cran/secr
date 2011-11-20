############################################################################################
## package 'secr'
## Dsurface.R
## 2011-10-21, modified through 2011-11-10
############################################################################################

predictD <- function (object, regionmask, group, session,
                      se.D = FALSE, cl.D = FALSE, alpha = 0.05) {

    ## For one session and group at a time
    ## not exported; used also by region.N()

    sessionlevels <- session(object$capthist)
    grouplevels <- group.levels(object$capthist,object$groups)
    if (is.null(session))
        session <- sessionlevels[1]
    if (is.null(group))
        group <- grouplevels[1]
    if (is.null(regionmask))
        regionmask <- object$mask
    if (ms(regionmask))
        regionmask <- regionmask[[session]]

    group <- match (group, grouplevels)   ## numeric

    ## must use mean and SD for scaling from original mask(s)
    ## this is a list for ms masks
    if (ms(object))
        meanSD <- attr(object$mask[[session]], 'meanSD')
    else
        meanSD <-attr(object$mask, 'meanSD')

    ## 2011-11-08 allow for 'mashing'
    if (ms(object))  {
        n.mash <- attr (object$capthist[[session]], 'n.mash')
        n.clust <- length(n.mash)
    }
    else {
        n.mash <- attr (object$capthist, 'n.mash')
        n.clust <- length(n.mash)
    }
    if (is.null(n.mash))
        n.clust <- 1

    ## no density model (conditional likelihood fit)
    if (object$CL == TRUE) {    ## implies is.null(object$model$D)
        temp <- derived(object, se.D = (se.D | cl.D)) ## inefficient as repeats for each sess
        if (!is.data.frame(temp))
            temp <- temp[[session]]
        D <- temp['D', 'estimate'] / n.clust
        if (se.D) {
            attr(D, 'SE') <- temp['D', 'SE'] / n.clust
        }
        if (cl.D) {
            z <- abs(qnorm(1-alpha/2))
            attr(D, 'lcl') <- temp['D', 'lcl'] / n.clust
            attr(D, 'ucl') <- temp['D', 'lcl'] / n.clust
        }
        return (D)
    }

    ## user-defined density model
    else if (userD(object)) {
        designD <- object$details$userDfn
        if (!is.function(designD))
            stop ("details$userDfn must be a function")
        if (se.D | cl.D)
            warning ("SE and confidence intervals may be unavailable with userDfn")
        if (n.clust>1)
            warning ("no adjustment for mashing when using userDfn")
        ## getD is in functions.R
        ## does not use link$D when calling user function userDfn
        D <- getD (designD, object$fit$par, regionmask, object$parindx,
                   object$link$D, object$fixed, grouplevels, sessionlevels)
        return(D[,group,session])
    }
    ## linear density model on link scale
    else {
        newdata <- D.designdata (regionmask, object$model$D,
             grouplevels, sessionlevels, sessioncov = object$sessioncov,
             meanSD = meanSD)
        dimD <- attr(newdata, "dimD")
        ## if newdata has more than one group or session...
        if (prod(dimD[2:3]) > 1) {
            ## select a single group
            groupOK <- (group == grouplevels) | (dimD[2]==1)
            groupOK <- rep(rep(groupOK, each = dimD[1]), dimD[3])
            ## select a single session
            sessionOK <- if (is.character(session))
                session == sessionlevels
            else
                session == 1:dimD[3]
            sessionOK <- rep(sessionOK, each = prod(dimD[1:2]))
            newdata <- newdata[sessionOK & groupOK,]
        }
        class(newdata) <- c('mask', 'data.frame')
        attr (newdata, 'area') <- attr(regionmask, 'area')
        indx <- object$parindx$D
        betaD <- object$fit$par[indx]
        if (object$model$D == ~1) {
            D <- untransform(betaD, object$link$D)
            D <- max(D,0) / n.clust
            return ( rep(D, nrow(regionmask) ))
        }
        else {
            vars <- all.vars(object$model$D)
            if (any(!(vars %in% names(newdata))))
                stop ("one or more model covariates not found")
            newdata <- as.data.frame(newdata)
            mat <- model.matrix(object$model$D, data=newdata)
            lpred <- mat %*% betaD
            temp <- untransform(lpred, object$link$D)
            temp <- pmax(temp, 0) / n.clust

            if (se.D | cl.D) {
                vcv <- vcov(object) [indx,indx]
                selpred <- sapply(1:nrow(mat), function(i)
                    mat[i,, drop=F] %*% vcv %*% t(mat[i,, drop=F]))^0.5
                if (se.D) {
                    attr(temp, 'SE') <- se.untransform (lpred, selpred, object$link$D) / n.clust
                    attr(temp, 'SE')[temp<=0] <- NA
                }
                if (cl.D) {
                    z <- abs(qnorm(1-alpha/2))
                    attr(temp, 'lcl') <- untransform (lpred - z * selpred, object$link$D) / n.clust
                    attr(temp, 'ucl') <- untransform (lpred + z * selpred, object$link$D) / n.clust
                    attr(temp, 'lcl')[temp<=0] <- NA
                    attr(temp, 'ucl')[temp<=0] <- NA
                }
            }
            return (temp)
        }
    }
}
############################################################################################

rectangularMask <- function (mask) {
    if (ms(mask)) {
        lapply(mask, rectangularMask)
    }
    else {
        temp <- expand.grid (x=sort(unique(mask$x)), y=sort(unique(mask$y)))
        temp <- read.mask(data = temp, spacing = spacing(mask))
        OK <- match(interaction(temp), interaction(mask))
        if (!is.null(covariates(mask))) {
            covariates(temp) <- covariates(mask)[OK,, drop = FALSE]
            rownames(covariates(temp)) <- 1:nrow(temp)
        }
        class(temp) <- class(mask)
        attr(temp, 'polygon') <-     attr(mask, 'polygon')
        attr(temp, 'poly.habitat') <-     attr(mask, 'poly.habitat')
        attr(temp, 'type') <- 'user'
        attr(temp, 'meanSD') <- attr(mask, 'meanSD')
        attr(temp, 'area') <-  attr(mask, 'area')
        attr(temp, 'boundingbox') <- attr(mask, 'boundingbox')
        attr(temp, 'OK') <- !is.na(OK)
        temp
    }
}
############################################################################################

predictDsurface <- function (object, mask = NULL, se.D = FALSE, cl.D = FALSE, alpha = 0.05) {
    sessionlevels <- session(object$capthist)
    grouplevels <- group.levels(object$capthist, object$groups)
    if (is.null(mask))
        mask <- object$mask
    densitylist <- vector('list')
    for (session in sessionlevels) {
        if (ms(mask))
            sessmask <- mask[[session]]
        else
            sessmask <- mask
        D <- data.frame (seq = 1:nrow(sessmask))
        for (group in grouplevels) {
            predicted <- predictD(object, sessmask, group, session, se.D, cl.D, alpha)
            D[,paste('D',group,sep='.')] <- as.numeric(predicted)
            if (se.D) {
                D[,paste('SE',group,sep='.')] <- as.numeric(attr(predicted, 'SE'))
            }
            if (cl.D) {
                D[,paste('lcl',group,sep='.')] <- as.numeric(attr(predicted, 'lcl'))
                D[,paste('ucl',group,sep='.')] <- as.numeric(attr(predicted, 'ucl'))
            }
        }
        if (is.null(covariates(sessmask)))
            tempcov <- D
        else
            tempcov <-  cbind (covariates(sessmask), D)
        tempcov$seq <- NULL ## drop dummy

        ## impose null observations
        if (!is.null(attr(sessmask, 'OK')))
             tempcov[!attr(sessmask,'OK'),] <- NA
        if (ms(mask))
            covariates(mask[[session]]) <- tempcov
        else
            covariates(mask) <- tempcov
    }
    if (ms(mask)) {
        ## drop 'data.frame' because fools ms(), but worrying...
        class (mask) <- c('list', 'Dsurface', 'mask')
        for (i in 1:length(mask))
            class(mask[[i]]) <- c('Dsurface','mask', 'data.frame')
    }
    else
        class (mask) <- c('Dsurface', 'mask', 'data.frame')
    mask
}
############################################################################################

plot.Dsurface <- function (x, covariate = 'D', group = NULL, plottype = 'shaded',
     scale = 1, ...) {
    if (ms(x)) {
        breaklist <- lapply(x, plot, covariate, group, plottype, ...)
        invisible(breaklist)
    }
    else {
        if (is.null(group))
            group <- 0
        if (length(covariate)>1)
            stop ("whoa... just one at a time")
        if (covariate %in% c('D','SE','lcl','ucl')) {
            covariate <- paste(covariate, group, sep='.')
        }
        if (!(covariate %in% names(covariates(x))))
            stop ("covariate ", covariate, " not found")
        covariates(x)[,covariate] <- covariates(x)[,covariate] * scale
        if (plottype %in% c('contour','persp')) {
            xval <- sort(unique(x$x))
            yval <- sort(unique(x$y))
            if (nrow(x) != length(xval)*length(yval)) {
                x <- rectangularMask(x)
                if(nrow(x) != length(xval)*length(yval))
                    stop ("failed to convert irregular mask to rectangle")
            }
            zmat <- matrix(covariates(x)[,covariate], nrow = length(xval))
            if (plottype == 'contour')
                contour(x=xval, y=yval, z=zmat, ...)
            else
                persp(x=xval, y=yval, z=zmat, ...)
        }
        else {
            class(x) <- c('mask','data.frame')
            covlevels <- plot(x, covariate = covariate, dots = (plottype == 'dots'), ...)
            if (!is.null(covlevels)) invisible(covlevels)
        }
    }
}
############################################################################################

Dsurface.as.data.frame <- function (x, scale = 1) {
    covnames <- names(covariates(x))
    OK <- (substring(covnames,1,2) == 'D.') |
        (substring(covnames,1,3) == 'SE.') |
        (substring(covnames,1,4) %in% c('lcl.', 'ucl.'))
    covnames <- covnames[OK]
    densities <- covariates(x)[,covnames] * scale
    df <- cbind(x, densities)
    names(df) <- c('x','y',covnames)
    df
}
############################################################################################

print.Dsurface <- function (x, scale = 1, ...) {
    if (ms(x)) {
        out <- vector('list')
        for (session in names(x)) {
            cat ('Session ', session, '\n')
            print(x[[session]], scale, ...)
            out[[session]] <- x[[session]]
        }
        names(out) <- names(x)
        out
    }
    else {
        df <- Dsurface.as.data.frame(x, scale)
        print(df, ...)
    }
    invisible(df)
}
############################################################################################

summary.Dsurface <- function (object, scale = 1, ...) {
    if (ms(object)) {
        temp <- lapply(object, summary, scale, ...)
        class(temp) <- c('summary.Dsurface', 'list')
      temp
    }
    else {
        covnames <- names(covariates(object))
        covnames <- covnames[substring(covnames,1,8) == 'D.']
        densities <- covariates(object)[,covnames,drop=FALSE]
        densities <- densities * scale
        if (dim(densities)[2] > 1)
            densities$Total <- apply(densities,1,sum)

         class(object) <- c('mask','data.frame')
         covariates(object) <- NULL
         tempmasksummary <- summary(object,...)

        list (
            mask = tempmasksummary,
            density = apply (densities, 2, summary))
    }
}
############################################################################################

spotHeight <- function (object, prefix = NULL, dec = 2, point = FALSE, text = TRUE,
                        sep = ', ', session = 1, scale = 1, ...) {
    if ((!inherits(object, 'mask')) | (is.null(covariates(object))))
        stop ("requires plotted Dsurface or mask with covariates")
    ## Esc or click outside to break
    decxy <- 0   ## decimal places for x-y coordinates
    if (is.null(prefix)) {
        if (inherits(object,'Dsurface'))
            prefix <- 'D.'
        else
            prefix <- ''
    }

    # can only deal with one session
    if (ms(object))
        object <- object[[session]]

    covnames <- names(covariates(object))
    if (all(prefix == ""))
        OK <- rep(TRUE, times=length(covnames))
    else
        OK <- pmatch(prefix, covnames)
     if (all(is.na(OK)))
        stop("prefix does not match any covariates")
    out <- vector('list')
    i <- 0
    repeat {
        xy <- unlist(locator(1))
        if (is.null(xy))
            break
        maskrow <- nearesttrap(xy, object)
        if (distancetotrap(xy,object[maskrow,]) > (2 * spacing(object)))
            break
        xy <- object[maskrow,,drop = FALSE]  ## centre
        if (point)
            points(xy, ...)
        D <- covariates(object)[maskrow,OK, drop = FALSE]
        for (j in 1:ncol(D)) {
            if (is.numeric(D[,j]))
                D[,j] <- format(round(D[,j]*scale, dec), nsmall = dec, trim = TRUE)
        }
        out[[i<-i+1]] <- cbind(xy,D)
        Dstr <- paste(D, collapse = sep)
        if (text)
            text(xy[1], xy[2], Dstr, ...)
        cat ('xy ', paste(round(xy,decxy), collapse = ', '), " : ", Dstr, '\n')
        flush.console()
    }
    invisible (do.call(rbind, out))
}
############################################################################################

getDensityArray <- function (x, paddedlength = NULL) {
    if (ms(x)) {
        ## find maximum mask points (may vary between sessions)
        maxnrow <- max(sapply(x, nrow))
        nsession <- length(x)
        densities <- lapply(x, getDensityArray, maxnrow)
        do.call(abind, densities)
    }
    else {
        covnames <- names(covariates(x))
        OK <- substring(covnames,1,2) == 'D.'
        covnames <- covnames[OK]
        densities <- covariates(x)[,covnames]
        if (!is.null(paddedlength))
            densities <- c(densities, rep(NA, paddedlength-length(densities)))
        array(densities, dim=c(nrow(x), length(covnames), 1))
    }
}

############################################################################################
