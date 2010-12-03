############################################################################################
## package 'secr'
## pdot.R
## return net detection probability in 'traps' for home ranges centred at X
## 2010 07 01 alpha detection fn
## 2010 10 09 extended for other detection functions
## 2010 10 10 esa.plot added
## 2010 10 11 esa.plot.secr added
## 2010 10 18 etc. pdot.contour, buffer.contour
## 2010 10 24 tweaking poly
## 2010 11 26 usage
## 2010 11 26 check for ms traps
############################################################################################

pdot <- function (X, traps, detectfn = 0, detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                  noccasions = 5) {

    ## X should be 2-column dataframe, mask, matrix or similar
    ## with x coord in col 1 and y coord in col 2

    ## added 2010-07-01
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)
    if ((detectfn > 9) & is.null(detectpar$cutval))
        stop ("requires 'cutval' for detectfn > 9")
    if (ms(traps))
        stop ("requires single-session traps")

    truncate <- ifelse(is.null(detectpar$truncate), 1e+10, detectpar$truncate)
    if (!is.null(usage(traps)))
        used <- unlist(usage(traps))
    else
        used <- rep(1, nrow(traps) * noccasions)

    X <- matrix(unlist(X),nc=2)
    if (detector(traps) == 'polygon') {
        k <- table(polyID(traps))
        temp <- .C('pdotpoly', PACKAGE = 'secr',
            as.double(X),
            as.integer(nrow(X)),
            as.double(unlist(traps)),
            as.integer(used),
            as.integer(ndetector(traps)),
            as.integer(k),
            as.integer(detectfn),   ## hn
            as.double(unlist(detectpar)),
            as.integer(noccasions),
            value = double(nrow(X))
        )
        temp$value
    }
    else if (detector(traps) == 'transect') {
        k <- table(transectID(traps))
        temp <- .C('pdottransect', PACKAGE = 'secr',
            as.double(X),
            as.integer(nrow(X)),
            as.double(unlist(traps)),
            as.integer(used),
            as.integer(ndetector(traps)),
            as.integer(k),
            as.integer(detectfn),
            as.double(unlist(detectpar)),
            as.integer(noccasions),
            value = double(nrow(X))
        )
        temp$value
    }
    else {
        temp <- .C('pdotpoint', PACKAGE = "secr",
            as.double(X),
            as.integer(nrow(X)),
            as.double(unlist(traps)),
            as.integer(used),
            as.integer(ndetector(traps)),
            as.integer(detectfn),
            as.double(unlist(detectpar)),
            as.integer(noccasions),
            as.double(truncate^2),
            value = double(nrow(X)))
        temp$value
    }
}

############################################################################################

esa.plot <- function (object, max.buffer = NULL, spacing = NULL, max.mask = NULL, detectfn,
                      detectpar, noccasions, thin = 0.1, poly = NULL, session = 1,
                      plt = TRUE, as.density = TRUE, n = 1, add = FALSE, overlay = TRUE, ...) {

    if (inherits(object, 'secr')) {
        esa.plot.secr (object, max.buffer, max.mask, thin, poly, session, plt,
                       as.density, add, overlay, ...)
    }
    else {

        if (!inherits(object, 'traps'))
            stop ("requires 'secr' or 'traps' object")
        args <- list(...)
        if(is.null(max.mask)) {
            if (is.null(spacing))
                 spacing <- spacing(object)/3
            max.mask <- make.mask (object, max.buffer, spacing,,  'trapbuffer', poly)
        }
        detectfn <- valid.detectfn(detectfn)
        a <- pdot (max.mask, object, detectfn, detectpar, noccasions)
        d <- distancetotrap(max.mask, object)
        ord <- order(d)
        cellsize <-  attr(max.mask, 'spacing')^2/10000
        a <- a[ord]
        output <- data.frame(buffer = d[ord], esa =  cumsum(a) * cellsize,
            density = n /  cumsum(a) / cellsize, pdot = a, pdotmin = cummin(a))
        maxesa <- max(output$esa)
        thinned <- seq(1,  nrow(max.mask), 1/thin)
        output <- output[thinned,]

        if (plt) {
            if (as.density) {
                if (add)
                    lines(output$buffer, n/output$esa, ...)
                else {
                    xlb <- 'Buffer width  m'
                    ylb <- expression(paste('n / esa(buffer)   ',ha^-1))
                    if ('ylim' %in% names(args))
                        plot(output$buffer, n/output$esa, type = 'l',
                            xlab = xlb, ylab = ylb, ...)
                    else  ## clunky!
                        plot(output$buffer, n/output$esa, type = 'l',
                            xlab = xlb, ylab = ylb, ...,
                            ylim= n / maxesa * c(0.9, 1.2))
                }
            }
            else {
                if (add)
                    lines(output$buffer, output$esa, ...)
                else
                    plot(output$buffer, output$esa, type = 'l',
                        xlab = 'Buffer width  m', ylab = 'esa(buffer)  ha', ...)
            }
            invisible(output)
        }
        else output
    }
}

############################################################################################

spatialscale <- function (object, detectfn, session = '') {
    if (inherits(object, 'secr')) {
        if (ms(object))
            detpar <- detectpar(object)[[session]]
        else
            detpar <- detectpar(object)
        cutval <- object$details$cutval
    }
    else {
        detpar <- object
        cutval <- object$cutval
    }
    if (!is.null(detpar$sigma)) detpar$sigma
    else if (detectfn == 10) {
        (cutval - detpar$beta0) / detpar$beta1
    }
    else if (detectfn == 11) {
        d11 <- function(d, beta0, beta1, c) beta0 + beta1 * (d-1) - 10 * log10(d^2) - c
        interval <- c(0,10 * (cutval - detpar$beta0) / detpar$beta1)
        uniroot (d11, interval, detpar$beta0, detpar$beta1, cutval)$root
    }
    else if (detectfn == 9) {
#        (0.5 - detpar$b0) / detpar$b1
        - 1 / detpar$b1   ## 2010-11-01
    }
    else stop ('unrecognised detectfn')
}

############################################################################################

esa.plot.secr <- function (object, max.buffer = NULL, max.mask = NULL, thin = 0.1,
                           poly = NULL, session = 1, plt = TRUE, as.density = TRUE,
                           add = FALSE, overlay = TRUE, ...) {

    if (!inherits(object,'secr'))
        stop('require secr object')

    MS <- ms(object)
    if (MS) {
        sessnames <- session(object$capthist)
        ## use alphanumeric session ID
        if (is.numeric(session))
            session <- sessnames[session]
    }

    ## recursive call
    if (MS & (length(session) > 1)) {
        esa.plot.outputs <- vector(mode='list')

        for (i in session) {
            addthisone <- ifelse (add | (overlay & (i != session[1])), TRUE, FALSE)
            esa.plot.outputs[[i]] <- esa.plot.secr (object, max.buffer, max.mask,
                thin, poly, i, plt, as.density, addthisone, overlay, ...)
        }
        if (plt)
            invisible(esa.plot.outputs)
        else
            esa.plot.outputs
    }
    ## not recursive
    else {
        if (MS) {
            ## select one session
            trps <- traps(object$capthist[[session]])
            n <- nrow(object$capthist[[session]])
            nocc <- ncol(object$capthist[[session]])
            spacg <- attr(object$mask[[session]], 'spacing')
            detpar <- detectpar(object)[[session]]
            spscale <- spatialscale(object, object$detectfn, session)
        }
        else {
            trps <- traps(object$capthist)
            n <- nrow(object$capthist)
            nocc <- ncol(object$capthist)
            spacg <- attr(object$mask, 'spacing')
            detpar <- detectpar(object)
            spscale <- spatialscale(object, object$detectfn)
        }
        if (is.null(max.mask)) {
            if (is.null(max.buffer)) {
                if (add)
                    max.buffer <- par()$usr[2]  ## width of existing plot
                else {
                    max.buffer <- 5 * spscale
                }
            }
        }
        esa.plot (trps, max.buffer, spacg, max.mask, object$detectfn, detpar,
                  nocc, thin, poly, session, plt, as.density, n, add, overlay, ...)
    }
}

############################################################################################

pdot.contour <- function (traps, border = NULL, nx = 64, detectfn = 0,
                          detectpar = list(g0 = 0.2, sigma = 25, z = 1),
                          noccasions = 5, levels = seq(0.1, 0.9, 0.1),
                          poly = NULL, plt = TRUE, add = FALSE, ...) {
    if (is.null(border))
        border <- 5 * spatialscale(detectpar, detectfn)
    tempmask <- make.mask (traps, border, nx = nx, type = 'traprect')
    xlevels <- unique(tempmask$x)
    ylevels <- unique(tempmask$y)
    z <- pdot(tempmask, traps, detectfn, detectpar, noccasions)
    if (!is.null(poly)) {
        OK <- insidepoly(tempmask, poly)
        z[!OK] <- 0
    }
    if (plt) {
        contour (xlevels, ylevels, matrix(z, nr = nx), add = add, levels = levels, ...)
        invisible(contourLines(xlevels, ylevels, matrix(z, nr = nx), levels = levels))
    }
    else
        contourLines(xlevels, ylevels, matrix(z, nr = nx), levels = levels)
}
############################################################################################

buffer.contour <- function (traps, buffer, nx = 64, convex = FALSE, ntheta = 100,
                            plt = TRUE, add = FALSE, poly = NULL, ...) {
    oneconvexbuffer <- function (buffer) {
        temp  <- data.frame(x = apply(expand.grid(traps$x, buffer * cos(theta)),1,sum),
                       y = apply(expand.grid(traps$y, buffer * sin(theta)),1,sum))
        temp <- temp[chull(temp), ]
        temp <- rbind(temp, temp[1,]) ## ensure closed
        rownames(temp) <- NULL
        if (plt)
            lines(temp,...)
        temp
    }
    if (!inherits(traps, 'traps'))
        stop("requires 'traps' object")
    if (convex) {
        if (!is.null(poly))
            warning("'poly' ignored when convex = TRUE")
        ## could use maptools etc. to get intersection?
        theta <- (2*pi) * (1:ntheta) / ntheta
        if (!add & plt)
            plot(traps, border = buffer)
        temp <- lapply(buffer, oneconvexbuffer)
        if (plt)
            invisible (temp)
        else
            temp
    }
    else {
        tempmask <- make.mask (traps, max(buffer)*1.2, nx = nx, type = 'traprect')
        xlevels <- unique(tempmask$x)
        ylevels <- unique(tempmask$y)
        z <- distancetotrap(tempmask, traps)
        if (!is.null(poly)) {
            OK <- insidepoly(tempmask, poly)
            z[!OK] <- 1e20
        }
        if (plt) {
            contour (xlevels, ylevels, matrix(z, nr = nx), add = add,
                 drawlabels = FALSE, levels = buffer,...)
            invisible(contourLines(xlevels, ylevels, matrix(z, nr = nx),
                levels = buffer))
        }
        else
            contourLines(xlevels, ylevels, matrix(z, nr = nx),
                levels = buffer)
    }
}
################################################################################
