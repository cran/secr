############################################################################################
## plot.secr.R
## Method for plotting detection function from fitted secr object
## 2009 02 12 2009 08 09 2009 09 13 (limits)
## 2009 11 04 signal strength
############################################################################################

plot.secr <- function (x, newdata=NULL, add = FALSE, 
    sigmatick = FALSE, rgr = FALSE, limits = TRUE, alpha = 0.05, xval = 0:200, 
    ylim = NULL, xlab = NULL, ylab = NULL, ...) 
{
    HN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]
        g0 * exp (-r^2 / 2 / sigma^2)
    }

    HZ <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * (1 - exp (-(r / sigma)^-z))
    }

    EX <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * exp (-r / sigma)
    }
    UN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]
        ifelse (r<=sigma, g0, 0)
    }
    CHN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        1 - (1 - g0 * exp (-r^2 / 2 / sigma^2)) ^ z
    }
    WEX <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
        ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
    }
    SS <- function (pars, r) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(x$details$cutval)) 
            stop ('require cut parameter for signal strength plot')
        mu <- beta0 + beta1 * r
        1 - pnorm (q=x$details$cutval, mean=mu, sd=sdS)
    }
    SSS <- function (pars, r) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(x$details$cutval)) 
            stop ('require cut parameter for signal strength plot')
        ## spherical so assume distance r measured from 1 m
        mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
        mu[r<1] <- beta0
        1 - pnorm (q=x$details$cutval, mean=mu, sd=sdS)
    }
  
    gline <- function (predicted, rowi = 1, eps = 1e-10) {
        ## eps is used to limit y to range where gradient() works
        ## may need later adjustment
        if (!is.data.frame(predicted)) {
            out <- list()
            for (i in 1:length(predicted)) 
                out[[i]] <- gline(predicted[[i]], i)
            names(out) <- names(predicted)
            out
        } 
        else {
        
            parnames <- switch (x$detectfn+1, 
                c('g0','sigma'),
                c('g0','sigma','z'),
                c('g0','sigma'),
                c('g0','sigma'),
                c('g0','sigma','z'),
                c('g0','sigma','w'),
                ,,,,
                c('beta0','beta1', 'sdS'),
                c('beta0','beta1', 'sdS')
            )
    
            pars <- predicted[parnames,'estimate']  
      
            pars[is.na(pars)] <- unlist(x$fixed)
            dfn <- switch (x$detectfn+1, HN, HZ, EX, UN, CHN, WEX,,,,,SS,SSS)   ## omits binary SS!
            if (sigmatick) {
              sigma <- pars[2]
              y <- dfn(pars, sigma)
              dy <- par()$usr[4]/20
              segments (sigma, y-dy, sigma, y+dy)
            }

            y <- dfn(pars, xval)
  
            if (rgr) {
              y <- xval * y
              ymax <- par()$usr[4]
              lines (xval, y * 0.8 * ymax / max(y), lty = 2, ...)  
            }
            else lines (xval, y, ...)
    
            if (limits & !rgr) {
                ## delta method variance of g()
    
                grad <- matrix(nr=length(xval), nc=length(x$fit$par))  ## beta pars
                if (is.null(newdata)) newdata <- secr.make.newdata (x)
                lkdfn <- function (beta, r) {
                    ## real g() from all beta pars and model.matrix 
                    real <- numeric(length(parnames))
                    names(real) <- parnames
                    for (rn in parnames) {
                         par.rn <- x$parindx[[rn]]
                         mat <- model.matrix(x$model[[rn]], data=newdata[rowi,,drop=F])
                         lp <- mat %*% matrix(beta[par.rn], nc = 1)
                         real[rn] <- untransform (lp, x$link[[rn]])
                    }
                    ## dfn(real, r)   # on natural scale
                    logit(dfn(real, r))
                } 

                for (i in 1:length(xval)) 

            ## Fast special cases: checking only
            ##    if ((x$detectfn==0) & (all(sapply(x$model, function(m) m == ~1)) ))
            ##    {
            ##        ## ASSUME DEFAULT LINK
            ##        g0 <- logit(pars[1])
            ##        sigma <- log(pars[2])
            ##        r <- xval[i]
            ##        ## D(expression(1/(1+exp(-g0)) * exp(-r^2/2/exp(sigma)^2)), 'g0')
            ##        ## D(expression(1/(1+exp(-g0)) * exp(-r^2/2/exp(sigma)^2)), 'sigma')
            ##        tempgrad <- c(exp(-g0)/(1 + exp(-g0))^2 * exp(-r^2/2/exp(sigma)^2),
            ##            -(1/(1 + exp(-g0)) * (exp(-r^2/2/exp(sigma)^2) * (-r^2/2 * (2 * (exp(sigma) 
            ##            * exp(sigma)))/(exp(sigma)^2)^2))))
            ##        if (!x$CL) tempgrad <- c(0,tempgrad)  ## for density
            ##        grad[i,] <- tempgrad
            ##    } 
            ##    else
            ##       if ((x$detectfn==1) & (all(sapply(x$model, function(m) m == ~1)) ))
            ##    {
            ##        ## ASSUME DEFAULT LINK
            ##        g0 <- logit(pars[1])
            ##        sigma <- log(pars[2])
            ##        z <- log(pars[3])
            ##           r <- xval[i]
            ##        tempgrad <- c(
            ##            exp(-g0)/(1 + exp(-g0))^2 * (1 - exp(-(r/exp(sigma))^(-exp(z)))),
            ##            -(1/(1 + exp(-g0)) * (exp(-(r/exp(sigma))^(-exp(z))) * ((r/exp(sigma))^((-exp(z)) - 
            ##            1) * ((-exp(z)) * (r * exp(sigma)/exp(sigma)^2))))),
            ##            -(1/(1 + exp(-g0)) * (exp(-(r/exp(sigma))^(-exp(z))) * ((r/exp(sigma))^(-exp(z)) * 
            ##            (log((r/exp(sigma))) * exp(z)))))
            ##        )
            ##        if (!x$CL) tempgrad <- c(0,tempgrad)  ## for density
            ##        grad[i,] <- tempgrad
            ##    } 
            ##    else

                # grad[i,] <- fdHess (pars = x$fit$par, fun = lkdfn, r = xval[i])$gradient
                # grad[i,] <- grad (func = lkdfn, x = x$fit$par, r = xval[i])  ## needs numDeriv
                grad[i,] <- gradient (pars = x$fit$par, fun = lkdfn, r = xval[i])  ## see 'functions.R'

                vc <- vcov (x)
                se <- apply(grad, 1, function(gg) { gg <- matrix(gg,nr=1); gg %*% vc %*% t(gg) })^0.5

                ## lcl <- pmax(y - z*se,0)  # on natural scale
                ## ucl <- pmin(y + z*se,1)

                ## limits on link scale
                lcl <- ifelse ((y>eps) & (y<(1-eps)), invlogit (logit(y) - z*se), NA) 
                ucl <- ifelse ((y>eps) & (y<(1-eps)), invlogit (logit(y) + z*se), NA)

                lines (xval, lcl, lty=2, ...)            
                lines (xval, ucl, lty=2, ...)            
            }
    
            if (limits & !rgr)
                data.frame(x=xval, y=y, lcl = lcl, ucl = ucl)
            else 
                data.frame(x=xval, y=y)
        }
    }

    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard-rate z!
    temp <- predict (x, newdata)
    if (is.null(ylim)) {
        if (x$detectfn %in% c(10,11)) {
            ylim <- c(0, 1)
        }
        else { 
            getmax <- function(x) {
                g0 <- x['g0','estimate']
                se.g0 <- x['g0','SE.estimate']
                if (limits) min(1, g0 + z * se.g0) else g0
            }
            if (is.data.frame(temp)) maxg0 <- getmax(temp)
            else maxg0 <- max(sapply (temp, getmax))

            if (is.na(maxg0)) maxg0 <- x$fixed$g0
            if (maxg0 > 0.75) maxg0 <- max(maxg0,1)          
            ylim <- c(0, maxg0)
        }
    }
    if (!add) {
        if (is.null(xlab)) 
            xlab <- 'Distance  (m)'
        if (is.null(ylab)) {
           binomN <- ifelse(is.null(x$details$binomN),0,x$details$binomN)
           dlambda <- (detector(traps(x$capthist)) %in% c('quadratbinary', 'polygon')) |
               ((detector(traps(x$capthist)) %in% c('count','quadratcount')) & (binomN==0))
           if (dlambda) 
               ylab <- 'Detection lambda'
           else
               ylab <- 'Detection probability'
        }      
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim, 
            xlab=xlab, ylab=ylab,...)
    }
    invisible(gline(temp))
}
############################################################################################

detectfnplot <- function (detectfn, pars, details = NULL,
    add = FALSE, sigmatick = FALSE, rgr = FALSE,  
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...) 
{
    HN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]
        g0 * exp (-r^2 / 2 / sigma^2)
    }

    HZ <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * (1 - exp (-(r / sigma)^-z))
    }

    EX <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        g0 * exp (-r / sigma)
    }
    UN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]
        ifelse (r<=sigma, g0, 0)
    }
    CHN <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; z <- pars[3]
        1 - (1 - g0 * exp (-r^2 / 2 / sigma^2)) ^ z
    }
    WEX <- function (pars, r) {
        g0 <- pars[1]; sigma <- pars[2]; w <- pars[3]
        ifelse(r<=w, g0, g0*exp (-(r-w) / sigma))
    }
    SS <- function (pars, r) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(details$cutval)) 
            stop ('require cut parameter for signal strength plot')
        mu <- beta0 + beta1 * r
        1 - pnorm (q=details$cutval, mean=mu, sd=sdS)
    }
  
    SSS <- function (pars, r) {
        beta0 <- pars[1]; beta1 <- pars[2]; sdS <- pars[3]
        if (is.null(details$cutval)) 
            stop ('require cut parameter for signal strength plot')
        ## spherical so assume distance r measured from 1 m
        mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
        mu[r<1] <- beta0
        1 - pnorm (q=details$cutval, mean=mu, sd=sdS)
    }
  
    gline <- function (pars) {
        ## here pars is a vector of parameter values

        dfn <- switch (detectfn+1, HN, HZ, EX, UN, CHN, WEX,,,,,SS,SSS)
        if (sigmatick) {
            sigma <- pars[2]
            y <- dfn(sigma, pars)
            dy <- par()$usr[4]/20
            segments (sigma, y-dy, sigma, y+dy)
        }
        y <- dfn(pars, xval)
        if (rgr) {
            y <- xval * y
            ymax <- par()$usr[4]
            lines (xval, y * 0.8 * ymax / max(y), lty = 2, ...)  
        }
        else lines (xval, y, ...)

        data.frame(x=xval, y=y)

    }

    ### mainline

    if (!is.matrix(pars)) pars <- matrix(pars, nr=1)

    needp <- c(2,3,2,2,3,3,0,0,0,0,3,3)[detectfn+1]
    nam <- c('halfnormal','hazard','exponential','uniform','compound halfnormal',
        'w exponential',"","","","", 'signal strength', 'binary signal strength')
    if (ncol(pars) != needp)
        stop(paste('require', needp, 'parameters for', nam[detectfn+1], 'detection function'))

    if (is.null(ylim)) {
        if (detectfn %in% c(10,11)) {
            ylim <- c(0, 1)
        }
        else { 
            ylim <- c(0, max(pars[,1]))   ## g0
        }
    }

    if (!add) {
        if (is.null(xlab)) 
            xlab <- 'Distance  (m)'
        if (is.null(ylab))
            ylab <- 'Detection'
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim, 
            xlab=xlab, ylab=ylab,...)
    }

    invisible( apply (pars, 1, gline) )

}
############################################################################################

attenuationplot <- function (pars, add = FALSE, spherical = TRUE,
    xval = 0:200, ylim = NULL, xlab = NULL, ylab = NULL, ...) {

    mufn <- function (pars, r) {
        beta0 <- pars[1]
        beta1 <- pars[2]
        ## if spherical, assume distance r measured from 1 m
        if (spherical) {
            mu <- beta0 - 10 * log ( r^2 ) / 2.302585 + beta1 * (r-1)
            mu[r<1] <- beta0
        }
        else
            mu <- beta0 + beta1 * r
        mu
    }

    aline <- function (pars, r) {
        y <- mufn (pars, xval)
        lines (xval, y, ...)
        data.frame(x=xval, y=y)
    }

    if (!is.matrix(pars)) pars <- matrix(pars, nr=1)

    if (is.null(ylim)) {
        lower <- min(apply(pars, 1, mufn, r = max(xval)))
        ylim <- c(lower, max(pars[,1]))
    }

    if (!add) {
        if (is.null(xlab)) 
            xlab <- 'Distance  (m)'
        if (is.null(ylab))
            ylab <- 'Acoustic power (dB)'
        plot (type ='n', 0,0, xlim=range(xval), ylim=ylim, 
            xlab=xlab, ylab=ylab,...)
    }
    invisible( apply (pars, 1, aline) )

}

