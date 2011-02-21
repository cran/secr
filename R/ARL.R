############################################################################################
## package 'secr'
## ARL.R
## last changed 2010 03 30
## polygonX, transectX, transect
############################################################################################

ARL <- function (capthist, min.recapt = 1, plt = FALSE, full = FALSE) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, ARL, plt = plt, full = full)   ## recursive
    }
    else {
        MMDMx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            if (sum(cx>0, na.rm=T) == 1) NA
            else {
              x <- traps$x[cx]
              y <- traps$y[cx]
              max(dist(cbind(x,y)))
            }
        }
        MMDMxy <- function (xy) {
            max(dist(cbind(xy$x, xy$y)))
        }
        traps <- traps(capthist)
        prox  <- length(dim(capthist)) > 2
        if (detector(traps) %in% c('polygon','transect','polygonX','transectX')) {
            lxy <- split (xy(capthist), animalID(capthist))
            maxd <- unlist(lapply (lxy, MMDMxy))
            n <- unlist(lapply (lxy, nrow))
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            maxd <- unlist(lapply(w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- try(coef(nls (maxd ~ aa * (1 - exp(bb * (n-1))),
            start= list (aa = max(maxd)*1.2, bb = -0.4))))

        if (inherits(temp, "try-error")) {
            warning ('nls failure')
            aa <- NA
            bb <- NA
        }
        else {
            aa <- temp[1]
            bb <- temp[2]
            if (plt) {
                plot (jitter(n, amount=0.1), maxd,
                    xlim=c(0,max(c(n,ncol(capthist)))),
                    xlab='Number of captures', ylab='ARL')
                xv <- seq(2,max(n),0.01)
                lines(xv, aa * (1 - exp(bb * (xv-1))))
            }
        }
        attr(aa,'names') <- NULL
        attr(bb,'names') <- NULL

        if (!full) aa
        else list (ARL = aa, b = bb, data = data.frame(maxd = maxd, n=n))
    }
}
############################################################################################

## source ('d:\\density secr 1.3\\secr\\r\\ARL.R')
##  data(Peromyscus)
##  ARL(Peromyscus.WSG, full=T)
