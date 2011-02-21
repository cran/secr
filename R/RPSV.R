############################################################################################
## package 'secr'
## RPSV.R
## last changed 2010 03 30 streamlined
## 2011 02 06 polygonX, transectX

############################################################################################

RPSV <- function (capthist)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, RPSV)   ## recursive
    }
    else {
        RPSVx <- function (cx) {
            cx <- abs(cx)  # no special trt for deads
            x <- traps$x[cx]
            y <- traps$y[cx]
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        RPSVxy <- function (xy) {
            x <- xy$x
            y <- xy$y
            n <- length(x)
            c(n = n-1, ssx = sum(x^2) - (sum(x))^2/n, ssy = sum(y^2) - (sum(y))^2/n)
        }
        traps <- traps(capthist)
        if (detector(traps) %in% c('polygon','transect','polygonX','transectX')) {
            lxy <- split ( xy(capthist), animalID(capthist))
            temp <- lapply (lxy, RPSVxy)
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            temp <- lapply(w,RPSVx)
        }
        temp <- matrix(unlist(temp), nrow = 3)
        temp <- apply(temp,1,sum, na.rm=T)
        temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
        attr(temp,'names') <- NULL
        temp
    }
}
############################################################################################
