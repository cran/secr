############################################################################################
## package 'secr'
## RPSV.R
## 2011 02 06 polygonX, transectX
## 2011 09 26 detector type check
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
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for RPSV")
        if (detector(traps) %in% .localstuff$polydetectors) {
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
