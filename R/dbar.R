############################################################################################
## package 'secr'
## dbar.R
## 2009 07 10 (multi-session)
## 2009 10 07 (prox)
## 2009 11 13 (polygon)
## 2010 03 30 streamlined
## 2010 06 04 added function moves
## 2011 02 04 extended for transect data
## 2011 02 06 polygonX, transectX
## 2011 09 26 detector type check
## 2011 12 20 try() to make more robust

############################################################################################

dbar <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, dbar)   ## recursive
    }
    else {
        traps <- traps(capthist)
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for dbar")
        dbarx    <- function (x) {
            x <- abs(unlist(x))
            sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
        }
        dbarxy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        if (detector(traps) %in% .localstuff$polydetectors) {
            lxy <- split (xy(capthist), animalID(capthist))
            d <- try(lapply(lxy,dbarxy), silent = TRUE)
            if (inherits(d, 'try-error'))
                d <- NA
            mean(unlist(d), na.rm=T)
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            d <- try(unlist(lapply(w,dbarx)), silent = TRUE)
            if (inherits(d, 'try-error'))
                d <- NA
            mean(d, na.rm=T)
        }
    }
}
############################################################################################

moves <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, moves)   ## recursive
    }
    else {
        movex    <- function (x) {
            x <- abs(unlist(x))
            sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
        }
        movexy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        traps <- traps(capthist)
        if (!(detector(traps) %in% .localstuff$individualdetectors))
            stop ("require individual detector type for moves")
        if (detector(traps) %in% .localstuff$polydetectors) {
            lxy <- split (xy(capthist), animalID(capthist))
            lapply (lxy, movexy)
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            lapply(w,movex)
        }
    }
}
############################################################################################

