############################################################################################
## package 'secr'
## dbar.R
## last changed 2009 07 10 (multi-session)  2009 10 07 (prox) 2009 11 13 (polygon)
## streamlined 2010 03 30
## added function moves 2010 06 04
## extended for transect data 2011 02 04
## polygonX, transectX 2011 02 06
############################################################################################

dbar <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, dbar)   ## recursive
    }
    else {
        traps <- traps(capthist)
        dbarx    <- function (x) {
            x <- abs(unlist(x))
            sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
        }
        dbarxy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        if (detector(traps) %in% c('polygon','polygonX', 'transect','transectX')) {
            lxy <- split (xy(capthist), animalID(capthist))
            mean(unlist(lapply (lxy, dbarxy)), na.rm=T)
        }
        else {
            w <- split(trap(capthist, names=F), animalID(capthist))
            mean(unlist(lapply(w,dbarx)), na.rm=T)
        }
    }
}
############################################################################################

moves <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, moves)   ## recursive
    }
    else {
        traps <- traps(capthist)
        movex    <- function (x) {
            x <- abs(unlist(x))
            sqrt(diff(traps$x[x])^2 + diff(traps$y[x])^2)
        }
        movexy    <- function (xy) {
            sqrt(diff(xy$x)^2 + diff(xy$y)^2)
        }
        if (detector(traps) %in% c('polygon', 'transect', 'polygonX', 'transectX')) {
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

