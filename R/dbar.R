############################################################################################
## package 'secr'
## dbar.R
## last changed 2009 07 10 (multi-session)  2009 10 07 (prox) 2009 11 13 (polygon)
## streamlined 2010 03 30
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
        if (detector(traps) == 'polygon') {
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

