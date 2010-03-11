############################################################################################
## package 'secr'
## dbar.R
## last changed 2009 07 10 (multi-session)  2009 10 07 (prox) 2009 11 13 (polygon)
############################################################################################

dbar <- function (capthist) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, dbar)   ## recursive
    }
    else {
        traps <- traps(capthist)
        prox  <- length(dim(capthist)) > 2
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
            if (prox) { # replace '1' with trap number
                    w <- capthist != 0   ## 2009 10 07 transform to 0,1
                    w <- apply(capthist,1:2,function(x) x * (1:length(x)))
                    w <- aperm(w, c(2,1,3))
                    mean(unlist(apply(w,1,dbarx)), na.rm=T)
            }
            else mean(unlist(apply(capthist,1,dbarx)), na.rm=T)
        }
    }
}
############################################################################################

