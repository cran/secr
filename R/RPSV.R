############################################################################################
## package 'secr'
## RPSV.R
## last changed 2009 12 01
############################################################################################

RPSV <- function (capthist)
{
    if (inherits (capthist, 'list')) {
        lapply(capthist, RPSV)   ## recursive
    }
    else {

        traps <- traps(capthist)
        prox  <- length(dim(capthist)) > 2    ## includes 'count' etc
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
        if (prox)
        {   
            if (detector(traps) %in% c('polygon','transect')) {
                lxy <- split ( xy(capthist), animalID(capthist))
                temp <- lapply (lxy, RPSVxy)
            }
            else {
                w <- capthist != 0   
                w <- apply(w, 1:2,function(x) x * (1:length(x)))
                w <- aperm (w, c(2,3,1)) 
                lw <- lapply (1:nrow(capthist), function (i) rep(unlist(w[i,,]), as.numeric(abs(capthist[i,,]))))
                temp <- lapply(lw, RPSVx)
            }
            temp <- matrix(unlist(temp), nr=3)
            temp <- apply(temp,1,sum, na.rm=T)
            temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
        }
        else {
            temp <- apply(capthist, 1, RPSVx)     
            temp <- apply(temp[,temp[1,,drop=F]>0, drop=F], 1, sum)
            temp <- sqrt((temp[2]+temp[3]) / (temp[1]-1))
        }

        attr(temp,'names') <- NULL
        temp
    }
}
############################################################################################

