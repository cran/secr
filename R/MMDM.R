############################################################################################
## package 'secr'
## MMDM.R
## last changed 2010 03 30
############################################################################################

MMDM <- function (capthist, min.recapt = 1, full = FALSE) {
    if (inherits (capthist, 'list')) {
        lapply(capthist, MMDM, full = full)   ## recursive
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
        MMDMxy    <- function (xy) {
                max(dist(cbind(xy$x, xy$y)))
        }
        traps <- traps(capthist)
        if (detector(traps) %in% c('polygon','transect','polygonX','transectX')) {
            lxy <- split (xy(capthist), animalID(capthist))
            maxd <- unlist(lapply (lxy, MMDMxy))
            n <- unlist(lapply (lxy, nrow))
        }
        else {
            ## streamlined 2010 03 30
            w <- split(trap(capthist, names=F), animalID(capthist))
            maxd <- unlist(lapply( w, MMDMx))
            n <- unlist(lapply(w, length))
        }
        maxd <- maxd[n>min.recapt]
        n <- n[n>min.recapt]

        temp <- mean (maxd, na.rm = TRUE)

        if (!full) temp
        else {
            SE <- function(x) sqrt(var(x, na.rm=T)/sum(!is.na(x)))
            summaryd <- data.frame (Ncapt = names(table(n)),
                           n = as.numeric(table(n)),
                           mean = tapply(maxd, n, mean, na.rm=T),
                           se = tapply(maxd, n, SE))
            summaryd$mean[is.na(summaryd$mean)] <- NA  ## tidier
            summaryd <- rbind(summaryd, data.frame(Ncapt='Total', n=sum(!is.na(maxd)),
                mean=temp, se=SE(maxd)))
            list (MMDM = temp, data = data.frame(maxd = maxd, n=n), summary = summaryd )
        }
    }
}
############################################################################################

## source ('d:\\density secr 1.3\\secr\\r\\mmdm.R')
## data(Peromyscus)

## MMDM(Peromyscus.WSG, full=T)$summary
##   Ncapt  n     mean        se
## 1     1  9       NA        NA
## 2     2  9 28.32839  9.434631
## 3     3 10 24.05921  9.335062
## 4     4  8 33.87949  5.471227
## 5     5  8 52.37655 15.470420
## 6     6  7 34.24929  5.961350
## 7 Total 42 33.93669  4.495646

## MMDM(Peromyscus.WSG, full=T)$summary[,3:4]/15.2
##       mean        se
## 1       NA        NA
## 2 1.863710 0.6206994
## 3 1.582843 0.6141488
## 4 2.228914 0.3599492
## 5 3.445826 1.0177908
## 6 2.253243 0.3921941
## 7 2.232677 0.2957662   <<<< 0.575?

## cf Otis et al 1978 p 87 Fig 23a
