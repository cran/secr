############################################################################################
## package 'secr'
## units.R
## started speculatively 2010-10-23
############################################################################################

## defaults
.localstuff$linearUnits <- 'm'
.localstuff$areaUnits <- 'ha'
.localstuff$validLinearUnits <- c('m', 'km')
.localstuff$validAreaUnits <- c('ha', 'km^2', '100km^2')

units <- function (object, ...) UseMethod("units")

units.default <- function (object, ...) {
    attr(object,'units')
}

'units<-' <- function (object, value)
    structure (object, units = value)

convertUnits <- function (x, linearUnits, areaUnits = NULL) {
    if (is.null(areaUnits))
        areaUnits <- switch(linearUnits, m = 'ha', km = 'km^2')
    if (!(linearUnits %in% .localstuff$validLinearUnits))
        stop ("unrecognised linearUnits")
    if (!(areaUnits %in% .localstuff$validAreaUnits))
        stop ("unrecognised areaUnits")

    ## scale numeric values

    ## change 'units' attribute
    units(x) <- list(linearUnits, areaUnits)
    x
}

compatibleUnits <- function (...) {
    objects <- list(...)  ##?
    if (length(objects)==1)
        TRUE
    else {
        lin <- sapply(objects, function(x) units(x)$linear)  ##? linearUnits
        are <- sapply(objects, function(x) units(x)$area)
        ! (any(diff(lin)>0) | (any(diff(are)>0)))
    }
}

## Affected objects
##
## capthist
## traps
## mask
## secr, secrlist
##
##
## Affected functions
##
## autoini D, sigma, mask
## print.secr
## plot.secr
## make.grid
## make.mask
## read.capthist
## secr.make.newdata
##
## ARL, MMDM, dbar, RPSV?
##
## esa.plot
## detectfnplot
## polyarea
##
##
## New functions needed
##
## Operational issues
##
## combining objects with incompatible units
##    automatically coerce with warning
