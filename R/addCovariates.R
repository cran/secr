##############################################################################
## package 'secr'
## addCovariates.R
## 2011-11-01
## 2013-01-19 handles missing values in character fields
###############################################################################

addCovariates <- function (object, spatialdata, columns = NULL) {
    if (!(inherits(object, 'mask') | inherits(object, 'traps')))
        stop ("require mask or traps object")
    if (!ms(object) & ms(spatialdata))
        stop ("mismatch of single session object, multisession spatialdata")

    if (ms(object)) {
        ## allow multiple sessions, and session-specific data sources
        nsession <- length(object)
        out <- object
        for (session in 1:nsession) {
            if (ms(spatialdata))
                out[[session]] <- addCovariates(out[[session]], spatialdata[[session]])
            else
                out[[session]] <- addCovariates(out[[session]], spatialdata)
        }
        out
    }
    else {

        if (is.character(spatialdata))
            type <- "shapefile"
        else if (inherits(spatialdata, "SpatialPolygonsDataFrame"))
            type <- "SPDF"
        else if (inherits(spatialdata, "SpatialGridDataFrame"))
            type <- "SGDF"
        else if (inherits(spatialdata, "mask"))
            type <- "mask"
        else if (inherits(spatialdata, "traps"))
            type <- "traps"
        else
            stop ("spatialdata type unrecognised or unsupported")

        if (type == "shapefile") {
            ## if (!require (maptools))
            ##    stop ("package 'maptools' required for ", type)
            polyfilename <- spatialdata  ## strip shp?
            spatialdata <- maptools::readShapePoly(polyfilename)
        }
        if (type %in% c("shapefile", "SPDF", "SGDF")) {
            xy <- matrix(unlist(object), ncol = 2)
            xy <- SpatialPoints(xy)
            df <- over (xy, spatialdata)
        }
        else {
            ## nearest point algorithm
            if (is.null(covariates(spatialdata)))
                stop ("spatialdata does not have covariates")
            index <- nearesttrap(object, spatialdata)
            df <- covariates(spatialdata)[index,, drop=FALSE]

        }

        ## select requested columns
        if (!is.null(columns))
            df <- df[,columns, drop = FALSE]

        ## check new covariates OK
        fn <- function(x) {
            if (is.numeric(x))
                !any(is.na(x))
            else
                !any((x == "") | is.na(x))
        }
        OK <- all(apply(df, 2, fn))
        if (!OK)
            warning ("missing values among new covariates")

        ## insert new covariates and return object
        rownames(df) <- 1:nrow(df)
        if (is.null(covariates(object)))
            covariates(object) <- df
        else
            covariates(object) <- cbind(covariates(object), df)
        object
    }
}

