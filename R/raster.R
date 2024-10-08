## register S3 classes as (virtual) S4 classes
setOldClass('mask')
setOldClass('Dsurface')

## define 'raster' methods for mask and Dsurface
setMethod("raster", signature(x = "mask"),
    ## define internal function for coercing mask to raster
    function(x,  covariate, values = 1, crs = NA) {
        if (missing(covariate)) {
            covariate <- 'tempcov'
            covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
        }
        mask <- rectangularMask(x)
        bbox <- attr(x, 'boundingbox')
        xmx <- max(bbox$x)
        xmn <- min(bbox$x)
        ymx <- max(bbox$y)
        ymn <- min(bbox$y)
        nx <- round((xmx - xmn) / spacing(mask))
        ny <- round((ymx - ymn) / spacing(mask))
        tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
        mat <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
        r <- raster(mat, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs, template = NULL)
        return(r)
    }
)

setMethod("raster", signature(x = "Dsurface"),
    ## define internal function for coercing Dsurface to raster
    function(x,  covariate, values = 1, crs = NA) {
        if (missing(covariate)) {
            covariate <- 'tempcov'
            covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
        }
        mask <- rectangularMask(x)
        bbox <- attr(x, 'boundingbox')
        xmx <- max(bbox$x)
        xmn <- min(bbox$x)
        ymx <- max(bbox$y)
        ymn <- min(bbox$y)
        nx <- round((xmx - xmn) / spacing(mask))
        ny <- round((ymx - ymn) / spacing(mask))
        tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
        x <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
        raster(x, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs)
    }
)

## checks
## getMethod('raster','mask')
## getMethod(raster, 'matrix')
## plot(raster(predictDsurface(fit.Dforest), 'D.0'), useRaster = F)

## 2022-01-24 terra SpatRaster
setMethod("rast", signature(x = "mask"),
    ## define internal function for coercing mask to SpatRaster
    function(x,  covariate, values = 1, crs = "") {
        if (missing(covariate)) {
            covariate <- 'tempcov'
            covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
        }
        mask <- rectangularMask(x)
        tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
        df <- cbind(mask, tmp)
        rast(df, type = 'xyz', crs = crs)
    }
)

setMethod("rast", signature(x = "Dsurface"),
    ## define internal function for coercing Dsurface to SpatRaster
    function(x,  covariate, values = 1, crs = "") {
        if (missing(covariate)) {
            covariate <- 'tempcov'
            covariates(x) <- data.frame(tempcov = rep(values, length.out=nrow(x)))
        }
        mask <- rectangularMask(x)
        tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
        df <- cbind(mask, tmp)
        rast(df, type = 'xyz', crs = crs)
    }
)

