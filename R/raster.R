
## register S3 classes as (virtual) S4 classes
setOldClass('mask')
setOldClass('Dsurface')

## define 'raster' methods for mask and Dsurface
setMethod("raster", signature(x = "mask"),
          ## define internal function for coercing mask to raster
          function(x,  covariate, values = 1, ...) {
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
            nx <- (xmx - xmn) / spacing(mask)
            ny <- (ymx - ymn) / spacing(mask)
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            x <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
            raster(x, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=NA)
          }
)
setMethod("raster", signature(x = "Dsurface"),
          ## define internal function for coercing Dsurface to raster
          function(x,  covariate, values = 1, ...) {
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
            nx <- (xmx - xmn) / spacing(mask)
            ny <- (ymx - ymn) / spacing(mask)
            tmp <- as.numeric(unlist(covariates(mask)[,covariate] ))
            x <- matrix (tmp, nrow = ny, ncol = nx, byrow = TRUE)[ny:1,]
            raster(x, xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=NA)
          }
)

## checks
## getMethod('raster','mask')
## getMethod(raster, 'matrix')
## plot(raster(predictDsurface(fit.Dforest), 'D.0'), useRaster = F)
