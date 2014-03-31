clone <- function (object, type, ...) UseMethod("clone")

clone.default <- function (object,  type, ...)       {
    if (length(dim(object)) != 2)
        stop ("requires 2-D object")
    type <- tolower (type)
    n <- nrow(object)
    if (n == 0) {
        out <- object
    }
    else {
        if (type == 'constant')
            freq <- rep(..., n)
        else if (type == 'poisson')
            freq <- rpois(n, ...)
        else if (type == 'nbinom')
            freq <- rnbinom(n, ...)
        else
            stop("unrecognised type")
        index <- rep(1:n, freq)
        object[index,]
    }
}

clone.popn <- function (object, type, ...) {
    if (ms(object)) {
        out <- lapply (object, clone, type, ...)
        class (out) <- c('popn','list')
        out
    }
    else {
        type <- tolower (type)
	n <- nrow(object)
        if (n == 0) {
            out <- object
        }
        else {
            if (type == 'constant')
                freq <- rep(..., n)
            else if (type == 'poisson')
                freq <- rpois(n, ...)
            else if (type == 'nbinom')
                freq <- rnbinom(n, ...)
            else
                stop("unrecognised type")
            index <- rep(1:n, freq)
            out <- object[index,]
            if (!is.null(covariates(object)))
                covariates(out) <- covariates(object)[index,]
            attr (out, 'freq') <- freq
        }
        out
    }
}
