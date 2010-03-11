############################################################################################
## package 'secr'
## write.traps.R
## last changed 2009 06 11 2009 11 17
## Write detector locations to text file in DENSITY format
## should remove conflict between row and ...
############################################################################################

write.traps <- function (object, file='', ..., deblank = TRUE, header = TRUE, ndec = 2) {
    objectname <-  deparse(substitute(object), control=NULL)
    if (!is(object, 'traps')) stop ('write.traps requires a traps object')
    n <- nrow(object)
    object$x <- round(object$x,ndec)
    object$y <- round(object$y,ndec)

    # purge blanks from names
    if (deblank) row.names(object) <- gsub(' ','',row.names(object))

    poly <- detector(object) %in% c('polygon')
    transect <- detector(object) %in% c('transect')
    if (poly) {
        temp <- cbind (polyID=polyID(object), x=object$x, y = object$y)
    }
    else if (transect) {
        temp <- cbind (transectID=transectID(object), x=object$x, y = object$y)
    }
    else {
        temp <- object
        if (!is.null(usage(object))) temp <- cbind(temp,usage(object))
    }

    if (header) {
        cat ("; Detector locations exported from '", objectname, "' \n", sep="", file=file)
        cat (';', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE, file=file)
        cat (';\n', append = TRUE, file=file)
        if (poly)
            cat ('; polyID  x  y\n ', append = TRUE, file=file)
        else
        if (transect)
            cat ('; transectID  x  y\n ', append = TRUE, file=file)
        else
            cat ('; detector  x  y\n', append = TRUE, file=file)
    }
    write.table(temp, file = file, append = header, 
        row.names = !poly & !transect, col.names = FALSE, quote = FALSE, ...)

}
############################################################################################

