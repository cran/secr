###############################################################################
## package 'secr'
## writeGPS.R
## last changed 2011-04-24
## Write detector locations to GPS format
##
###############################################################################

## writeGPS basics
## o = 'garmin' is only for Garmin usb/serial protocol
## use o = 'gdb' for MapSource

writeGPS <- function (xy, o = "garmin", F = "usb:", proj = '+proj=nzmg')  {

    if (proj=="") {
        latlon <- xy
    }
    else {
        ## express as lat long
        if (!require (rgdal))
            stop ("package 'rgdal' required for writeGPS")
        xy <- as.matrix(xy)  ## required by rgdal, not proj4
        latlon <- project (xy, inv = TRUE, proj)
    }

    latlon <- data.frame(latlon)[,c(2,1)] ## need lat first
    latlon$label <- rownames(xy)
    tempf <- tempfile('waypts')
    old <- options(digits=12)  ## ensure plenty of digits
    write.table(latlon, quote = FALSE, file = tempf,
        col = FALSE, row = FALSE, sep = ',')

    ## adapted from maptools `readGPS'
    GB <- Sys.which("gpsbabel")
    if (nchar(GB) == 0 || !file.exists(GB))
        stop ("gpsbabel not found")
    cmd <- paste(GB, "-D9 -w -i csv -f ", tempf, " -o ", o, " -F ", F)
    if (.Platform$OS.type == "windows")
        gpsdata <- system(cmd, intern = TRUE, invisible = TRUE)
    else
        gpsdata <- system(cmd, intern = TRUE)

    file.remove(tempf)
    options(old)
    invisible()

}