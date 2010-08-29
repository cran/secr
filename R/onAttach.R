############################################################################################
## package 'secr'
## onAttach.R
## last changed 2010 05 28
############################################################################################

.onAttach <- function (libname, pkgname) {
    ## with apologies to mgcv
    version <- library(help=secr)$info[[1]]
    version <- version[pmatch('Version',version)]
    um      <- strsplit(version, ' ')[[1]]
    version <- um [nchar(um) > 0][2]
    cat(paste("This is secr ", version,
              ". For overview type RShowDoc('secr-overview', package='secr')\n",
        sep=""))
}
