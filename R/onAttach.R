############################################################################################
## package 'secr'
## onAttach.R
## last changed 2010 12 02
############################################################################################

.onAttach <- function (libname, pkgname) {
    version <- packageVersion('secr')
    cat(paste("This is secr ", version,
              ". For overview type RShowDoc('secr-overview', package='secr')\n",
        sep=""))
}
