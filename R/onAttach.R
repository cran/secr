###############################################################################
## package 'secr'
## onAttach.R
## last changed 2011-06-16
###############################################################################

.onAttach <- function (libname, pkgname) {
    version <- packageVersion('secr')
    packageStartupMessage( "This is secr ", version,
        ". For overview type ?secr" )
# deprecated on r-devel 2011-06-12
#    cat(paste("This is secr ", version,
#       ". For overview type RShowDoc('secr-overview', package='secr')\n",
#        sep=""))
}
