############################################################################################
## package 'secr'
## write.capthist.R
## last changed 2009 03 30, 2009 06 11, 2009 07 08 2009 11 17
## Write capture histories to text file in DENSITY format
############################################################################################

write.capthist <- function (object, file='', ..., deblank = TRUE, header = TRUE, sess = '1', ndec = 2) 

{
    objectname <-  deparse(substitute(object), control=NULL)

    if (!is(object, 'capthist')) stop ('write.capthist requires a capthist object')

    if (inherits(object, 'list')) {
        lapply (object, write.capthist, file = file, deblank = deblank, 
            header = FALSE, ...)  ## not finished
    }
    else {
        n <- nrow(object)
        S <- ncol(object)
        det <- detector(traps(object))
        if (det %in% c('polygon','transect')) {
            XY <- xy(object)            
            session <- rep(sess,nrow(XY))
            ID <- animalID(object) 
            occ <- occasion(object)
            temp <- data.frame (session=session, ID=ID, occ=occ, 
                x=round(XY$x,ndec), y=round(XY$y,ndec))
        }
        else if (det %in% c('signal')) {
            signal <- signal(object)
            session <- rep(sess, length(signal))            ## assumes 1 session 
            ID <- animalID(object) 
            occ <- occasion(object)
            trap <- trap(object)
            temp <- data.frame (session=session, ID=ID, occ=occ, trap=trap, signal=signal)
        }
        else if (det %in% c('proximity','count','quadratbinary','quadratcount')) {
                K <- dim(object)[3]
                session <- rep(rep(sess, n*S*K), abs(object))
                trap <- rep(rep(row.names(traps(object)), rep(n*S, K)), abs(object))
                ID <- rep(rep( row.names(object), S*K ), abs(object))
                occ <- rep(rep (rep (1:S, K), rep(n, S*K)), abs(object))
            temp <- data.frame (session=session, ID=ID, occ=occ, trap=trap)
        }
        else {    ## single, multi
            K <- 1
            OK <- abs(object) > 0
            trap <- row.names(traps(object))
            if (deblank) trap <- gsub(' ', '', trap)
            trap <- trap[abs(object)]
            session <- rep(sess, n*S*K)[OK]             ## assumes 1 session 
            ID <- rep( 1:n, S*K )[OK]
            occ <- rep (rep (1:S, K), rep(n, S*K))[OK]
            temp <- data.frame (session=session, ID=ID, occ=occ, trap=trap)
        }
   
        if (header) {
            cat ("; Capture histories exported from '", objectname, "' \n", sep="", file=file)
            cat (';', format(Sys.time(), "%a %b %d %X %Y"), '\n', append = TRUE, file=file)
            cat (';\n', append = TRUE, file=file)
            if ( detector(traps(object)) %in% c('polygon') )
                cat ('; Session ID Occasion   x   y \n', append = TRUE, file=file)
            else 
            if ( detector(traps(object)) %in% c('signal') )
                cat ('; Session ID Occasion Detector Signal\n', append = TRUE, file=file)
            else 
                cat ('; Session ID Occasion Detector\n', append = TRUE, file=file)
        }
        write.table(temp, file = file, append = header, 
            row.names = FALSE, col.names = FALSE, quote = FALSE, ...)
    
    }
}
############################################################################################

