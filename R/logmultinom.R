############################################################################################
## package 'secr'
## logmultinom.R
## multinomial term in SECR likelihood
## last changed 2009 02 14, 2009 03 10, 2009 06 11, 2009 10 07
## verified vs DENSITY for 4-session dataset, island fox data set 2009 02 14
############################################################################################

logmultinom <- function (capthist, grp) {
## grp must not be NULL
    if (is.null(grp)) stop ('NULL grp in logmultinom')
    if (inherits (capthist, 'list')) {
        # assume grp also a list, by session
        nsession <- length(capthist) 
        lmn <- 0
        for (i in 1:nsession) { 
            lmn <- lmn + logmultinom (capthist[[i]], grp[[i]])
        }
        lmn
    }
    else {
        dim3 <- length(dim(capthist)) > 2  # 2009 10 07
        if (nrow(capthist)==0) stop(paste('no data for session'))
        nc <- nrow(capthist)

        # following distinguishes CH that end in death from CH that end with animal alive...
        # Count = number per unique capture history

        if (dim3) capthist <- matrix(capthist, nr=nc)  # here, 'nc' = number caught (rows)
        groupeddata <- split.data.frame(capthist, grp)
        ## need to protect make.lookup from bad data...
        count <- function(x) table(make.lookup(x)$index)
        counts <- sapply(groupeddata, count)
        sum(lgamma(table(grp)+1)) - sum(lgamma(unlist(counts)+1))
    }
}
############################################################################################

