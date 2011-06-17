############################################################################################
## package 'secr'
## derived.R
## 2010-10-22 allow object to be secrlist
## 2011-03-27 adjustments for zero capthist
############################################################################################

derived <- function (object, sessnum = NULL, groups=NULL, alpha=0.05, se.esa = FALSE,
    se.D = TRUE, loginterval = TRUE, distribution = NULL) {

## Generate table of derived parameters from fitted secr object

## modified 2009 07 21 for multiple sessions
## modified 2009 08 27 to report variance components

## multi-session -> separate each session, for the present
## is.null(sessnum) implies all possible sessions
## is.null(groups) implies all possible individuals

## multi-session -> pooled sessions not implemented

## groups within sessions
## groups to be found in covariates(object$capthist) cf


    if (inherits(object, 'secrlist')) {
        lapply(object, derived, sessnum, groups, alpha, se.esa, se.D,
               loginterval, distribution)
    }
    else {

        if (!is.null(distribution)) {
            if (tolower(distribution) %in% c('poisson','binomial'))
                object$details$distribution <- distribution
            else stop ("distribution not recognised")
        }

        if (is.list(object$capthist) & is.null(sessnum)) {
            ## recursive call if MS
            sessnames <- session(object$capthist)
            nsess <- length(sessnames)
            output <- vector('list', nsess )
            for (i in 1:nsess) {
                j <- match (sessnames[i], session(object$capthist))
                output[[i]] <- derived(object, sessnum = j, groups = groups,
                    alpha = alpha, se.esa = se.esa, se.D = se.D,
                    loginterval = loginterval, distribution = distribution)
            }
            names(output) <- sessnames
            output
        }
        else {

            maskarea <- function (mask, asess) {
               if (is.data.frame(mask)) nrow(mask) * attr(mask,'area')
               else nrow(mask[[asess]]) * attr(mask[[asess]],'area')
            }
            se.deriveD <- function (selection, selected.a, asess) {
                s2 <- switch (tolower(object$details$distribution),
                   poisson  = sum (1/selected.a^2),
                   binomial = sum (( 1 - selected.a / maskarea(object$mask, asess)) / selected.a^2))
                CLg  <- CLgradient (object, selection, asess)
                varDn <- CLg %*% object$beta.vcv %*% CLg
                list(SE=sqrt(s2 + varDn), s2=s2, varDn=varDn)
            }
            se.deriveesa <- function (selection, asess) {
                CLesa  <- esagradient (object, selection, asess)
                sqrt(CLesa %*% object$beta.vcv %*% CLesa)
            }
            weighted.mean <- function (a) {
                ## allows for varying a 2010-11-30
                length(a) / sum(1/a)
            }

            getderived <- function (selection) {
                    if (is.null(sessnum)) sessnum <- 1
                    selected.a <- esa(object, sessnum)[selection]
                    derivedmean <- c(weighted.mean(selected.a), sum(1/selected.a) )
                    derivedSE <- c(NA, NA)
                    varcomp1 <- c(NA, NA)
                    varcomp2 <- c(NA, NA)
                    if (se.esa) derivedSE[1] <- se.deriveesa(selection, sessnum)
                    if (se.D) {
                        varDlist <- se.deriveD(selection, selected.a, sessnum)
                        derivedSE[2] <- varDlist$SE
                        varcomp1[2] <- varDlist$s2
                        varcomp2[2] <- varDlist$varDn
                    }

                    temp <- data.frame (
                        row.names = c('esa','D'),
                        estimate = derivedmean,
                        SE.estimate = derivedSE)
                    nmash <- attr(capthist, 'n.mash')
                    temp <- add.cl(temp, alpha, loginterval)

                    temp$CVn <- varcomp1^0.5 / temp$estimate
                    temp$CVa <- varcomp2^0.5 / temp$estimate
                    temp$CVD <- temp$SE.estimate / temp$estimate
                    temp$CVD[1] <- NA   ## not for esa

                    if (!is.null(nmash)) {
                        temp[2,1:4] <- temp[2,1:4] / length(nmash)
                        ## message ("D was adjusted for ", length(nmash), " mashed clusters\n")
                    }

                    temp
            }

            if (is.null(sessnum))
                capthist <- object$capthist
            else
                capthist <- object$capthist[[sessnum]]
            grp <- group.factor(capthist, groups)  ## see functions.R
            if (nrow(capthist)>0)
                individuals <- split (1:nrow(capthist), grp)
            else
                individuals <-  split (numeric(0), grp) ## list of empty grp levels
            n <- length(individuals)
            if ( n > 1)
                lapply (individuals, getderived)
            else
                if (n == 1)
                    getderived(individuals[[1]])
                else
                    getderived(numeric(0))
        }
    }
}
############################################################################################

