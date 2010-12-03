########################################################################################
## package 'secr'
## sim.popn.R
## simulate spatially distributed population
## transferred from methods.R 2010-06-08
## last changed 2010-06-14  added Ndist = 'specified', using D to carry N
#########################################################################################

toroidal.wrap <- function (pop) {
    bb <- attr(pop, 'boundingbox')
    xmin <- min(bb$x)
    xmax <- max(bb$x)
    ymin <- min(bb$y)
    ymax <- max(bb$y)
    xrange <- xmax-xmin
    yrange <- ymax-ymin
    remainder <- function (x,y) x - (x %/% y * y)
    pop$x <- ifelse (pop$x>xmax, xmin + remainder (pop$x-xmax, xrange), pop$x)
    pop$x <- ifelse (pop$x<xmin, xmax - remainder (xmin - pop$x, xrange), pop$x)
    pop$y <- ifelse (pop$y>ymax, ymin + remainder (pop$y-ymax, yrange), pop$y)
    pop$y <- ifelse (pop$y<ymin, ymax - remainder (ymin - pop$y, yrange), pop$y)
    pop
}

sim.popn <- function (D, core, buffer = 100, model2D = 'poisson',
    buffertype = 'rect', covariates = list(sex=c(M=0.5,F=0.5)),
    number.from = 1, Ndist = 'poisson', nsession = 1, details = NULL,
    seed = NULL)
{
    if (nsession > 1) {
        discrete <- function(x) {
            fr <- x-trunc(x)
            sample (c(trunc(x), trunc(x)+1), size=1, prob=c(1-fr, fr))
        }
        session.popn <- function (s) {
            ## independent population
            sim.popn (D, core, buffer, model2D, buffertype,
                covariates, number.from, Ndist, nsession = 1, details, seed)
        }
        turnover <- function (oldpopn) {
            ## project existing population
            ## assume normal movement kernel
            ## assume lambda lacks process variance
            ## ideally lambda lognormal
            ## need 'wrap' option for toroidal wrapping of 'rect' locations
            newstart <- max(as.numeric(rownames(oldpopn))) + 1
            if (turnoverpar$survmodel=='binomial') {
                survive <- sample (c(FALSE, TRUE), nrow(oldpopn), replace = TRUE,
                    c(1-turnoverpar$phi,turnoverpar$phi))
                nsurv <- sum(survive)
            }
            else {   ## assume 'discrete'
                nsurv <- discrete (turnoverpar$phi * nrow(oldpopn))
                survive <- sample (nrow(oldpopn), replace = FALSE, size = nsurv)
                survive <- sort(survive)   ## numeric indices
            }
            newpopn <- subset.popn(oldpopn, subset=survive)
            newpopn[,] <- newpopn[,] + rnorm (2*nsurv, mean = 0, sd =
                turnoverpar$sigma.m)
            if (turnoverpar$wrap)
                newpopn <- toroidal.wrap(newpopn)

            gam <- turnoverpar$lambda - turnoverpar$phi
            if (gam<0)
                stop ("invalid gamma in turnover")
            nrecruit <- switch (turnoverpar$recrmodel,
                constantN = nrow(oldpopn) - nsurv,
                discrete = discrete(gam * nrow(oldpopn)),
                poisson = rpois (1, gam * nrow(oldpopn)))
                ## under Pradel model members of superpopulation have binomial
                ## probability of entry at this point?
                ## cf Schwarz & Arnason betas
            if (nrecruit>0) {
                recruits <- sim.popn(D = nrecruit, core = core, buffer = buffer,
                    model2D = model2D, buffertype = buffertype, covariates =
                    covariates, number.from = newstart, Ndist = 'specified',
                    nsession = 1, details = details)  ## danger: resets random seed
                newpopn <- rbind.popn(newpopn, recruits, renumber = FALSE)
            }
            attr(newpopn, 'losses') <- nrow(oldpopn)-nsurv
            attr(newpopn, 'recruits') <- nrecruit
            newpopn
        }
        turnoverpar <- list(lambda = NULL, phi = 0.7, sigma.m = 20, wrap = TRUE,
                            survmodel = 'binomial', recrmodel = 'poisson')
        turnoverpar <- replace (turnoverpar, names(details), details)
        if (is.null(details$lambda)) {
            ## independent
            MSpopn <- lapply (1:nsession, session.popn)
        }
        else {
            ## projected
            MSpopn <- vector(nsession, mode = 'list')
            MSpopn[[1]] <- session.popn(1)
            for (i in 2:nsession) {
                MSpopn[[i]] <- turnover(MSpopn[[i-1]])
            }
        }
        class(MSpopn) <- c('list','popn')
        names(MSpopn) <- 1:nsession
        MSpopn
    }
    else {
        ##########################
        ## set random seed
        ## copied from simulate.lm

        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        ##########################

        if (model2D == 'IHP') {
            nr <- nrow(core)
            if (!inherits(core, 'mask'))
                stop("for model2D = IHP, 'core' should be a habitat mask")
            if (Ndist != 'poisson')
                stop ("IHP not implemented for fixed or specified N")
            nm <- rpois(nr, D * attr(core,'area'))   ## 'area' is cell area, D vector, 1 per mask cell
            N <- sum(nm)
            jitter <- matrix ((runif(2*sum(nm))-0.5) * attr(core,'spacing'), nc=2)
            animals <- core[rep(1:nr, nm),] + jitter
            animals <- as.data.frame(animals)
            xl <- range(animals[,1])
            yl <- range(animals[,2])
        }
        else {
            if (buffertype != 'rect')
                stop (dQuote("rect"), " is only buffertype implemented")
            # population in arena +/- buffer from traps
            buff <- c(-buffer,+buffer)
            xl <- range(core$x) + buff
            yl <- range(core$y) + buff
            area <- diff(xl) * diff(yl) * 0.0001  # ha
            allowedNdist <- c('poisson','fixed','specified')
            if (!(Ndist %in% allowedNdist))
                stop ("'Ndist' should be one of ",
                      paste(sapply(allowedNdist, dQuote), collapse=","))
            N  <- switch (Ndist,
                poisson = rpois(1, lambda=D[1] * area),
                fixed = discreteN (1, D[1] * area),
                specified = D[1])

            if (model2D=='poisson') {
                animals <- data.frame (x = runif(N)*diff(xl)+xl[1], y =
                    runif(N)*diff(yl)+yl[1])
            }
            else if (model2D=='cluster') {
                ## Neyman-Scott distribution with wrapping
                xrange <- diff(xl)
                yrange <- diff(yl)
                if (details$mu<=0) {
                    nparent <- N   ## not clustered
                    offspr <- sweep(matrix(runif(2*nparent),nc=2),2,c(xrange,yrange),'*')
                }
                else {
                     nparent <- switch (Ndist,
                         poisson = rpois(1, lambda=D[1] * area/details$mu),
                         fixed = discreteN (1, D[1] * area / details$mu),
                         specified = discreteN (1, D[1] / details$mu))  ## here arg D is N
                     N <- nparent * details$mu
                     if (nparent==0)
                         warning ("zero clusters")
                     parent <-  sweep(matrix(runif(2*nparent),nc=2),2,c(xrange,yrange),'*')

                     offspr <- matrix(rnorm(2*N),nc=2) * details$hsigma
                     parentn <- rep(1:nparent, details$mu)
                     offspr <- offspr + parent[parentn,]
                     while (any ((offspr[,1]<0) | (offspr[,1]>xrange) | (offspr[,2]<0) |
                                 (offspr[,2]>yrange))) {
                       offspr[,1] <- ifelse (offspr[,1]<0, offspr[,1]+xrange, offspr[,1])
                       offspr[,1] <- ifelse (offspr[,1]>xrange, offspr[,1]-xrange, offspr[,1])
                       offspr[,2] <- ifelse (offspr[,2]<0, offspr[,2]+yrange, offspr[,2])
                       offspr[,2] <- ifelse (offspr[,2]>yrange, offspr[,2]-yrange, offspr[,2])
                     }
                }
                animals <- as.data.frame(sweep(offspr,2,c(xl[1],yl[1]),'+'))
            }
            else stop ("unrecognised 2-D distribution")
        }
        names(animals) <- c('x','y')
        row.names (animals) <- number.from : (nrow(animals)+number.from-1)
        attr(animals,'covariates') <- NULL
        if (!is.null(covariates)) {
            tempcov <- list()
            for (i in 1:length(covariates)) {
               covi <- sample (names(covariates[[i]]), replace=T, size=N, prob=covariates[[i]])
               temptxt <- paste ('tempcov$', names(covariates[i]), '<- covi', sep='')
               eval(parse(text=temptxt))
            }
            attr(animals,'covariates') <- as.data.frame(tempcov)
        }
        attr(animals, 'seed') <- RNGstate   ## save random seed
        attr(animals, 'Ndist') <- Ndist
        attr(animals, 'model2D') <- model2D
        attr(animals, 'boundingbox') <- expand.grid (x=xl,y=yl)[c(1,3,4,2),]
        class(animals) <- c('popn', 'data.frame')
        animals
    }
}

