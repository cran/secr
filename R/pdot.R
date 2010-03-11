############################################################################################
## package 'secr'
## pdot.R
## return net detection probability in 'traps' for home ranges centred at X
## last changed 2009 07 03
############################################################################################

## TO DO: implement other df using form of plot.secr

pdot <- function (X, traps, detectfn = 0, detectpar = list(g0 = 0.2, sigma = 20), noccasions = 5) {
    ## X should be 2-column dataframe, mask, matrix or similar 
    ## with x coord in col 1 and y coord in col 2
    dk <- function (X) apply(traps,1,function(xk) sum((X - xk)^2)^0.5)
    pdothn <- function (X) 1 - prod(1 - detectpar$g0 * exp(-dk(X)^2 / 2 / detectpar$sigma^2)) ^ noccasions
    pdothz <- function (X) 1 - prod(1 - detectpar$g0 * (1 - exp (-(dk(X) / detectpar$sigma)^-detectpar$z))) ^ noccasions
    
    if (detectfn != 0) stop ('Only hn implemented for now...')
    X <- matrix(unlist(X),nc=2)
    if (detector(traps) == 'polygon') {
       k <- table(polyID(traps))
       temp <- .C('pdotpoly', PACKAGE = 'secr', 
           as.double(X),
           as.integer(nrow(X)),
           as.double(unlist(traps)),
           as.integer(ndetector(traps)),
           as.integer(k),
           as.integer(detectfn),   ## hn
           as.integer(noccasions),
           as.double(c(detectpar$g0, detectpar$sigma,1)),
           value = double(nrow(X))
       )
       temp$value
    }
    else if (detector(traps) == 'transect') {
       k <- table(transectID(traps))
       temp <- .C('pdottransect', PACKAGE = 'secr', 
           as.double(X),
           as.integer(nrow(X)),
           as.double(unlist(traps)),
           as.integer(ndetector(traps)),
           as.integer(k),
           as.integer(detectfn),   ## hn
           as.integer(noccasions),
           as.double(c(detectpar$g0, detectpar$sigma,1)),
           value = double(nrow(X))
       )
       temp$value
    }
    else {
        apply(X, 1, pdothn)
    }  
}

# pdot(c(0,0), traps=temptrap)
# tempmask <- make.mask(temptrap)
# plot(subset(tempmask, pdot (tempmask, traps=temptrap, detectpar=list(g0=0.2, sigma=20), noccasions=5) > 0.01))

