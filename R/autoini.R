############################################################################################
## package 'secr'
## autoini.R
## last changed 2009 04 21, 2009 07 12 (failure R2.9)
## 2009 07 17 add f.lower, f.upper to uniroot call for small speed improvement
## 2009 10 27 changes to dbar for robustness with counts etc.
## 2009 11    RPSV
## 2009 12 11 update integralprw1 call for nmix=1
## 2010 02 14 update integralprw1 call for spherical (no longer used)
############################################################################################

autoini <- function (capthist, mask, detectfn = 0, thin = 0.2) 

# obtain approximate fit of HN SECR model
# for use as starting values in MLE
# uses external C code

{
    naivedcall <- function (sigma) 
    {      
        temp <- .C ("naived",  PACKAGE = 'secr',
          as.double (sigma),               # Parameter : detection scale
          as.integer(k),               
          as.integer(m),
          as.double (unlist(traps)),       # x,y locations of traps (first x, then y) 
          as.double (unlist(mask)),        # x,y points on mask (first x, then y) 
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value 
        )
        db-temp$value
    }
    naiveRPSVcall <- function (sigma) 
    {      
        temp <- .C ("naiveRPSV",  PACKAGE = 'secr',
          as.double (sigma),               # Parameter : detection scale
          as.integer(k),               
          as.integer(m),
          as.double (unlist(traps)),       # x,y locations of traps (first x, then y) 
          as.double (unlist(mask)),        # x,y points on mask (first x, then y) 
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value 
        )
        if (temp$value > 0)
            obsRPSV - temp$value
        else
            obsRPSV                        # dummy value  
    }
    naivecap2 <- function (g0, sigma, cap) 
    # using new algorithm for number of captures per animal 2009 04 21
    {      
        temp <- .C ("naivecap2",  PACKAGE = 'secr',
          as.integer(prox),
          as.double (g0),                  # Parameter : detection magnitude
          as.double (sigma),               # Parameter : detection scale
          as.integer(s),                   # number of occasions 
          as.integer(k),                   # number of traps
          as.integer(m),                   # number of points in mask
          as.double (unlist(traps)),       # x,y locations of traps (first x, then y) 
          as.double (unlist(mask)),        # x,y points on mask (first x, then y) 
          as.integer (detectfn),           # code 0 = halfnormal
          value = double(1)                # return value 
        )
        cap - temp$value
    }
    naiveesa <- function (g0, sigma) 
    {      
      nc <- 1
      g0sigma0 <- matrix(rep(c(g0,sigma),c(2,2)),nr=2)
      gs0 <- rep(1,2*s*k)
      area <- attr(mask,'area')  # area of single mask cell

      temp <- try ( .C("integralprw1",  PACKAGE = 'secr',
          as.integer(dettype),
          as.double(g0sigma0),
          as.integer(nc), 
          as.integer(s), 
          as.integer(k), 
          as.integer(m), 
          as.integer(1),
          as.double(unlist(traps)), 
          as.double(unlist(mask)), 
          as.integer(nrow(g0sigma0)),  # rows in lookup
          as.integer(gs0), # index of nc+1,S,K to rows in g0sigma0
          as.integer(1),
          as.double(area), 
          as.integer(detectfn), 
          as.integer(0),     # binomN
          as.double(0),      # cut
          a=double(nc),
          resultcode=integer(1))
      )

      if (temp$resultcode != 0) stop ('Error in external function integralprw1; possibly the mask is too large', call.=F)
      temp$a
    }

    if (nrow(capthist)<5) stop ('Too few values in session 1 to determine start; set manually')
    if (detectfn!=0) stop ('Only halfnormal detection function (0) implemented in autoini')
    traps    <- traps(capthist)
    dettype  <- detectorcode(traps)
    if (!(dettype %in% c(-1:5,8))) list(D=NA,g0=NA,sigma=NA)   ## 2009 11 18
    else {
        prox     <- detector(traps) %in% c('proximity', 'count','quadratbinary','quadratcount','signal','times')
        n        <- nrow(capthist)    # number of individuals
        s        <- ncol(capthist)    # number of occasions
        k        <- nrow(traps)
        # optionally thin mask 
        if (thin>0 & thin < 1)
            mask <- mask[runif(length(mask$x)) < thin,]
        else 
            thin <- 1.0
        m        <- length(mask$x)    # number of points in mask
        if (length(dim(capthist))>2)
            cpa     <- sum(abs(capthist))/n      ## captures per animal    
        else
            cpa     <- sum(abs(capthist)>0)/n    ## captures per animal    

        # if (cpa<=1.0) stop ('No recaptures')

        obsRPSV <- RPSV(capthist)

        if (is.na(obsRPSV) | (obsRPSV<1e-10)) {    ## try db
            db <- dbar(capthist)
            if (!is.null(attr(traps,'spacing'))) {
                if (is.na(db)) {
                    warning ('could not calculate dbar; using detector spacing')       
                    db <- attr(traps, 'spacing')
                }
                if (db < (attr(traps, 'spacing')/100)) {
                    warning ('dbar close to zero; using detector spacing instead')       
                    db <- attr(traps, 'spacing')
                }
            }
            if (is.na(db) | is.nan(db) | (db<1e10) ) 
                return(list(D=NA,g0=NA,sigma=NA))   
            else  
                tempsigma <- uniroot (naivedcall, lower = db/10, upper = db*10, tol=0.001)$root
        }    
        else {
            tempsigma <- uniroot (naiveRPSVcall, lower = obsRPSV/10, upper = obsRPSV*10, tol=0.001)$root
        }

        low <- naivecap2(0.001, sigma=tempsigma, cap=cpa)
        upp <- naivecap2(0.999, sigma=tempsigma, cap=cpa)
        badinput <- FALSE
        if (is.na(low) | is.na(upp)) badinput <- TRUE
        else if (sign(low) == sign(upp)) badinput <- TRUE
        if (badinput) {
            # not sure what conditions cause this 28/4/08
            warning ('autoini failed to find g0; setting initial g0 = 0.1')
            tempg0 <- 0.1
        }
        else tempg0 <- uniroot (naivecap2, lower=0.001, upper=0.999, 
                 f.lower = low, f.upper = upp, tol=0.001,  
                 sigma=tempsigma, cap=cpa)$root
        esa       <- naiveesa (tempg0, tempsigma)
        list(D= n / esa * thin, g0 = tempg0, sigma = tempsigma)
    }    
}
##################################################################################
