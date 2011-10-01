
predictD.secr <- function (object, newdata = NULL, se.fit = FALSE, alpha = 0.05,
    savenew = FALSE, ...) {

    if (is.null(object$fit)) {
        warning ("empty (NULL) object")
        return(NULL)
    }
    if (object$CL) {
        stop ("requires full-likelihood model")
    }
    vars <- object$vars

    if (is.null(newdata)) {
        nr <- nrow(object$mask)
        newdata <- data.frame(rep(1,nr))
        for (v in vars) {
            if (v=='x') newdata$x <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
            if (v=='y') newdata$y <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'y']
            if (v=='x2') newdata$x2 <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
            if (v=='y2') newdata$y2 <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'y']
            if (v=='xy') newdata$xy <- rep(0,nr)   # mean attr(mask,'meanSD')[1,'x']
        }
    }

    ## allow for fixed beta parameters
    fb <- object$details$fixedbeta
    if (!is.null(fb)) {
        nbeta <- length(fb)
        beta.vcv <- matrix(0, nrow = nbeta, ncol = nbeta)
        beta.vcv[is.na(fb[row(beta.vcv)]) & is.na(fb[col(beta.vcv)])] <- object$beta.vcv
        fb[is.na(fb)] <- object$fit$par
        beta <- fb    ## complete
    }
    else {
        beta <- object$fit$par
        beta.vcv <- object$beta.vcv
    }

    model <- object$model$D
    indx <- object$parindx$D

    vars <- all.vars(model)
    if (any(!(vars %in% names(newdata))))
        stop ("one or more model covariates not found in 'newdata'")
    newdata <- as.data.frame(newdata)

    lpred <- matrix(ncol = 2, nrow = nrow(newdata),
        dimnames = list(NULL,c('estimate','se')))
    mat <- model.matrix(model, data=newdata)
    lpred[,1] <- mat %*% beta[indx]

    predict  <- cbind(newdata,lpred)
    if (se.fit) {
        vcv <- beta.vcv[indx,indx]    ## OR maybe all betas?
        nrw <- 1:nrow(mat)
        vcv <- apply(expand.grid(nrw, nrw), 1, function(ij)
            mat[ij[1],, drop=F] %*% vcv %*% t(mat[ij[2],, drop=F]))  # link scale
        vcv <- matrix (vcv, nrow = nrw)
        lpred[,2] <- diag(vcv)^0.5
        predict <- cbind(newdata,lpred)
        attr(predict, 'vcv') <- vcv
    }
    z <- abs(qnorm(1-alpha/2))   ## beware confusion with hazard z!

    lpred  <- predict[,'estimate']
    Xlpred <- untransform(lpred, object$link$D)

#    if (se.fit) {
#        out <- list(nrow(newdata))
#        selpred <- predict[,'se']
#        temp <- data.frame (
#          row.names = object$realnames,
#          link = unlist(object$link[object$realnames]),
#              estimate = Xlpred,
#              SE.estimate = se.Xuntransform (lpred, selpred, object$link, object$realnames),
#              lcl = Xuntransform(lpred-z*selpred, object$link, object$realnames),
#              ucl = Xuntransform(lpred+z*selpred, object$link, object$realnames)
#              )
#            # truncate density at zero
#            temp['D', -1][temp['D',-1]<0] <- 0
#            if (nrow(newdata)==1) out <- temp
#            else {
#                out[[new]] <- temp
#                names(out)[new] <- paste (
#                        paste(names(newdata),'=', unlist(lapply(newdata[new,],as.character)),
#                        sep=' ',collapse=', '),
#                    sep=',')
#            }
#        }
#        else { # no SE; terse format

    out <- data.frame(D=Xlpred)


    if (savenew) attr(out, 'newdata') <- newdata
    out
}
