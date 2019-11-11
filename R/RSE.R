## extract relative standard error from fitted secr or openCR model
## default is vector for all parameters in model
RSE <- function (fit, parm = NULL, newdata = NULL) {
    if (is.null(parm)) {
        parm <- names(fit$parindx)
    }
    if (length(parm) > 1) {
        sapply(parm, RSE, fit = fit, newdata = newdata)
    }
    else {
        ## condition ensures 
        ## (i) parm in model (ii) no modelled variation in parm
        pindex <- length(fit$parindx[[parm]])
        if (pindex == 0) {
            warning ("parameter ", parm, " not in model")
        }
        else if (pindex == 1 && fit$link[[parm]] == 'log') {
            if (!is.null(newdata)) {
                warning ("non-null newdata ignored for ", parm)
            }
            ## Efford & Boulanger MEE 2019
            sqrt(exp(vcov(fit)[parm, parm])-1)    
        }
        else {
            if (pindex > 1 && is.null(newdata)) {
                warning ("parameter ", parm, " varies in model; consider specifying newdata")
            }
            pred <- predict(fit, newdata = newdata)
            if (parm %in% rownames(pred)) {
                pred[parm, 'SE.estimate'] / pred[parm, 'estimate']
            }
            else if (parm %in% rownames(pred[[1]])) {
                pred[[1]][parm, 'SE.estimate'] / pred[[1]][parm, 'estimate']
            }
            else if (parm %in% names(pred)) {
                pred[[parm]][1,'SE.estimate'] / pred[[parm]][1,'estimate']
            }
            else NA
        }
    }
}