############################################################################################
## package 'secr'
## LRtest.R
## likelihood ratio test for two models (assumed to be nested)
## last changed 2011 12 13 to include openCR
############################################################################################

LR.test <- function (model1, model2) {
    if (!((inherits(model1, 'secr') | inherits(model1,'openCR')) &
        (inherits(model2, 'secr') | inherits(model1,'openCR'))))
        stop ("requires two fitted 'secr' or 'openCR' models")
    statistic <- 2 * abs(model1$fit$value - model2$fit$value)
    if (length(statistic) != 1)
        stop ("problem with 'model1' or 'model2'")
    parameter <- abs(length(model1$fit$par) - length(model2$fit$par))
    names(statistic) <- 'X-square'
    names(parameter) <- 'df'
    s1name <- deparse(substitute(model1))
    s2name <- deparse(substitute(model2))
    temp <- list (
        statistic = statistic,
        parameter = parameter,
        p.value = 1 - pchisq(statistic, parameter),
        method = 'Likelihood ratio test for two models',
        data.name = paste (s1name, 'vs', s2name)
    )
    class(temp) <- 'htest'
    temp
}
############################################################################################

