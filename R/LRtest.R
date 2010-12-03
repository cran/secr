############################################################################################
## package 'secr'
## LRtest.R
## likelihood ratio test for two secr models (assumed to be nested)
## last changed 2009 06 11
############################################################################################

LR.test <- function (secr1, secr2) {
    if (class(secr1) != 'secr' | class(secr2) != 'secr')
        stop ("requires two fitted 'secr' models")
    statistic <- 2 * abs(secr1$fit$value - secr2$fit$value)
    if (length(statistic) != 1)
        stop ("problem with 'secr1' or 'secr2'")
    parameter <- abs(length(secr1$fit$par) - length(secr2$fit$par))
    names(statistic) <- 'X-square'
    names(parameter) <- 'df'
    s1name <- deparse(substitute(secr1))
    s2name <- deparse(substitute(secr2))
    temp <- list (
        statistic = statistic,
        parameter = parameter,
        p.value = 1 - pchisq(statistic, parameter),
        method = 'Likelihood ratio test for two SECR models',
        data.name = paste (s1name, 'vs', s2name)
    )
    class(temp) <- 'htest'
    temp
}
############################################################################################

