## Started 2022-01-05

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

####################################################################

set.seed(123)
detectors <- make.grid (nx = 6, ny = 6, detector = "count")
CH1 <- sim.capthist (detectors, popn = list(D = 0.5, buffer = 100), 
    detectfn = 'HHN', detectpar = list(lambda0 = 0.2, sigma = 25), 
    noccasions = 4, nsessions = 3)

test_that("reduce.capthist accepts zero, one, or more animals per session", {
    sum1 <- summary(CH1, terse = TRUE)
    CH2 <- reduce(CH1, outputdetector = "count", by = "all", 
        verify = FALSE, dropunused = FALSE)
    sum2 <- summary(CH2, terse = TRUE)
    expect_equal(sum1[2:4,], sum2[2:4,], tolerance = 1e-4, check.attributes = FALSE)
})

####################################################################
