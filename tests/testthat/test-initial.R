## Started 2021-01-30, 2021-05-18, 2021-05-28

library(secr)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

# create small working datasets
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, 
    dropnullocc = FALSE)))
nonspatialCH <- reduce(ovenCH, outputdetector = 'nonspatial', verify = FALSE)
msk <- make.mask(traps(smallCH), buffer = 200, nx = 20, type = 'trapbuffer')

argssecr <- list(capthist = smallCH, mask = msk, detectfn = 'HHN',
    start = list(lambda0 = 0.037, sigma = 65.3),
    details = list(LLonly = TRUE, fastproximity = FALSE))

test_that("correct likelihood (multi-session CL)", {
    argssecr$CL <- TRUE
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -301.1396, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct likelihood (multi-session CL fastproximity)", {
    argssecr$CL <- TRUE
    argssecr$details$fastproximity <- TRUE
    LL <- do.call(secr.fit, argssecr)[1]
    expect_equal(LL, -109.61114, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct likelihood (single session 'single')", {
    args <- list(capthist = captdata, buffer = 100, detectfn = 'HN',
        start = list(D = 5.4798, g0 = 0.2732, sigma = 29.3658),
        details = list(LLonly = TRUE))
    # 'multi' likelihood used for single-catch traps
    expect_warning(LL <- do.call(secr.fit, args)[1])  
    expect_equal(LL, -759.02575, tolerance = 1e-4, check.attributes = FALSE)
})

###############################
## Check bug fixes 2021 onwards

# make.capthist bug reported by Richard Glennie 2021-01-30
test_that("correct rejection of duplicates at exclusive detectors", {
    captures <- data.frame(session = c(1, 1, 1), 
        ID = c(1, 11, 1), 
        occasion = c(12, 2, 2),
        trap = c("A1", "A2", "A2"))
    traps <- make.grid(detector = "multi")
    expect_equal(nrow(make.capthist(captures, traps)), 2)
})


###############################
## Multi-polygon bug 2021-05-18

datadir <- system.file("extdata", package = "secr")
polyexample1 <- read.traps(file = paste0(datadir, '/polygonexample1.txt'), 
    detector = 'polygon')
polygonCH <- sim.capthist(polyexample1, popn = list(D = 1, buffer = 200),
    detectfn = 'HHN', detectpar = list(lambda0 = 5, sigma = 50),
    noccasions = 1, seed = 123)

test_that("simulated polygon data have correct RPSV", {
    rpsv <- RPSV(polygonCH, CC = TRUE)
    expect_equal(rpsv, 45.16401, tolerance = 1e-4, check.attributes = FALSE)
})

test_that("correct likelihood (multi-detector polygon data)", {
    args <- list(capthist = polygonCH, buffer = 200, detectfn = 'HHN',
        start = list(D=1, lambda0=5, sigma = 50), verify = FALSE, 
        details = list(LLonly = TRUE))
    LL <- do.call(secr.fit, args)[1]
    expect_equal(LL, -2026.92169, tolerance = 1e-4, check.attributes = FALSE)
})

#################################
## join nonspatial bug 2021-05-28

nonspatialCH <- reduce(ovenCH, outputdetector = 'nonspatial', verify = FALSE)
test_that("join works on nonspatial data", {
    sumcapt <- sum(join(nonspatialCH))
    expect_equal(sumcapt, 190, tolerance = 1e-4, check.attributes = FALSE)
})


###############################
## future tests

## data manipulation

# subset
# reduce
# join
# rbind

## for data types
# point  
# polygon
# transect
# signal

## using RPSV as test criterion?

## check alongtransect for multi-transect data

## 2021-05-19 known bug: sim.capthist with multiple transects
