## Started 2021-01-30

library(secr)

# create small working dataset
suppressWarnings(smallCH <- join(subset(ovenCHp, sessions = 1:3, traps = 1:20, 
    dropnullocc = FALSE)))
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
    LL <- do.call(secr.fit, args)[1]
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
