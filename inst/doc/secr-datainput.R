## ----eval = FALSE-------------------------------------------------------------
#  read.capthist(captfile, trapfile, detector = "multi", fmt = c("trapID", "XY"),
#      noccasions = NULL, covnames = NULL, trapcovnames = NULL, cutval = NULL,
#      verify = TRUE, noncapt = "NONE", ...)

## ----message=FALSE------------------------------------------------------------
library(secr)
captfile <- system.file("extdata", "stoatcapt.txt", package = "secr")
trapfile <- system.file("extdata", "stoattrap.txt", package = "secr")
stoatCH <- read.capthist(captfile, trapfile, detector = "proximity")
summary(stoatCH)

## ----envvar, echo = FALSE-----------------------------------------------------
## Following is not needed as no multithreaded operations in this vignette 
## To avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

## ----readxl, warning = FALSE--------------------------------------------------
xlsname <- system.file("extdata", "stoat.xlsx", package = "secr")
CH <- read.capthist (xlsname, sheet = c("stoatcapt", "stoattrap"), skip = 1,
                  detector = "proximity")
summary(CH)

## -----------------------------------------------------------------------------
write.capthist(signalCH, "temp")  ## export data for demo
tempCH <- read.capthist("tempcapt.txt", "temptrap.txt", detector = "signal", cutval = 52.5)

## ----eval = FALSE-------------------------------------------------------------
#  read.capthist("captXY.txt", "perimeter.txt", fmt = "XY", detector = "polygon")

## ----eval = FALSE-------------------------------------------------------------
#  temppoly <- read.traps(file = "clipboard", detector = "polygon")
#  tempcapt <- sim.capthist(temppoly, popn = list(D = 1, buffer = 1000), detectpar =
#                             list(g0 = 0.5, sigma = 250))
#  plot(tempcapt, label = TRUE, tracks = TRUE, title = "Simulated detections within polygons")

## -----------------------------------------------------------------------------
summary(subset(stoatCH, traps = 1:47, occasions = 1:5))

