\name{addTelemetry}
\alias{addTelemetry}

\title{Combine Telemetry and Detection Data}

\description{
  
Animal locations determined by radiotelemetry can be used to augment
capture--recapture data. The procedure in \pkg{secr} is first to form a
capthist object containing the telemetry data and then to combine this
with true capture--recapture data (e.g. detections from hair-snag DNA)
in another capthist object. \code{secr.fit} automatically detects the
telemetry data in the new object.

}

\usage{

addTelemetry (detectionCH, telemetryCH)

}

\arguments{
  \item{detectionCH}{single-session capthist object, detector type `proximity' or `count'}
  \item{telemetryCH}{single-session capthist object, detector type `telemetry'}
}

\details{

  It is assumed that a number of animals have been radiotagged in the
vicinity of the detector array, and their telemetry data
(xy-coordinates) have been input to \code{telemetryCH}, perhaps using
\code{\link{read.capthist}} with \code{detector = "telemetry"} and \code{fmt =
  "XY"}, or with \code{\link{read.telemetry}}.

A new capthist object is built comprising all the detection
histories in \code{detectionCH}, plus empty (all-zero) histories for
every telemetered animal not in \code{detectionCH}. The telemetry
locations are carried over from telemetryCH as attribute `xylist' (each
component of xylist holds the coordinates of one animal).

}

\value{
  A single-session capthist object with the same detector type as
  \code{detectionCH}, but possibly with empty rows and an `xylist' attribute.
}

\note{

Telemetry provides independent data on the location and presence of a
sample of animals. These animals may be missed in the main sampling that
gives rise to detectionCH i.e., they may have all-zero detection
histories.

The `telemetry' detector type is like a `polygon' detector (detections
have x-y coordinates). Although perimeter coordinates are required they
are not at present used in analyses.
  
Combining telemetry and detection data is new in \pkg{secr} 2.4.0, and
not yet fully documented.

}

\seealso{
  
  \code{\link{capthist}}, \code{\link{make.telemetry}}, \code{\link{read.telemetry}}
  
}

\examples{

\dontrun{

# Generate some detection and telemetry data, combine them using
# addTelemetry, and perform analyses

# detectors
te <- make.telemetry()
tr <- make.grid(detector = "proximity")

# simulated population and 50\% telemetry sample
totalpop <- sim.popn(tr, D = 20, buffer = 100)
tepop <- subset(totalpop, runif(nrow(totalpop)) < 0.5)

# simulated detection histories and telemetry
trCH <- sim.capthist(tr,  popn = totalpop, renumber = FALSE)
teCH <- sim.capthist(te, popn = tepop, renumber=FALSE,
    detectpar = list(g0 = 3, sigma = 25))

combinedCH <- addTelemetry(trCH, teCH)

# summarise and display
summary(combinedCH)
plot(combinedCH, border = 150)
ncapt <- apply(combinedCH,1,sum)
points(totalpop[row.names(combinedCH)[ncapt==0],], pch = 1)
points(totalpop[row.names(combinedCH)[ncapt>0],], pch = 16)

fit.tr <- secr.fit(trCH, CL = TRUE)                  ## trapping alone
fit.te <- secr.fit(teCH, CL = TRUE, start = log(20)) ## telemetry alone
fit2   <- secr.fit(combinedCH, CL = TRUE)            ## combined
fit2a   <- secr.fit(combinedCH, CL = TRUE,           ## combined, using info
    details = list(telemetrysigma = TRUE))           ## on sigma from telemetry

# improved precision when focus on realised population
# (compare CVD)
derived(fit.tr, distribution = "binomial")
derived(fit2, distribution = "binomial")

# may also use CL = FALSE

}

}

\keyword{ manip }
