\name{capthist.parts}
\alias{animalID}
\alias{occasion}
\alias{trap}
\alias{alive}
\alias{xy}
\alias{signalframe}
\alias{signal}
\alias{noise}
\alias{alongtransect}
\alias{xy<-}
\alias{signalframe<-}
\alias{signal<-}
\alias{noise<-}

\title{ Dissect Spatial Capture History Object }

\description{
  Extract parts of an object of class `capthist'.
}

\usage{

animalID(object, names = TRUE)
occasion(object)
trap(object, names = TRUE)
alive(object)
alongtransect(object, tol = 0.01)
xy(object)
xy(object) <- value
signalframe(object)
signalframe(object) <- value
signal(object)
signal(object) <- value
noise(object)
noise(object) <- value
}

\arguments{

  \item{object}{ a `capthist' object }
  \item{names}{ if FALSE the values returned are numeric indices rather
    than names }
  \item{tol}{ tolerance for snapping to transect line (m) }
  \item{value}{ replacement value (see Details) }

}
\details{

These functions extract data on detections, ignoring occasions when an
animal was not detected. Detections are ordered by occasion, animalID
and trap.

\code{trap} returns polygon or transect numbers if \code{traps(object)}
has detector type `polygon' or `transect'.

\code{alongtransect} returns the distance of each detection from the
start of the transect with which it is associated.

Replacement values must precisely match \code{object} in number of
detections and in their order. \code{xy<-} expects a dataframe of x and y
coordinates for points of detection within a `polygon' or `transect'
detector.

}

\value{
For \code{animalID} and \code{trap} a vector of numeric or character values, one per detection.

For \code{alive} a vector of logical values, one per detection.

For \code{occasion}, a vector of numeric values, one per detection.

For \code{xy}, a dataframe with one row per detection and columns `x' and `y'.

For \code{signalframe} , a dataframe containing signal data and
covariates, one row per detection. The data frame has one row per
detection. See \code{\link{signalmatrix}} for a matrix with one row per
cue and columns for different microphones.

For \code{signal} and \code{noise}, a numeric vector with one element per detection.

If \code{object} has multiple sessions, the result is a list with one
component per session.

}

\seealso{ \code{\link{capthist}}, \code{\link{polyID}}, \code{\link{signalmatrix}} }

\examples{

## `captdata' is a demonstration dataset
animalID(captdata)

temp <- sim.capthist(popn = list(D = 1), make.grid(detector
    = "count"))
cbind(ID = as.numeric(animalID(temp)), occ = occasion(temp),
    trap = trap(temp))

## ovensong dataset has very simple signalframe
head(signalframe(signalCH))

}

\keyword{ manip }

