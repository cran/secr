\name{read.telemetry}
\alias{read.telemetry}

\title{ Import Radio Fixes }

\description{

  A shortcut function for constructing a telemetry capthist object from a file of
telemetry fixes. Telemetry data are generally similar in format to polygon data
(see \code{\link{addTelemetry}}).

}

\usage{
read.telemetry(file = NULL, data = NULL, noccasions = NULL, covnames = NULL,
verify = TRUE, ...)
}

\arguments{
  \item{file}{character name of text file }
  \item{data}{data.frame containing fixes }
  \item{noccasions}{number of occasions (optional) }
  \item{covnames}{character vector of names for individual covariates }
  \item{verify}{logical for whether to check input }
  \item{\dots}{other arguments passed to countfields, read.table etc. }
}

\details{

Input data may be in a text file (argument \code{file}) or a dataframe
(argument \code{data}). Data should be in the XY format for that
function i.e. the first 5 columns should be Session, ID, Occasion, X,
Y. Further columns are treated as individual covariates.

No `traps' input is required. A traps object is generated automatically
from the convex hull of the fixes, inflated by a tiny fraction (1e-8) in
all directions to ensure all fixes lie inside.

}

\value{
  
An secr capthist object including attribute `detectedXY' with the x-y
coordinates, and a `traps' object with detector type = `telemetry'

}

\seealso{
\code{\link{addTelemetry}}, \code{\link{read.capthist}}
}

\examples{

\dontrun{
gps2008CH <- read.telemetry('gps2008.txt')
plot( gps2008CH, gridsp = 10000)
head(gps2008CH)
secr.fit(gps2008CH, start = c(log(4000)))
}

}

\keyword{ IO }
