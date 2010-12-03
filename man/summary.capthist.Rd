\name{summary.capthist}
\alias{summary.capthist}
\alias{print.summary.capthist}
\alias{counts}

\title{ Summarise Detections }

\description{
  Concise description of \code{capthist} object.
}
\usage{
\method{summary}{capthist}(object, ...)

\method{print}{summary.capthist}(x, ...)

counts(CHlist, counts = 'M(t+1)')

}

\arguments{
  \item{object }{ \code{\link{capthist}} object }
  \item{x }{ \code{summary.capthist} object }
  \item{\dots}{ arguments passed to other functions }

  \item{CHlist}{ capthist object, especially a multi-session object }
  \item{counts}{ character vector of count names }

}

\value{
An object of class \code{summary.capthist}, a list with at least these components
  \item{detector }{detector type in \{'single', 'multi', 'proximity'\} }
  \item{ndetector }{number of detectors }
  \item{xrange }{range of x coordinates of detectors }
  \item{yrange }{range of y coordinates of detectors }
  \item{spacing }{mean distance from each trap to nearest other trap }
  \item{counts }{matrix of summary counts (rows) by occasion (columns).  See Details.}
  \item{dbar }{mean recapture distance }
  \item{RPSV }{root pooled spatial variance }
}

\details{
These counts are reported by \code{summary.capthist}  
\tabular{ll}{
  n \tab number of individuals detected on each occasion \cr
  u \tab number of individuals detected for the first time on each occasion \cr
  f \tab number of individuals detected exactly f times \cr
  M(t+1) \tab cumulative number of individuals detected \cr
  losses \tab number of individuals reported as not released on each occasion \cr
  detections \tab number of detections, including within-occasion 'recaptures' \cr
  traps visited \tab number of detectors at which at least one detection was recorded \cr
  traps set \tab number of detectors, excluding any `not set' in \code{usage} attribute of \code{traps} attribute \cr
}

\code{counts} may be used to return the specified counts in a compact 
session x occasion table. If more than one count is named then a list is returned
with one component for each type of count.

}

\seealso{ \code{\link{dbar}}, \code{\link{RPSV}}, \code{\link{capthist}} }

\examples{
temptrap <- make.grid(nx = 5, ny = 3)
summary(sim.capthist(temptrap))
summary(sim.capthist(temptrap))$counts['n',]
}

\keyword{ models }

