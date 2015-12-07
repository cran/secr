\name{discretize}
\alias{discretize}

\title{Rasterize Area Search or Transect Data}

\description{

  It is sometimes useful to re-cast area search (polygon) data as if it was from a
  set of closely spaced point detectors, i.e. to rasterize the detection locations. This function makes that conversion. Each polygon detector in the input is replaced by a number of point detectors, each representing a square pixel. Detections are mapped to the new detectors on the basis of their x-y coordinates.
  
If \code{object} contains transect data the problem is passed to \code{\link{snip}} and \code{reduce.capthist}.
  
}

\usage{

discretize(object, spacing = 5, outputdetector = c("proximity", "count", "multi"), 
    tol = 0.001, ...)

}

\arguments{
  \item{object}{secr capthist or traps object}
  \item{spacing}{numeric spacing between point detectors in metres }
  \item{outputdetector}{character output detector type}
  \item{tol}{numeric fractional inflation of perimeter (see Details)}
  \item{\dots}{other arguments passed to \code{\link{snip}} if object is transect}
}

\details{
  
  The input should have detector type `polygon' or `polygonX'.

  Centres of a new array of equally spaced detectors are placed within
  each polygon of the input, inflated radially by 1 + tol to avoid some
  inclusion problems. The origin of the superimposed grid is fixed automatically.

  `tol' may be negative, in which case the array(s) will be shrunk relative
  to the polygon(s). 

}

\value{
  
  A capthist or traps object of the requested detector type, but otherwise
  carrying forward all attributes of the input. The embedded traps object has a 
  factor covariate `polyID' recording the polygon to which each point detector relates.
  
}

\note{
  
  Consider the likely number of detectors in the output before you start. 
  
}

\seealso{
  
  \code{\link{reduce.capthist}}
  
}

\examples{

## generate some polygon data
pol <- make.poly()
CH <- sim.capthist(pol, popn=list(D=30), detectpar=list(g0 = 0.3))
plot(CH, border = 10, gridl = FALSE, varycol = FALSE)

## discretize and plot
CH1 <- discretize(CH, spacing = 10, output = 'count')
plot(CH1, add = TRUE, cappar = list(col = 'orange'), varycol =
    FALSE, rad = 0)
plot(traps(CH1), add = TRUE)

## show how detections are snapped to new detectors
newxy <- traps(CH1)[nearesttrap(xy(CH),traps(CH1)),]
segments(xy(CH)[,1], xy(CH)[,2], newxy[,1], newxy[,2])

}

\keyword{ manip }

