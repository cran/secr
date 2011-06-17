\name{plot.popn}
\alias{plot.popn}

\title{ Plot popn Object }

\description{
Display animal locations from a \code{popn} object.
}

\usage{
\method{plot}{popn}(x, add = FALSE, frame = TRUE, circles = NULL, ...)
}

\arguments{
  \item{x}{ object of class \code{popn} }
  \item{add}{ logical to add points to an existing plot }
  \item{frame}{ logical to add frame within which points were simulated
 }
  \item{circles}{ vector giving the radii if circles are to be plotted }
  \item{\dots}{ arguments passed to \code{eqscplot} and \code{points} or
  \code{symbols}}
}

\details{
  If \code{circles} is provided then a circle of the given radius is
  plotted for each animal using the \code{symbols} function. The
  arguments \code{fg} and \code{bg} may be used to control the colour of
  the perimeter and the fill of each circle (see Examples).
}

\seealso{ \code{\link{popn}}, \code{\link{sim.popn}} }

\examples{

temppopn <- sim.popn(D = 5, expand.grid(
    x = c(0,100), y = c(0,100)))
plot(temppopn, pch = 16, col = 'blue')

plot(temppopn, circles = 20, bg = 'tan', fg =
    'white')
plot(temppopn, pch = 16, cex = 0.5, add = TRUE)

}
\keyword{ hplot }

