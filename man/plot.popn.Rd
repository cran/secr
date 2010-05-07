\name{plot.popn}
\alias{plot.popn}

\title{ Plot popn Object }

\description{
Display animal locations from a \code{popn} object.
}

\usage{
\method{plot}{popn}(x, add = FALSE, frame = TRUE, ...)
}

\arguments{
  \item{x}{ object of class \code{popn} }
  \item{add}{ logical to add points to an existing plot }
  \item{frame}{ logical to add frame within which points were simulated }
  \item{\dots}{ arguments passed to \code{eqscplot} and \code{points} }
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{popn}}, \code{\link{sim.popn}} }

\examples{

temppopn <- sim.popn(D = 5, expand.grid(x = c(0,100), y = c(0,100)))
plot(temppopn, pch = 16, col = 'blue')

}
\keyword{ hplot }

