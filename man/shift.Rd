\name{shift}
\alias{shift}

\title{ Shift Points }

\description{
Translate an array of points.
}

\usage{
shift (object, shiftxy, ...)
}

\arguments{
  \item{object}{ a 2-column matrix or object that can be coerced to a matrix }
  \item{shiftxy}{ vector of x and y displacements }
  \item{\dots}{ other arguments (not used) }
}

\value{
A matrix with the location of each point shifted by the desired amount.
}

\details{
This is a generic function. The default method is redundant, but the method for \code{traps} objects may be useful.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{rotate}}, \code{\link{flip}} }

\examples{
temp <- matrix(runif (20) * 2 - 1, nc = 2)
temp2 <- shift(temp, c(0.1, 0.1))
plot(temp, xlim=c(-1.5,1.5), ylim = c(-1.5,1.5), pch = 16)
points (0,0, pch=2)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)	
}

\keyword{ manip }

