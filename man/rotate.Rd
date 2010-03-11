\name{rotate}
\alias{rotate}
\title{ Rotate Points }

\description{
Rotate a set of points. 
}

\usage{
rotate (object, degrees, centrexy = NULL, \dots)
}

\arguments{
  \item{object}{ object that may be coerced to a numeric matrix }
  \item{degrees}{ clockwise angle of rotation in degrees }
  \item{centrexy}{ vector with xy coordinates of rotation centre }
  \item{\dots}{ other arguments (not used) }
}

\details{
The first column of \code{object} holds the x coordinates of the points and the second holds the y coordinates. If \code{centrexy} is NULL then rotation is about (0,0).

A generic method, introduced for the class-specific method \code{\link{rotate.traps}}.
}

\value{
A matrix with the location of each point rotated about the centre.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{shift}} }

\examples{
temp <- matrix(runif (20) * 2 - 1, nc = 2)
temp2 <- rotate(temp, 25)
plot(temp, xlim=c(-1.5,1.5), ylim = c(-1.5,1.5), pch = 16)
points (0,0, pch=2)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)	
}

\keyword{ manip }

