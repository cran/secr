\name{flip}
\alias{flip}

\title{ Flip Points }

\description{
Flip an array of points about a vertical or horizontal axis.
}

\usage{
flip (object, lr = F, tb = F, ...)
}

\arguments{
  \item{object}{ a 2-column matrix or object that can be coerced to a matrix }
  \item{lr}{ either logical for whether array should be flipped left-right, or numeric value for x-coordinate of axis about which it should be flipped left-right }
  \item{tb}{ either logical for whether array should be flipped top-bottom, or numeric value for y-coordinate of axis about which it should be flipped top-bottom }
  \item{\dots}{ other arguments (not used) }
}

\details{
Logical values for \code{lr} or \code{tb} indicate that points should be flipped about the mean on the relevant axis. 

Numeric values indicate the particular axis value(s) about which points should be flipped.  The default arguments result in no change.

This is a generic function. A method is provided for \code{traps} objects.
}

\value{
A matrix with the coordinates of each point reflected about the desired axis or axes.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{traps}}, \code{\link{rotate.traps}} }

\examples{
temp <- matrix(runif (20) * 2 - 1, nc = 2)
temp2 <- flip(temp, lr = 1)
plot(temp, xlim=c(-1.5,4), ylim = c(-1.5,1.5), pch = 16)
points (temp2, pch = 1)
arrows (temp[,1], temp[,2], temp2[,1], temp2[,2], length = 0.1)	
}

\keyword{ manip }

