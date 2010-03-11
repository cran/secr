\name{flip.traps}
\alias{flip.traps}

\title{ Flip Detector Array }

\description{
Flip a detector array about a vertical or horizontal axis.
}

\usage{
\method{flip}{traps} (object, lr = F, tb = F, ...)
}

\arguments{
  \item{object}{ a 2-column matrix or object that can be coerced to a matrix }
  \item{lr}{ either logical for whether array should be flipped left-right, or numeric value for x-coordinate of axis about which it should be flipped left-right }
  \item{tb}{ either logical for whether array should be flipped top-bottom, or numeric value for y-coordinate of axis about which it should be flipped top-bottom }
  \item{\dots}{ other arguments (not used) }
}

\details{
Logical values for \code{lr} or \code{tb} indicate that points should be flipped about the mean on the relevant axis. 

Numeric values indicate the particular axis value(s) about which points should be flipped. 

The default arguments result in no change.
}

\value{
Object of class \code{traps} with the coordinates of each point reflected about the desired axis or axes.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{traps}}, \code{\link{rotate.traps}}, \code{\link{shift.traps}} }

\examples{
par(mfrow=c(1,2), xpd = TRUE)
traps1 <- make.grid(nx = 8, ny = 6, ID = 'numxb')
traps2 <- flip (traps1, lr = TRUE)
plot(traps1, border = 5, lab = TRUE, offset = 7, gridl = FALSE)
plot(traps2, border = 5, lab = TRUE, offset = 7, gridl = FALSE)
}

\keyword{ manip }

