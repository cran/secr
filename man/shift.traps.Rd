\name{shift.traps}
\alias{shift.traps}

\title{ Shift Detectors }

\description{
Translate detectors while retaining other attributes.
}

\usage{
\method{shift}{traps}(object, shiftxy, \dots)
}

\arguments{
  \item{object}{ object of class \code{traps} }
  \item{shiftxy}{ vector with displacements in x and y directions }
  \item{\dots}{ other arguments (not used) }
}

\value{
An object of class \code{traps} with the location of each detector shifted by the desired amount.
}

\details{
May be used with \code{\link{rbind.traps}} and \code{\link{rotate.traps}} to create complex geometries.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{traps}}, \code{\link{rotate.traps}}, \code{\link{flip.traps}} }

\examples{
hollow1 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
hollow2 <- shift(make.grid(nx = 6, ny = 6, hollow = TRUE), c(20, 20))
nested <- rbind (hollow1, hollow2)
plot(nested, gridlines = FALSE, lab = TRUE)
}

\keyword{ manip }

