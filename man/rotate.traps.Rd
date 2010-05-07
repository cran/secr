\name{rotate.traps}
\alias{rotate.traps}

\title{ Rotate Detectors }

\description{
Rotate detectors while retaining other attributes.
}

\usage{
\method{rotate}{traps}(object, degrees, centrexy = NULL, \dots)
}

\arguments{
  \item{object}{ object of class \code{traps} }
  \item{degrees}{ clockwise angle of rotation in degrees }
  \item{centrexy}{ vector with x,y coordinates of point about which to rotate }
  \item{\dots}{ other arguments (not used). }
}

\value{
An object of class \code{traps} with the location of each detector rotated about the centre.
}

\details{
May be used with \code{flip.traps}, \code{\link{rbind.traps}} and \code{\link{shift.traps}} to create complex geometries.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{traps}}, \code{\link{shift.traps}} }

\examples{
hollow1 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
hollow2 <- make.grid(nx = 8, ny = 8, hollow = TRUE)
nested <- rbind (hollow1, rotate(hollow2, 45, c(70, 70)))
plot(nested, gridlines = FALSE)
}

\keyword{ manip }

