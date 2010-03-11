\name{write.capthist}
\alias{write.capthist}
\alias{write.traps}

\title{ Write Data to Text File }
\description{
  Export detections or detector layout to a text file in form suitable for DENSITY
}

\usage{

write.capthist(object, file = "", ..., deblank = TRUE, header = TRUE,
    sess = '1', ndec = 2)

write.traps(object, file = "", ..., deblank = TRUE,  header = TRUE, 
    ndec = 2)

}
\arguments{
  \item{object}{ \code{capthist} or \code{traps} object }
  \item{file}{ name of output file (character string) }
  \item{\dots}{ other arguments passed to \code{write.table} }
  \item{deblank}{ remove blanks from character string used to identify detectors? }
  \item{header}{ output descriptive header? }
  \item{sess}{ session identifier }
  \item{ndec}{ number of digits after decimal point for x,y coordinates }
}
\details{
  Existing file will be replaced without warning.
}
\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\examples{
  data (captdata)
  write.capthist (captdata)
}
\keyword{ IO }
