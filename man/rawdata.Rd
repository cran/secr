\name{rawdata}
\alias{rawdata}
\alias{captXY}
\alias{trapXY}

\docType{data}

\title{ Demonstration Data }

\description{
  Two dataframes containing the demonstration data from program Density 4.1 (Efford 2007).
}

\usage{
data(rawdata)
}

\details{
Fictitious data on 235 captures of 76 individuals over 5 occasions.
100 traps 30 metres apart on a square grid with origin (365 365).

Dataframe \code{trapXY} contains the data from the Density input file
`trap.txt', and \code{captXY} contains the data from `capt.txt'.

}

\source{
  Efford, M. G. (2007) Density 4.1: software for spatially
explicit capture-recapture.  Department of Zoology, University of Otago,
Dunedin, New Zealand. \url{http://www.otago.ac.nz/density}.
}

\seealso{ \code{\link{read.traps}}, \code{\link{make.grid}},
    \code{\link{make.capthist}} 
}

\examples{
## load demonstration data
data(rawdata)  ## trapXY, captXY
## construct a traps object from raw trap data
temptrap <- read.traps(data = trapXY, detector = 'single')
plot(temptrap)
## construct a capthist object
tempcapt <- make.capthist(captXY, temptrap, fmt='XY')
plot(tempcapt, tracks = TRUE)
## fmt='trapID' is usually simpler 

}

\keyword{datasets}
