\name{captdata}
\alias{captdata}
\docType{data}

\title{ Detection Dataset }

\description{

An object of class \code{\link{capthist}} containing the demonstration
data from program Density 4.1 (Efford 2007), specifically the trap
layout data in 'trap.txt' and the captures in 'capt.txt'.

}

\usage{data(captdata)}

\details{

These data are entirely fictitious. They comprise 235 captures of 76
animals in 100 single-catch traps over 5 occasions.

}

\source{

Efford, M.G. (2007) Density 4.1: software for spatially explicit
capture-recapture. Department of Zoology, University of Otago, Dunedin,
New Zealand. \url{http://www.otago.ac.nz/density}.

}

\examples{
data(captdata)
plot(captdata)
secr.fit(captdata)
}

\keyword{datasets}
