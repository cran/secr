\name{capthist.parts}
\alias{animalID}
\alias{occasion}
\alias{trap}
\alias{xy}
\alias{signal}
\alias{xy<-}
\alias{signal<-}
%% \alias{times}
%% \alias{times<-}

\title{ Dissect Spatial Capture History Object }

\description{
  Extract parts of an object of class 'capthist'.
}

\usage{

animalID(object, names = TRUE)
occasion(object)
trap(object, names = TRUE)
xy(object)
signal(object)
xy(object) <- value
signal(object) <- value
%% times(object)
%% times(object) <- value

}

\arguments{

  \item{object}{ a 'capthist' object }
  \item{names}{ if FALSE the values returned are numeric indices rather than names }
  \item{value}{ replacement value (see Details) }

}
\details{

These functions extract data on detections, ignoring occasions when an
animal was not detected. 

\code{trap} returns polygon or transect numbers if \code{traps(object)}
has detector type 'polygon' or 'transect'.

Replacement values must precisely match \code{object} in number of
detections in their order. \code{xy<-} expects a dataframe of x and y
coordinates for points of detection within a 'polygon' or 'transect'
detector. 

}

\value{
For \code{animalID} and \code{trap} a vector of numeric or character values, one per detection.

For \code{occasion}, a vector of numeric values, one per detection.

For \code{xy}, a dataframe with one row per detection and columns 'x' and 'y'.

%%For \code{signal} and \code{times}, a numeric vector with one element per detection.
For \code{signal}, a numeric vector with one element per detection.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{capthist}}, \code{\link{polyID}} }

\examples{

data(captdata)
animalID(captdata)

temp <- sim.capthist(popn=list(D=1), make.grid(detector='count'))
cbind(ID=as.numeric(animalID(temp)), occ=occasion(temp), trap=trap(temp))
 
}

\keyword{ models }

