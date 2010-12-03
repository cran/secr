\name{traps.info}
\alias{searcharea}
% \alias{searcharea<-}
\alias{polyID}
\alias{polyID<-}
\alias{transectID}
\alias{transectID<-}
\alias{transectlength}

\title{ Detector Attributes }

\description{
  Extract or replace attributes of an object of class 'traps'.
}

\usage{
polyID(object)
polyID(object) <- value
transectID(object)
transectID(object) <- value
searcharea(object)
% searcharea(object) <- value
transectlength(object)
}

\arguments{
  \item{object}{ a 'traps' object }
  \item{value}{ replacement value (see Details) }
}

\details{

The 'polyID' and 'transectID' functions assign and extract the attribute of a 'traps'
object that relates vertices (rows) to particular polygons or transects. The
replacement value should be a factor of length equal to nrow(object).

% The 'searcharea' attribute of a 'quadratbinary' or 'quadratcount' traps
% object is the area in hectares searched at each detector point
% (quadrat). Usually, this is the area of a rectangular pixel determined
% by the detector spacing (spacex, spacey). Replacement creates square
% pixels with dimensions spacex = spacey = value^0.5 * 100. 

The 'searcharea' of a 'polygon' traps object is a vector of the areas of
the component polygons.
%This is a read-only value (i.e. 'searcharea<-' does not apply).

The 'transectlength' of a 'transect' traps object is a vector of the
lengths of the component transects in metres. This is a read-only value.

}

\value{
\code{polyID} - a factor with one level per polygon.
\code{searcharea} - numeric value of polygon areas, in hectares.
\code{transectlength} - numeric value of transect lengths, in metres.
}

\seealso{ 
  \code{\link{traps}}
}

\examples{

## default is a single polygon
temp <- make.grid(detector = 'polygon', hollow = TRUE)
polyID(temp)
plot(temp)

## split in two
temp <- make.grid(detector = 'polygon', hollow = TRUE)
polyID(temp) <- factor(rep(c(1,2),rep(10,2)))
plot(temp)
 
}

\keyword{ manip }

