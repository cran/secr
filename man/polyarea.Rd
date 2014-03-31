\name{polyarea}
\alias{polyarea}

\title{
Area of Polygon(s)
}
\description{
Area of a single closed polygon (simple x-y coordinate input) or of
multiple polygons, possibly with holes.
}
\usage{
polyarea(xy, ha = TRUE)
}

\arguments{
  \item{xy}{dataframe or list with components `x' and `y', or a
    SpatialPolygons or SpatialPolygonsDataFrame object from package \pkg{sp}}
  \item{ha}{logical if TRUE output is converted from square metres to hectares}
}
\details{
  For SpatialPolygons or SpatialPolygonsDataFrame objects, the packages
  \pkg{sp} and \pkg{rgeos} are used.
}

\value{
  A scalar.
}

\seealso{
  \code{\link{buffer.contour}}
}

\examples{
polyarea(make.grid(hollow = TRUE))
}

\keyword{ manip }
