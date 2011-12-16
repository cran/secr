\name{pointsInPolygon}
\alias{pointsInPolygon}
\title{Points Inside Polygon}

\description{
Determines which of a set of  points lie inside a closed polygon or at
least one of a set of polygons
}

\usage{
pointsInPolygon(xy, poly)
}

\arguments{
  
  \item{xy}{2-column matrix or dataframe of x-y coordinates for points to assess}

  \item{poly}{2-column matrix or dataframe containing perimeter points of polygon, or
    a SpatialPolygonsDataFrame object from package \pkg{sp}}
  
}

\details{

  If \code{poly} is a SpatialPolygonsDataFrame object then the function
\code{overlay} is used from \pkg{sp}. This allows multiple polygons and
polygons with holes.

Otherwise, the algorithm is adapted from some code posted on the S-news list
by Peter Perkins (23/7/1996). The polygon should be closed (last point
same as first).

}

\value{
  Vector of logical values, one for each row in \code{xy}
}

\seealso{
\code{\link[sp]{overlay}} 
}

\examples{

## 100 random points in unit square
xy <- matrix(runif(200), ncol = 2)
## triangle centred on (0.5, 0.5)
poly <- data.frame(x = c(0.2,0.5,0.8,0.2), y = c(0.2,0.8,0.2,0.2))
plot(xy, pch = 1 + pointsInPolygon(xy, poly))
lines(poly)

}

\keyword{ manip }