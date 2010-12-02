\name{polyarea}
\alias{polyarea}

\title{
Area of simple polygon
}
\description{
Area of a closed polygon
}
\usage{
polyarea(xy, ha = TRUE)
}
\arguments{
  \item{xy}{dataframe or list with components `x' and `y'}
  \item{ha}{logical if TRUE output is converted from square metres to hectares}
}
\value{
  A scalar
}
\seealso{
  \code{\link{buffer.contour}}
}

\examples{
polyarea(make.grid(hollow = TRUE))
}

\keyword{ manip }
