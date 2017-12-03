\name{make.systematic}
\alias{make.systematic}

\title{
Construct Systematic Detector Design
}

\description{
A rectangular grid of clusters within a polygonal region.
}

\usage{

make.systematic(n, cluster, region, spacing = NULL, origin = NULL, ...)

}


\arguments{
  
  \item{n}{integer approximate number of clusters (see Details) }
  
  \item{cluster}{traps object defining a single cluster }
  
  \item{region}{dataframe or SpatialPolygonsDataFrame with coordinates
    of perimeter }

  \item{spacing}{ scalar distance between cluster centres }

  \item{origin}{vector giving x- and y-cooordinates of fixed grid origin
    (origin is otherwise random)}
   
  \item{\dots}{ arguments passed to \code{trap.builder} }
}

\details{

  \code{region} may be any shape. The \pkg{sp} class
  SpatialPolygonsDataFrame is useful for complex shapes and input from
  shapefiles using \pkg{rgdal} (see Examples). Otherwise,
  \code{region} should be a dataframe with columns `x' and `y'.

  \code{spacing} may be a vector with separate values for spacing in x-
  and y- directions. If \code{spacing} is provided then \code{n} is ignored.
  
  If \code{n} is a scalar, the spacing of clusters is determined from
  the area of the bounding box of \code{region} divided by the requested
  number of clusters (this does not necessarily result in exactly n
  clusters). If \code{n} is a vector of two integers these are taken to be
  the number of columns and the number of rows.

  After preparing a frame of cluster centres, \code{make.systematic}
  calls \code{\link{trap.builder}} with method = `all'; \dots allows the
  arguments `rotation', `edgemethod', `plt', and `detector' to be
  passed. Setting the \code{trap.builder} arguments \code{frame},
  \code{method}, and \code{samplefactor} has no effect.
  
}

\value{
A single-session `traps' object.
}

\seealso{
  \code{\link{trap.builder}}, \code{\link{cluster.centres}},
  \code{\link[rgdal]{readOGR}}
}

\note{

    Do not confuse with the simpler function \code{\link{make.grid}},
    which places single detectors in a rectangular array.

}

\examples{

mini <- make.grid(nx = 2, ny = 2, spacing = 100)
region <- cbind(x=c(0,2000,2000,0), y=c(0,0,2000,2000))
temp <- make.systematic(25, mini, region, plt = TRUE)
temp <- make.systematic(c(6, 6), mini, region, plt = TRUE,
    rotation = -1)

## Example using shapefile "possumarea.shp" in
## "extdata" folder. By default, each cluster is 
## a single multi-catch detector

\dontrun{

datadir <- system.file("extdata", package = "secr")
possumarea <- rgdal::readOGR(dsn = datadir, layer = "possumarea")

possumgrid <- make.systematic(spacing = 100, region =
    possumarea, plt = TRUE)

## or with 2 x 2 clusters
possumgrid2 <- make.systematic(spacing = 300,
    cluster = make.grid(nx = 2, ny = 2, spacing = 100),
    region = possumarea, plt = TRUE, edgemethod =
    "allinside")
## label clusters
text(cluster.centres(possumgrid2), levels(clusterID
    (possumgrid2)), cex=0.7)

## If you have GPSBabel installed and on the Path
## then coordinates can be projected and uploaded
## to a GPS with `writeGPS', which also requires the
## package `proj4'. Defaults are for a Garmin GPS
## connected by USB.

writeGPS(possumgrid, proj = "+proj=nzmg")

}

}

\keyword{ manip }

