\name{trap.builder}
\alias{trap.builder}
\alias{mash}
\alias{cluster.counts}
\alias{cluster.centres}

\title{
Complex Detector Layouts
}

\description{
  
Construct detector layouts comprising small arrays (clusters) replicated
across space, possibly at a probability sample of points.

}

\usage{

trap.builder (n = 10, cluster, region = NULL, frame =
    NULL, method = 'SRS', edgemethod = 'clip', samplefactor = 2,
    ranks = NULL, rotation = NULL, detector, plt = FALSE,
    add = FALSE) 

mash (object, origin = c(0,0), clustergroup = NULL)

cluster.counts (object)

cluster.centres (object)

}

\arguments{
  \item{n}{ integer number of clusters (ignored if method = "all") }
  \item{cluster}{ traps object }
  \item{region}{ bounding polygon}
  \item{frame}{ data frame of points used as a finite sampling frame }
  \item{method}{ character string (see Details) }
  \item{edgemethod}{character string (see Details) }
  \item{samplefactor}{ oversampling to allow for rejection of edge clusters
    (multiple of n) }
  \item{ranks}{ vector of relative importance (see Details)}
  \item{rotation}{ angular rotation of each cluster about centre (degrees) }
  \item{detector}{ character detector type (see \code{\link{detector}}) }
  \item{plt}{ logical: should array be plotted? }
  \item{add}{ logical: add to existing plot }
  \item{object}{ single-session multi-cluster capthist object, or traps
    object for \code{cluster.centres} }
  \item{origin}{ new coordinate origin for detector array }
  \item{clustergroup}{ list of vectors subscripting the clusters to be
  mashed }
}

\details{

  The detector array in \code{cluster} is replicated \code{n}
  times and translated to centres sampled from the area sampling frame
  in \code{region} or the finite sampling frame in \code{frame}. Each
  cluster may be rotated about its centre either by a fixed number of
  degrees (\code{rotation} positive), or by a random angle (\code{rotation}
  negative).
  
  If the \code{cluster} argument is not provided then single detectors of
  the given type are placed according to the design.

  The sampling frame is finite (the points in \code{frame}) whenever
  \code{frame} is not NULL. If \code{region} and \code{frame} are both
  specified, sampling uses the finite frame but sites may be clipped
  using the polygon.

  \code{region} may be a two-column matrix or dataframe of x-y coordinates
  for the boundary, or a SpatialPolygonsDataFrame object from \pkg{sp}.
  
  \code{method} may be `SRS', `GRTS', `all' or `rank'. `SRS' takes a simple
  random sample (without replacement in the case of a finite sampling
  frame). `GRTS' takes a spatially representative sample using the
  `generalized random tessellation stratified' (GRTS) method of Stevens
  and Olsen (2004). `all' replicates \code{cluster} across all points in
  the finite sampling frame. `rank' selects \code{n} sites from
  \code{frame} on the basis of their ranking on the vector `ranks',
  which should have length equal to the number of rows in
  \code{frame}; ties are resolved by drawing a site at random.
  
  \code{edgemethod} may be `clip' (reject individual detectors),
  `allowoverlap' (no action) or `allinside' (reject whole cluster if any
  component is outside \code{region}). Sufficient additional samples
  (\code{(samplefactor--1) * n}) must be drawn to allow for replacement
  of any rejected clusters; otherwise, an error is reported (`not enough
  clusters within polygon').
  
  The package \pkg{sp} is required. GRTS samples require function
  \code{grts} in package \pkg{spsurvey} of Olsen and Kincaid. Much more
  sophisticated sampling designs may be specified by using \code{grts}
  directly.

  \code{mash} collapses a multi-cluster capthist object as if all
  detections were made on a single cluster. The new detector coordinates
  in the `traps' attribute are for a single cluster with (min(x),
  min(y)) given by \code{origin}. \code{clustergroup} optionally selects
  one or more groups of clusters to mash; if \code{length(clustergroup)
  > 1} then a multisession capthist object will be generated, one
`session' per clustergroup. By default, all clusters are mashed.

  \code{mash} discards detector-level covariates and occasion-specific
  `usage', with a warning.

  \code{cluster.counts} returns the number of \emph{distinct}
  individuals detected per cluster in a single-session multi-cluster
  capthist object.

}

\value{

  \code{trap.builder} produces an object of class `traps'. A covariate
  `cluster' retains the sequence number of the cluster from which a
  location is derived. A covariate `trapnum' retains the within-cluster
  sequence number of the detector.

  \code{method = 'GRTS'} causes messages to be displayed regarding the
    stratum (always `None'), and the initial, current and final number
    of levels from the GRTS algorithm.

  \code{plt = TRUE} causes a plot to be displayed, including the polygon
  or finite sampling frame as appropriate.

  \code{mash} produces a capthist object with the same number of rows as
  the input but different detector numbering and `traps'. An attribute
  `n.mash' is a vector of the numbers recorded at each cluster; its
  length is the number of clusters. An attribute `centres' is a a
  dataframe containing the x-y coordinates of the cluster centres. The
  \code{predict} method for secr objects and the function \code{derived}
  both recognise and adjust for mashing.

  \code{cluster.counts} returns a vector with the number of individuals
  detected at each cluster.

  \code{cluster.centres} returns a dataframe of x- and y-coordinates.

}

\references{
  
  Stevens, D. L., Jr., and Olsen, A. R. (2004) Spatially-balanced
  sampling of natural resources. \emph{Journal of the American
  Statistical Association} \bold{99}, 262--278.

}

\note{

  The function \code{\link{make.systematic}} should be used to generate
  systematic random layouts.

}

\seealso{
\code{\link{make.grid}}, \code{\link{traps}},
\code{\link{make.systematic}},
\code{\link{clusterID}}
}

\examples{

## solitary detectors placed randomly within a rectangle
tempgrid <- trap.builder (n = 10, method = 'SRS',
    region = cbind(x = c(0,1000,1000,0),
    y = c(0,0,1000,1000)), plt = TRUE)

## GRTS sample of mini-grids within a rectangle
## edgemethod = 'allinside' avoids truncation at edge
minigrid <- make.grid(nx = 3, ny = 3, spacing = 50,
    detector = 'proximity')
tempgrid <- trap.builder (n = 20, cluster = minigrid,
    method = 'GRTS', edgemethod = 'allinside', region =
    cbind(x = c(0,6000,6000,0), y = c(0,0,6000,6000)),
    plt = TRUE)

## one detector in each 100-m grid cell -
## a form of stratified simple random sample
origins <- expand.grid(x = seq(0, 900, 100),
    y = seq(0, 1100, 100))
XY <- origins + runif(10 * 12 * 2) * 100
temp <- trap.builder (frame = XY, method = 'all',
    detector = 'multi')
plot(temp, border = 0)  ## default grid is 100 m

## regular lattice of mini-arrays
tempgrid <- trap.builder (cluster = minigrid , method =
    'all', frame = expand.grid(x = seq(1000, 5000, 2000),
    y = seq(1000, 5000, 2000)), plt = TRUE)

## simulate some data
tempcapt <- sim.capthist(tempgrid, popn = list(D=10))
cluster.counts(tempcapt)
cluster.centres(tempgrid)

## 'mash' the CH
summary(mash(tempcapt))

## compare timings (estimates are near identical)
\dontrun{
tempmask1 <- make.mask(tempgrid, type = 'clusterrect',
    buffer = 200, spacing = 10)
secr.fit(tempcapt, mask = tempmask1)         ## 241.52 s

tempmask2 <- make.mask(minigrid, spacing = 10)
secr.fit(mash(tempcapt), mask = tempmask2)   ## 2.85 s
## density estimate is adjusted automatically
## for the number of mashed clusters (9)
}

## two-phase design: preliminary sample across region,
## followed by selection of sites for intensive grids
\dontrun{
arena <- data.frame(x = c(0,2000,2000,0), y = c(0,0,2500,2500))
t1 <- make.grid(nx = 1, ny = 1)
t4 <- make.grid(nx = 4, ny = 4, spacing = 50)
singletraps <- make.systematic (n = c(8,10), cluster = t1,
    region = arena)
CH <- sim.capthist(singletraps, popn = list(D = 2))
plot(CH, type='n.per.cluster', title='Number per cluster')
temp <- trap.builder(10, frame = traps(CH), cluster = t4,
    ranks = cluster.counts(CH), method = 'rank',
    edgemethod = 'allowoverlap', plt = TRUE, add = TRUE)
}
}

\keyword{ manip }