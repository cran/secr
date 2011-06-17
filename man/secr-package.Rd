\name{secr-package}
\alias{secr-package}
\alias{secr}

\docType{package}

\title{Spatially Explicit Capture--Recapture Models }

\description{
  Functions to estimate the density of a spatially distributed animal
population sampled with an array of passive detectors, such as traps, or
by searching polygons or transects. }

\details{
\tabular{ll}{
Package: \tab secr\cr
Type: \tab Package\cr
Version: \tab 2.1.0\cr
Date: \tab 2011-06-17\cr
License: \tab GNU General Public License Version 2 or later\cr
}

Data comprise the locations of detectors (traps, searched areas, etc.)
in an object of class `traps', and the detection histories of
individually marked animals. Individual histories are stored in an
object of class `capthist' that includes the relevant `traps'
object. Models for population density (animals per hectare) and
detection are defined using symbolic formula notation. Density models
may include spatial or temporal trend. Possible predictors for detection
probability include both pre-defined variables (t, b, etc.)
corresponding to `time', `behaviour' and other effects), and
user-defined covariates of several kinds. Habitat is distinguished from
nonhabitat with an object of class `mask'.

Models are fitted by maximizing either the full likelihood or the
likelihood conditional on the number of individuals observed
(\eqn{n}). Conditional likelihood models are limited to homogeneous
Poisson density, but allow continuous individual covariates for
detection. Fitting (function \code{\link{secr.fit}}) creates an object
of class \code{secr}. Generic methods (plot, print, summary, etc.) are
provided for each object class.

A link at the bottom of this and every help page takes you to the
help index. Several vignettes complement the help pages:

\tabular{ll}{
  \url{../doc/secr-overview.pdf} \tab general introduction \cr
  \url{../doc/secr-datainput.pdf} \tab data formats and input functions \cr
  \url{../doc/secr-polygondetectors.pdf} \tab using polygon and transect detector types \cr
  \url{../doc/secr-sound.pdf} \tab analysing data from microphone arrays\cr
  \url{../doc/secr-finitemixtures.pdf} \tab mixture models for individual heterogeneity \cr
}

The analyses in \pkg{secr} extend those available in the software
Density (see
\href{http://www.otago.ac.nz/density/index.html}{www.otago.ac.nz/density}
for the most recent version of Density).  Help is available on the
`DENSITY | secr' forum at
\href{http://www.phidot.org/forum/index.php}{www.phidot.org}.  Feedback on the
software is also welcome, including suggestions for additional
documentation or new features consistent with the overall design.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\section{Acknowledgements}{

David Borchers made these methods possible with his work on the
likelihood, and I'm grateful for his continuing advice. Jeff Laake
provided encouragement and reviewed an early version. Ray Brownrigg got
my Windows code running under Unix. Deanna Dawson edited some of the
documentation (the cleaner bits!) and her support and collaboration were
important throughout. Tiago Marques suggested many improvements to the
help pages.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area or transect searches. Unpublished
manuscript.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27},
 217--228.

}

\keyword{ package }

\seealso{
\code{\link{read.capthist}},
\code{\link[secr]{secr.fit}},
\code{\link[secr]{traps}},
\code{\link[secr]{capthist}},
\code{\link[secr]{mask}}
}
\examples{
\dontrun{

## generate some data & plot
detectors  <- make.grid (nx = 10, ny = 10, spacing = 20,
    detector = 'multi')
plot(detectors, label = TRUE, border = 0, gridspace = 20)
detections <- sim.capthist (detectors, noccasions = 5,
    popn = list(D = 5, buffer = 100),
    detectpar = list(g0 = 0.2, sigma = 25))
session(detections) <- 'Simulated data'
plot(detections, border = 20, tracks = TRUE, varycol = TRUE)

## generate habitat mask
mask <- make.mask (detectors, buffer = 100, nx = 48)

## fit model and display results
secr.model <- secr.fit (detections, model = g0~b, mask = mask)
secr.model

}
}
