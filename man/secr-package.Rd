\name{secr-package}
\alias{secr-package}
\alias{secr}

\docType{package}

\title{Spatially Explicit Capture--Recapture Models }

\description{ Analyse data from a spatially distributed animal
population sampled with an array of passive detectors, such as traps. }

\details{
\tabular{ll}{
Package: \tab secr\cr
Type: \tab Package\cr
Version: \tab 1.3.0\cr
Date: \tab 2010-03-11\cr
License: \tab GNU General Public License Version 2 or later\cr
}

Warning: Version 1.3.0 should be viewed as a beta release: some functions
may not work with all documented settings. Feedback is very welcome,
including suggestions for additional documentation or new features
consistent with the overall design.

Data comprise the locations of detectors (traps) in an object of class
'traps' and the detection histories of individually marked animals in an
object of class 'capthist'. Models for population density and detection
are defined using symbolic formula notation. Possible predictors for
detection probability include several pre-defined variables (t, b etc.)
corresponding to 'time', 'behaviour' and other effects. Habitat is
distinguished from nonhabitat with an object of class 'mask'. Models are
fitted by maximizing either the full likelihood or the likelihood
conditional on the number of individuals (\eqn{n}). Conditional
likelihood models, while limited to homogeneous Poisson density,
allow continuous individual covariates for detection. Fitting creates an
object of class \code{secr}. Generic methods (plot, print, summary etc.)
are available for each object class. 

A more extensive overview can be got by typing

RShowDoc ('secr-overview', package='secr') 

at the R prompt after the package has been loaded.

The analyses in \pkg{secr} extend those available in the software
Density (see \url{www.otago.ac.nz/density} for the most recent version
of Density).

}

\author{
Murray Efford \email{murray.efford@otago.ac.nz}
}

\section{Acknowledgements}{

David Borchers made these methods possible with his work on the
likelihood, and I'm grateful for his continuing advice. Jeff Laake
provided encouragement and reviewed an early version. Ray Brownrigg got
my Windows code running under Unix. Deanna Dawson editted some of the
documentation (the cleaner bits!) and her support and collaboration were
important throughout.

}

\references{ 

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

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
