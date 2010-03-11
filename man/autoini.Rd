\name{autoini}
\alias{autoini}

\title{ Initial Parameter Values for SECR }

\description{
  Find plausible initial parameter values for \code{\link{secr.fit}}.  A simplified model is fitted by a fast ad hoc method.
}

\usage{
autoini(capthist, mask, detectfn = 0, thin = 0.2)
}

\arguments{

  \item{capthist}{\code{capthist} object}
  \item{mask}{\code{mask} object compatible with the detector layout in \code{capthist}}
  \item{detectfn}{shape of detection function 0 = half-normal}
  \item{thin}{proportion of points to retain in mask }}

\details{ 

Plausible starting values are needed to avoid numerical
problems when fitting SECR models. Actual models
to be fitted will usually have more than the three basic parameters
output by \code{autoini}; other initial values can usually be set to
zero for \code{secr.fit}. If the algorithm encounters problems obtaining
a value for g0, the default value of 0.1 is returned. 

Only the half-normal detection function is available in \code{autoini} (cf
other options in e.g. \code{\link{sim.capthist}}).

\code{autoini} implements a modified version of the algorithm proposed
by Efford et al. (2004). In outline, the algorithm is
\enumerate{ 
\item Find value of sigma that predicts the 2-D dispersion of individual locations (see \code{\link{RPSV}})
\item Find value of g0 that, with sigma, predicts the observed mean number of captures per individual (by algorithm of Efford et al. (2009, Appendix 2))
\item Compute the effective sampling area from g0, sigma, using thinned mask (see \code{\link{esa}})
\item Compute D = \eqn{n}/esa(g0, sigma), where \eqn{n} is the number of individuals detected
}

Here 'find' means solve numerically for zero difference between the observed and predicted values, using \code{\link{uniroot}}.

If \code{RPSV} cannot be computed the algorithm tries to use observed
mean recapture distance \eqn{\bar{d}}{d-bar}. Computation of
\eqn{\bar{d}}{d-bar} fails if there no recaptures, and all returned
values are NA.

A proportion 1--\code{thin} of the points in the mask may be discarded at random to speed execution.
}

\value{
  A list of parameter values :
  \item{ D }{Density (animals per hectare)}
  \item{ g0 }{Magnitude (intercept) of detection function }
  \item{ sigma }{Spatial scale of detection function (m) }
}

\references{ 

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software for analysing capture--recapture data from passive detector arrays. \emph{Animal Biodiversity and Conservation} \bold{27},
 217--228. 

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz}}

\note{

\code{autoini} may in future include an option to use \code{\link{RPSV}} instead of \code{\link{dbar}}.

}

\seealso{ \code{\link{capthist}}, \code{\link{mask}}, \code{\link{secr.fit}}, \code{\link{dbar}}  }

\examples{
demotraps <- make.grid()
demomask <- make.mask(demotraps)
demoCH <- sim.capthist (demotraps, popn = list(D = 5, buffer = 100))
autoini (demoCH, demomask) 
}
\keyword{ models }
