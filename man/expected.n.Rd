\name{expected.n}
\alias{expected.n}

\title{
Expected Number of Individuals
}

\description{

  Computes the expected number of individuals detected at each detector
or at each cluster of detectors.

}

\usage{

expected.n(object, session = NULL, group = NULL, bycluster
    = TRUE, splitmask = FALSE)

}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  \item{session}{ character session vector }
  \item{group}{ group -- for future use }
  \item{bycluster}{ logical to output the expected number for
  clusters of detectors rather than separate detectors }
  \item{splitmask}{logical for computation method (see Details)}
}

\details{

  The expected number of individuals detected is \eqn{E(n) = \int p.(X)
  D(X) dX}{E(n) = integral p.(X) D(X) dX} where the integration is a
  summation over \code{object$mask}. \eqn{p.(X)} is the probability an
  individual at \eqn{X} will be detected at least once either on the
  whole detector layout (\code{bycluster = FALSE}) or on the detectors
  in a single cluster (see \link{pdot} for more on \eqn{p.}). \eqn{D(X)}
  is the expected density at \eqn{X}, given the model. \eqn{D(X)} is
  constant (i.e. density surface flat) if \code{object$CL == TRUE} or
  \code{object$model$D == ~1}, and for some other possible models.

If the \code{bycluster} option is selected and detectors are not,
in fact, clustered then each detector will be treated as a cluster, with a
warning.

By default, a full habitat mask is used for each detector or
cluster. This is the more robust option. Alternatively, the mask may be
split into subregions defined by the cells closest to each detector or
cluster.

The calculation takes account of any fitted continuous model for spatial
variation in density (note warning).

}

\section{Warning}{
The prediction of density at present considers only the base level of
density covariates, such as pixel-specific habitat variables.
}

\value{

Vector of expected count with length equal to the number of detectors or
the number of clusters. For multi-session data, a list of such vectors.

}

\seealso{
\code{\link{region.N}}
}

\examples{

expected.n(secrdemo.0)
expected.n(ovenbird.model.D)

\dontrun{
mini <- make.grid(nx = 3, ny = 3, spacing = 50, detector =
    'proximity')
tempgrids <- trap.builder (cluster = mini , method = 'all',
    frame = expand.grid(x = seq(1000, 9000, 2000),
    y = seq(1000, 9000, 2000)), plt = TRUE)
capt <- sim.capthist(tempgrids, popn = list(D = 2))
tempmask <- make.mask(tempgrids, buffer = 100,
    type = "clusterbuffer")
fit <- secr.fit(capt, mask = tempmask, trace = FALSE)
En <- expected.n(fit)

## GoF or overdispersion statistic
p <- length(fit$fit$par)
y <- cluster.counts(capt)
## scaled by n-p
sum((y - En)^2 / En) / (length(En)-p)
sum((y - En)^2 / En) / sum(y/En)

}

}

\keyword{ model }
