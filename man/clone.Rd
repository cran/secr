\name{clone}
\alias{clone}
\alias{clone.default}
\alias{clone.popn}

\title{
Replicate Rows
}

\description{
  
  Clone rows of a 2-D object a constant or random number of times

}

\usage{
  \method{clone}{default}(object, type, ...)
  \method{clone}{popn}(object, type, ...)
}

\arguments{
  \item{object}{any object}
  \item{type}{character 'constant', 'poisson' or 'nbinom'}
  \item{\dots}{other arguments for distribution function }
}
\details{
  
The \dots argument specifies the number of times each row should be
repeated. For random distributions (Poisson or negative binomial) \dots
provides the required parameter values: \code{lambda} for Poisson,
\code{size, prob} or \code{size, mu} for negative binomial.

One application is to derive a population of cues from a popn object,
where each animal in the original popn generates a number of cues from
the same point.

}

\value{

Object of same class as \code{object} but with varying number of
rows. For \code{clone.popn} an attribute `freq' is set, a vector of
length equal to the original number of rows giving the number of repeats
(including zeros).

If \code{popn} is a multi-session population the returned value will be
a multi-session population of the same length.

}

\seealso{
\code{\link{sim.popn}}
}

\examples{
## population of animals at 1 / hectare generates random
## Poisson number of cues, lambda = 5
mics4 <- make.grid( nx = 2, ny = 2, spacing = 44, detector = "signal")
pop <- sim.popn (D = 1, core = mics4, buffer = 300, nsession = 66)
pop <- clone (pop, "poisson", 5)
attr(pop[[1]],"freq")
}

\keyword{ manip }
