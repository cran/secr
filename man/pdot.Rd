\name{pdot}
\alias{pdot}

\title{ Net Detection Probability }

\description{

Compute spatially explicit net probability of detection for individual(s) 
at given coordinates. 

}

\usage{
pdot(X, traps, detectfn = 0, detectpar = list(g0 = 0.2, sigma = 20),
    noccasions = 5) 

}

\arguments{
  \item{X}{ coordinates }
  \item{traps}{ \code{traps} object }
  \item{detectfn}{ integer code for detection function q.v. }
  \item{detectpar}{ a list giving a value for each named parameter of detection function }
  \item{noccasions}{ number of intervals (occasions) }
}
\details{

The probability computed is \eqn{p.(\mathbf{X}) = 1 - \prod\limits _{k}
\{1 - p_s(\mathbf{X},k)\}^{S}}{p.(X) = 1 - (1 - prod(p_s(X,k))^S} where
the product is over the detectors in \code{traps}. The
per-occasion detection function \eqn{p_s}{p\_s} is by default
half-normal, and is assumed not to vary over the \eqn{S} occasions.

The calculation is not valid for single-catch traps because
\eqn{p.(\mathbf{X})}{p.(X)} is reduced by competition between animals.

}

\value{
A vector of probabilities, one for each row in X.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{secr}}, \code{\link{make.mask}}, \code{\link{detection functions}} }

\examples{
  temptrap <- make.grid()
  ## per-session detection probability for an individual centred
  ## at a corner trap. By default, noccasions = 5.
  pdot (c(0,0), temptrap, detectpar = list(g0 = 0.2, sigma = 20))
}

\keyword{ manip }

