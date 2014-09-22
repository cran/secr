\name{pdot}
\alias{pdot}

\title{ Net Detection Probability }

\description{

Compute spatially explicit net probability of detection for individual(s) 
at given coordinates. 

}

\usage{
pdot(X, traps, detectfn = 0, detectpar = list(g0 = 0.2,
    sigma = 25, z = 1), noccasions = NULL, binomN = NULL,
    userdist = NULL) 

}

\arguments{
  \item{X}{ vector or 2-column matrix of coordinates }
  \item{traps}{ \code{traps} object }
  \item{detectfn}{ integer code for detection function q.v. }
  \item{detectpar}{ a list giving a value for each named parameter of detection function }
  \item{noccasions}{ number of sampling intervals (occasions) }
  \item{binomN}{ integer code for discrete distribution (see
    \code{\link{secr.fit}}) }
  \item{userdist}{user-defined distance function or matrix (see \link{details})}
}
\details{

  If \code{traps} has a \link{usage} attribute then \code{noccasions} is
  set accordingly; otherwise it must be provided.
  
The probability computed is \eqn{p.(\mathbf{X}) = 1 - \prod\limits _{k}
\{1 - p_s(\mathbf{X},k)\}^{S}}{p.(X) = 1 - (1 - prod(p_s(X,k))^S} where
the product is over the detectors in \code{traps}, excluding any not
used on a particular occasion. The per-occasion detection function
\eqn{p_s}{p_s} is halfnormal (0) by default, and is assumed not to vary
over the \eqn{S} occasions.

For detection functions (10) and (11) the signal threshold `cutval' should be
included in \code{detectpar}, e.g., \code{detectpar = list(beta0 = 103, beta1
= -0.11, sdS = 2, cutval = 52.5)}.

The calculation is not valid for single-catch traps because
\eqn{p.(\mathbf{X})}{p.(X)} is reduced by competition between animals.

\code{userdist} cannot be set if `traps' is any of polygon, polygonX,
transect or transectX.

}

\value{
A vector of probabilities, one for each row in X.
}

\seealso{ \code{\link{secr}}, \code{\link{make.mask}},
  \code{\link{Detection functions}}, \code{\link{pdot.contour}} }

\examples{
  temptrap <- make.grid()
  ## per-session detection probability for an individual centred
  ## at a corner trap. By default, noccasions = 5.
  pdot (c(0,0), temptrap, detectpar = list(g0 = 0.2, sigma = 25),
    noccasions = 5)
}

\keyword{ manip }

