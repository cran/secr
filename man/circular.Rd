\name{circular}
\alias{circular.r}
\alias{circular.p}

\title{
Circular Probability
}

\description{
  
Functions to answer the question "what radius is expected to include
proportion p of points from a circular bivariate distribution
corresponding to a given detection function", and the reverse. These
functions may be used to relate the scale parameter(s) of a detection
function (e.g., \eqn{\sigma}) to home-range area (specifically, the area
within an activity contour for the corresponding simple home-range
model) (see Note).

}

\usage{

circular.r (p = 0.95, detectfn = 0, sigma = 1, detectpar = NULL) 

circular.p (r = 1, detectfn = 0, sigma = 1, detectpar = NULL) 

}

\arguments{
  \item{p}{vector of probability levels for which radius is required}
  \item{r}{vector of radii for which probability level is required}
  \item{detectfn}{integer code or character string for shape of
      detection function 0 = halfnormal, 1 = hazard rate etc. -- see
      \link{detectfn}}
  \item{sigma}{spatial scale parameter of detection function }
  \item{detectpar}{named list of detection function parameters }
}

\details{

\code{circular.r} is the quantile function of the specified circular
bivariate distribution (analogous to \code{qnorm}, for example). The
quantity calculated by \code{circular.r} is sometimes called `circular
error probable' (see Note).
  
For detection functions with two parameters (intercept and scale) it is
enough to provide \code{sigma}. Otherwise, \code{detectpar} should be a
named list including parameter values for the requested detection
function (g0 may be omitted, and order does not matter).

}

\value{
Vector of values for the required radii or probabilities.
}

\seealso{
  \code{\link{detectfn}},  \code{\link{detectfnplot}}
}

\references{ 
Calhoun, J. B. and Casby, J. U. (1958) Calculation of home range and
density of small mammals. Public Health Monograph No. 55. United States
Government Printing Office.
}

\note{

  The term `circular error probable' has a military origin. It is
commonly used for GPS accuracy with the default probability level set to
0.5 (i.e. half of locations are further than CEP from the true
location). A circular bivariate normal distriubution is commonly assumed
for the circular error probable; this is equivalent to setting
\code{detectfn = "halfnormal"}.

Closed-form expressions are used for the normal and uniform cases; in
the bivariate normal case, the relationship is \eqn{r =
\sigma \sqrt{-2\mbox{ln}(1-p)} }{r = (-2*log(1-p))^0.5 * sigma}. Otherwise,
the probability is computed numerically by integrating the radial
distribution. Numerical integration is not foolproof, so check
suspicious or extreme values.

When \code{circular.r} is used with the default \code{sigma = 1}, the result
may be interpreted as the factor by which sigma needs to be inflated to
include the desired proportion of activity (e.g., 2.45 sigma for 95\%
of points from a circular bivariate normal distribution).

}

\examples{

## Calhoun and Casby (1958) p 3.
## give p = 0.3940, 0.8645, 0.9888
circular.p(1:3)

## halfnormal, hazard-rate and exponential
circular.r ()
circular.r (detectfn = 1, detectpar = list(sigma = 1, z = 4))
circular.r (detectfn = 2)

plot(seq(0, 5, 0.01), circular.p(r = seq(0, 5, 0.01)),
    type = "l", xlab = "Radius (multiples of sigma)", ylab = "Probability")
lines(seq(0, 5, 0.01), circular.p(r = seq(0, 5, 0.01), detectfn = 2),
    type = "l", col = "red")
lines(seq(0, 5, 0.01), circular.p(r = seq(0, 5, 0.01), detectfn = 1,
    detectpar = list(sigma = 1,z = 4)), type = "l", col = "blue")
abline (h = 0.95, lty = 2)

legend (2.8, 0.3, legend = c("halfnormal","hazard-rate, z = 4", "exponential"),
    col = c("black","blue","red"), lty = rep(1,3))

## in this example, a more interesting comparison would use
## sigma = 0.58 for the exponential curve.

}

\keyword{ models }
