\name{LLsurface.secr}
\alias{LLsurface.secr}
\title{
Plot likelihood surface
}
\description{
Calculate log likelihood over a grid of values of two beta parameters
from a fitted secr model and optionally make an approximate contour plot of the log
likelihood surface.
}
\usage{
LLsurface.secr(object, betapar = c("g0", "sigma"), xval = NULL, yval = NULL,
centre = NULL, realscale = TRUE, plot = TRUE, plotfitted = TRUE, ...)
}
\arguments{
  \item{object}{
\code{secr} object output from \code{secr.fit}
}
  \item{betapar}{
character vector giving the names of two beta parameters
}
  \item{xval}{
vector of numeric values for x-dimension of grid
}
  \item{yval}{
vector of numeric values for y-dimension of grid
}
\item{centre}{
  vector of central values for all beta parameters
}
\item{realscale}{
logical. If TRUE input and output of x and y is on the
untransformed (inverse-link) scale.
}
  \item{plot}{
logical. If TRUE a contour plot is produced
}
  \item{plotfitted}{
logical. If TRUE the MLE from \code{object} is shown on the plot (+)
}
  \item{\dots}{
other arguments passed to \code{\link{contour}}
}
}

\details{

\code{centre} is set by default to the fitted values of the beta
parameters in \code{object}. This has the effect of holding parameters
other than those in \code{betapar} at their fitted values.

If \code{xval} or \code{yval} is not provided then 11 values are set at
equal spacing between 0.8 and 1.2 times the values in \code{centre} (on
the `real' scale if \code{realscale} = TRUE and on the `beta' scale
otherwise).

Contour plots may be customized by passing graphical parameters through
the \dots argument.

}

\value{
Invisibly returns a matrix with the log likelihood evaluated at each
grid point
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\note{

  \code{LLsurface.secr} works for named `beta' parameters rather than
  `real' parameters. The default \code{realscale = TRUE} only works for
  beta parameters that share the name of the real parameter to which
  they relate i.e. the beta parameter for the base level of the real
  parameter. This is because link functions are defined for real
  parameters not beta parameters.

  The contours are approximate because they rely on
  interpolation. See Examples for a more reliable way to compare the
  likelihood at the MLE with nearby points on the surface.
}

\examples{
\dontrun{
data(secrdemo)
LLsurface.secr(secrdemo.CL, xval = seq(0.16,0.40,0.02),
    yval = 25:35, nlevels = 20)

## now verify MLE
## click on MLE and apparent 'peak'
xy <- locator(2)
temp <- LLsurface.secr(secrdemo.CL, xval = xy$x,
     yval = xy$y, plot = FALSE)
temp
}
}

\keyword{ models }
