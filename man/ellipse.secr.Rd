\name{ellipse.secr}
\alias{ellipse.secr}
\title{
Confidence ellipse
}
\description{
Plot joint confidence ellipse for two parameters of secr model
}
\usage{
ellipse.secr(object, par = c("g0", "sigma"), alpha = 0.05,
    npts = 100, plot = TRUE, linkscale = TRUE, add = FALSE,
    col = palette(), ...)
}

\arguments{
  \item{object}{
\code{secr} object output from \code{secr.fit}
}
  \item{par}{
character vector of length two, the names of two 'beta' parameters
}
  \item{alpha}{
alpha level for confidence intervals 
}
  \item{npts}{
number of points on perimeter of ellipse
}
  \item{plot}{
logical for whether ellipse should be plotted
}
  \item{linkscale}{
logical; if FALSE then coordinates will be backtransformed from the link scale
}
  \item{add}{
logical to add ellipse to an existing plot 
}
  \item{col}{
vector of one or more plotting colours
}
  \item{\dots}{
arguments to pass to plot functions
}
}
\details{
A confidence ellipse is calculated from the asymptotic variance-covariance matrix of the beta parameters (coefficients), and optionally plotted.

If \code{linkscale} == FALSE, the inverse of the appropriate link transformation is applied to the coordinates of the ellipse, causing it to deform.

If \code{object} is a list of secr models then one ellipse is constructed for each model. Colours are recycled as needed.
}
\value{
A list containing the x and y coordinates is returned invisibly
}
\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\examples{
data(secrdemo)
ellipse.secr(secrdemo.0)
}
\keyword{ hplot }
