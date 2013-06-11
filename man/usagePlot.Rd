\name{usagePlot}
\alias{usagePlot}

\title{
Plot Usage
}

\description{

This function displays variation in effort (usage) over detectors as a
bubble plot (circles with radius scaled so that area is proportional to
effort).

}

\usage{

usagePlot(object, add = FALSE, occasion = NULL, col = "black", fill =
FALSE, scale = 2, metres = TRUE, rad = 5, ...)

}

\arguments{
  \item{object}{traps object with usage attribute}
  \item{add}{logical; if FALSE plot.traps is called to create a base plot}
  \item{occasion}{integer number of the occasion for which effort is
  plotted, or NULL}
  \item{col}{character or integer colour value}
  \item{fill}{logical; if TRUE the circle is filled with the line colour}
  \item{scale}{numeric value used to scale radius}
  \item{metres}{logical; if TRUE scale is a value in metres (see Details)}
  \item{rad}{numeric; radial displacement of symbol centre for each occasion from
  true detector location (metres) }
  \item{\dots}{other arguments passed to plot.traps}
}

\details{

By default (\code{occasion = NULL}) circles representing usage on each
  occasion are plotted around the detector location at distance
  \code{rad}, as in the petal plot of
  \code{\link{plot.capthist}}. Otherwise, the usage on a single
  specified occasion is plotted as a circle centred at the detector
  location.
  
The \code{metres} argument switches between two methods. If \code{metres
= TRUE}, the \code{symbols} function is used with \code{inches = FALSE} to plot
circles with radius scaled in the units of \code{object} (i.e. metres;
\code{scale} is then the radius in metres of the symbol for a detector
with usage = 1.0).  Otherwise, plotting uses \code{points}; this has the
advantage of producing better filled circles, but a suitable value of
scale must be found by trial and error.
  
Package \pkg{sp} provides an alternative (see Examples).

}

\value{

  No value is returned.
  
}

\seealso{
  
  \code{\link{usage}}, \code{\link{symbols}}, code{\link[sp]{bubble}}
  
}

\examples{

simgrid <- make.grid(nx = 10, ny = 10, detector = 'proximity')
usage(simgrid) <- matrix(rep(1:10, 50), nrow = 100, ncol = 5)
usagePlot(simgrid, border = 20, scale = 1.5, fill = FALSE,
  metres = FALSE)

# It is hard to get the legend just right
# here is one attempt
legend (x = -50, y = 185, legend = c(1,2,5,10), pch = 1, pt.cex =
  c(1,2,5,10)^0.5 * 1.5, x.intersp = 3, y.intersp = 1.8, adj = 1,
  bty = 'n', title = 'Usage')

usagePlot(simgrid, occasion = NULL, border = 20, scale = 1.5, fill = FALSE,
  metres = FALSE)

\dontrun{
# bubble plot in package 'sp'
library(sp)
simgrid$usage <- usage(simgrid)[,1]  ## occasion 1
class(simgrid) <- 'data.frame'
coordinates(simgrid) <- c('x','y')
bubble(simgrid)   
}

}


\keyword{ hplot }
