\name{plot.capthist}
\alias{plot.capthist}

\title{ Plot Detection Histories }

\description{
  Display a plot of detection (capture) histories over a map of the detectors.
}

\usage{
\method{plot}{capthist}(x, rad = 5, hidetraps = FALSE, tracks = FALSE,
    title = TRUE, subtitle = TRUE, add = FALSE, varycol = TRUE, 
    icolours = NULL, randcol = FALSE, 
    lab1cap = FALSE, laboffset = 4, ncap = FALSE, 
    splitocc = NULL, col2 = "green", 
    cappar = list(cex = 1.3, pch = 16, col = "blue"), 
    trkpar = list(col = "blue", lwd = 1), 
    labpar = list(cex = 0.7, col = "black"), ...)
}

\arguments{
  \item{x}{ an object of class \code{capthist} }
  \item{rad}{ radial displacement of dot indicating each capture event from the detector location (used to separate overlapping points) }
  \item{hidetraps}{ logical indicating whether trap locations should be displayed }
  \item{tracks}{ logical indicating whether consecutive locations of individual animals should be joined by a line }
  \item{title}{ logical or character string for title }
  \item{subtitle}{ logical or character string for subtitle }
  \item{add}{ logical for whether to add to existing plot }
  \item{varycol}{ logical for whether to distinguish individuals by colour }
  \item{icolours}{ vector of individual colours (when \code{varycol} = TRUE) }
  \item{randcol}{ logical to use random colours (\code{varycol} = TRUE) }
  \item{lab1cap}{ logical for whether to label the first capture of each animal }
  \item{laboffset}{ distance by which to offset labels from points }
  \item{ncap}{ logical to display the number of detections per trap per occasion }
  \item{splitocc}{ optional occasion from which second colour is to be used  }
  \item{col2}{ second colour (used with \code{splitocc}) }
  \item{cappar}{ list of named graphical parameters for detections (passed to \code{par}) }
  \item{trkpar}{ list of named graphical parameters for tracks (passed to \code{par}) }
  \item{labpar}{ list of named graphical parameters for labels (passed to \code{par}) }
  \item{\dots}{ arguments to be passed to \code{lines} if tracks are plotted }
}

\details{ 

A plot is generated in the style of Density (Efford 2007)
using \code{eqscplot} from the MASS library. 

If \code{title} = FALSE no title is displayed; if \code{title} = TRUE,
the session identifer is used for the title.

If \code{subtitle} = FALSE no subtitle is displayed; if \code{subtitle}
= TRUE, the subtitle gives the numbers of occasions, detections and
individuals.

If \code{icolours = NULL} and \code{varycol = TRUE} then a vector of
colours is generated automatically as terrain.colors((nrow(x)+1) * 1.5).
If there are too few values in \code{icolours} for the number of
individuals then colours will be re-used. 

If \code{x} is a multi-session capthist object then a separate plot is
produced for each session. Use \code{par(mfrow = c(nr, nc))} to allow a
grid of plots to be displayed simultaneously (nr rows x nc columns).

}

\value{
  The number of detections in \code{x}.
}

\references{ 
Efford, M. G. (2007) \emph{Density 4.1: software for spatially explicit capture--recapture}.  Department of Zoology, University of Otago, Dunedin, New Zealand. \url{http://www.otago.ac.nz/density}
}

\seealso{ \code{\link{capthist}} }

\examples{
demotrap <- make.grid()
tempcapt <- sim.capthist(demotrap, 
    popn = list(D = 5, buffer = 50), 
    detectpar = list(g0 = 0.15, sigma = 30))
plot(tempcapt, border = 10, rad = 3, tracks = TRUE, 
    lab1cap = TRUE, laboffset = 2.5)
}

\keyword{ hplot }
