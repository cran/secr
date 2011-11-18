\name{plot.capthist}
\alias{plot.capthist}

\title{ Plot Detection Histories }

\description{
  Display a plot of detection (capture) histories over a map of the detectors.
}

\usage{
\method{plot}{capthist}(x, rad = 5,
    hidetraps = FALSE, tracks = FALSE,
    title = TRUE, subtitle = TRUE, add = FALSE, varycol = TRUE, 
    icolours = NULL, randcol = FALSE, 
    lab1cap = FALSE, laboffset = 4, ncap = FALSE, 
    splitocc = NULL, col2 = "green",
    type = "petal",
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
  \item{icolours}{ vector of individual colours (when \code{varycol} =
  TRUE), or colour scale (non-petal plots) }
  \item{randcol}{ logical to use random colours (\code{varycol} = TRUE) }
  \item{lab1cap}{ logical for whether to label the first capture of each animal }
  \item{laboffset}{ distance by which to offset labels from points }
  \item{ncap}{ logical to display the number of detections per trap per occasion }
  \item{splitocc}{ optional occasion from which second colour is to be used  }
  \item{col2}{ second colour (used with \code{splitocc}) }
  \item{type}{ character string ("petal", "n.per.detector" or "n.per.cluster") }
  \item{cappar}{ list of named graphical parameters for detections (passed to \code{par}) }
  \item{trkpar}{ list of named graphical parameters for tracks (passed to \code{par}) }
  \item{labpar}{ list of named graphical parameters for labels (passed to \code{par}) }
  \item{\dots}{ arguments to be passed to \code{lines} if tracks are plotted }
}

\details{ 

By default, a `petal' plot is generated in the style of Density (Efford 2007)
using \code{eqscplot} from the MASS library. If \code{type =
"n.per.detector"} or \code{type = "n.per.cluster"} the result is a
colour-coded plot of the number of individuals at each unit, pooled over
occasions.

If \code{title} = FALSE no title is displayed; if \code{title} = TRUE,
the session identifer is used for the title.

If \code{subtitle} = FALSE no subtitle is displayed; if \code{subtitle}
= TRUE, the subtitle gives the numbers of occasions, detections and
individuals.

If \code{x} is a multi-session capthist object then a separate plot is
produced for each session. Use \code{par(mfrow = c(nr, nc))} to allow a
grid of plots to be displayed simultaneously (nr rows x nc columns).

These arguments are used only for petal plots: \code{rad},
\code{tracks}, \code{varycol}, \code{randcol}, \code{lab1cap},
\code{laboffset}, \code{ncap}, \code{splitocc}, \code{col2},
\code{trkpar}, and \code{labpar}.

If \code{icolours = NULL} and \code{varycol = TRUE} then a vector of
colours is generated automatically as terrain.colors((nrow(x)+1) * 1.5).
If there are too few values in \code{icolours} for the number of
individuals then colours will be re-used. 

}

\value{
  
  For \code{type = "petal"}, the number of detections in \code{x}.
  For \code{type = "n.per.detector"} or \code{type = "n.per.cluster"}, a
  dataframe with data for a legend (see Examples).

}

\references{ 

  Efford, M. G. (2007) \emph{Density 4.1: software for spatially
  explicit capture--recapture}.  Department of Zoology, University of
  Otago, Dunedin, New Zealand. \url{http://www.otago.ac.nz/density}

}

\seealso{ \code{\link{capthist}} }

\examples{
demotrap <- make.grid()
tempcapt <- sim.capthist(demotrap, 
    popn = list(D = 5, buffer = 50), 
    detectpar = list(g0 = 0.15, sigma = 30))
plot(tempcapt, border = 10, rad = 3, tracks = TRUE, 
    lab1cap = TRUE, laboffset = 2.5)


## type = n.per.cluster

## generate some captures
testregion <- data.frame(x = c(0,2000,2000,0),
    y = c(0,0,2000,2000))
popn <- sim.popn (D = 10, core = testregion, buffer = 0,
    model2D = "hills", details = list(hills = c(-2,3)))
t1 <- make.grid(nx = 1, ny = 1)
t1.100 <- make.systematic (cluster = t1, spacing = 100,
    region = testregion)
capt <- sim.capthist(t1.100, popn = popn, noccasions = 1)

## now plot captures ...
temp <- plot(capt, title = "Individuals per cluster",
    type = "n.per.cluster", hidetraps = FALSE,
    gridlines = FALSE, cappar = list(cex = 1.5))

\dontrun{
## add legend; click on map to place top left corner
legend (locator(1), pch = 21, pt.bg = temp$colour,
    pt.cex = 1.3, legend = temp$legend, cex = 0.8)
}

\dontrun{
## try varying individual colours - requires RColorBrewer
require(RColorBrewer)
plot(infraCH[[2]], icolours = brewer.pal(12, 'Set3'), tracks = T,
    bg = 'black', cappar = list(cex = 2), border = 10, rad = 2,
    gridlines=F)
}

}

\keyword{ hplot }
