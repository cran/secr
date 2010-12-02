\name{contour}
\alias{pdot.contour}
\alias{buffer.contour}

\title{
Contour Detection Probability
}
\description{

  Display contours of the net probability of detection p.(X), or the
area within a specified distance of detectors.  \code{buffer.contour}
adds a conventional `boundary strip' to a detector (trap) array, where
\code{buffer} equals the strip width.

}
\usage{

pdot.contour(traps, border = NULL, nx = 64, detectfn = 0,
    detectpar = list(g0 = 0.2, sigma = 25, z = 1), noccasions = 5,
    levels = seq(0.1, 0.9, 0.1), poly = NULL, plt = TRUE,
    add = FALSE, ...)

buffer.contour(traps, buffer, nx = 64, convex = FALSE, ntheta = 100,
     plt = TRUE, add = FALSE, poly = NULL, ...)

}

\arguments{
  \item{traps}{ \code{traps} object }
  \item{border}{ width of blank margin around the outermost detectors }
  \item{nx}{ dimension of interpolation grid in x-direction }
  \item{detectfn}{ integer code or character string for shape of
    detection function 0 = halfnormal etc. -- see \link{detectfn}} 
  \item{detectpar}{ list of values for named parameters of detection
  function }
  \item{noccasions}{ number of sampling occasions}
  \item{levels}{ vector of levels for p.(X)}
  \item{poly}{ matrix of two columns, the x and y
    coordinates of a bounding polygon (optional) }
  \item{plt}{ logical to plot contours }
  \item{add}{ logical to add contour(s) to an existing plot }
  \item{\dots}{ other arguments to pass to \code{contour} }
  \item{buffer}{ vector of buffer widths }
  \item{convex}{ logical, if TRUE the plotted contour(s) will be convex }
  \item{ntheta}{ integer value for smoothness of convex contours }
}
\details{
\code{pdot.contour} constructs a rectangular mask and applies\code{\link{pdot}} to
compute the p.(X) at each mask point.

if \code{convex = FALSE}, \code{buffer.contour} constructs a mask and
contours the points on the basis of distance to the nearest trap at the
levels given in \code{buffer}.

if \code{convex = TRUE}, \code{buffer.contour} constructs a set of
potential vertices by adding points on a circle of radius =
\code{buffer} to each detector location; the desired contour is the
convex hull of these points (this algorithm derives from Efford, 2009).

}
\value{
  
  Coordinates of the plotted contours are returned as a list with one
  component per polygon. The list is returned invisibly if \code{plt =
  TRUE}.  
  
}

\references{

  Efford, M. G. (2009) \emph{DENSITY 4.4: software for spatially
explicit capture--recapture}. Department of Zoology, University of
Otago, Dunedin, New Zealand \url{http://www.otago.ac.nz/density}.

}

\note{
 
  The precision (smoothness) of the fitted line in \code{buffer.contour}
  is controlled by \code{ntheta} rather than \code{nx} when \code{convex
    = TRUE}.

  To suppress contour labels, include the argument \code{drawlabels =
  FALSE} (this will be passed via \dots to \code{contour}). Other useful
  arguments of \code{contour} are \code{col} (colour of contour lines)
  and \code{lwd} (line width).
  
}

\seealso{

  \code{\link{pdot}}, \code{\link{make.mask}}
  
}
\examples{

possumtraps <- traps(possumCH)

plot(possumtraps, border = 270)
pdot.contour(possumtraps, detectfn = 0, nx = 128, detectpar =
    detectpar(possum.model.1), levels = c(0.1, 0.01, 0.001),
    add = TRUE)

## clipping to polygon
olddir <- setwd(system.file("extdata", package = "secr"))
possumarea <- read.table('possumarea.txt', header = TRUE)
oldpar <- par(xpd = TRUE, mar = c(1,6,6,6))
plot(possumtraps, border = 400, gridlines = FALSE)
pdot.contour(possumtraps, detectfn = 0, nx = 256, detectpar =
    detectpar(possum.model.1), levels = c(0.1, 0.01, 0.001),
    add = TRUE, poly = possumarea, col = 'blue')
lines(possumarea)
setwd(olddir)
par(oldpar)

## convex and concave buffers
plot(possumtraps, border = 270)
buffer.contour(possumtraps, buffer = 100, add = TRUE, col = 'blue')
buffer.contour(possumtraps, buffer = 100, convex = TRUE, add = TRUE)

## areas
buff.concave <- buffer.contour(possumtraps, buffer = 100,
    plt = FALSE)
buff.convex  <- buffer.contour(possumtraps, buffer = 100,
    plt = FALSE, convex = TRUE)
sum (sapply(buff.concave, polyarea)) ## sum over parts
sapply(buff.convex, polyarea)

## effect of nx on area
buff.concave2 <- buffer.contour(possumtraps, buffer = 100,
    nx = 128, plt = FALSE)
sum (sapply(buff.concave2, polyarea))

}

\keyword{ hplot }
