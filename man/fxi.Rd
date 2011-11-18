\name{fxi}
\alias{fxi.contour}
\alias{fxi.mode}
\alias{fxi.secr}
\alias{fxi}

\title{ Probability Density of Home Range Centre }

\description{
  Display contours of the probability density function for the estimated
  location of one or more range centres (f(Xi|wi)), compute values for
  particular points X, or compute mode of pdf.
}

\usage{
fxi.contour (object, i = 1, sessnum = 1, border = 100, nx = 64,
    levels = NULL, p = seq(0.1,0.9,0.1), plt = TRUE, add = FALSE,
    fitmode = FALSE, plotmode = FALSE, normal = TRUE, ...)
fxi.secr(object, i = 1, sessnum = 1, X, normal = TRUE)
fxi.mode(object, i = 1, sessnum = 1, start = NULL, ...)
}

\arguments{
  \item{object}{ a fitted secr model }
  \item{i}{ integer or character vector of individuals for which to plot contours,
  or a single individual as input to other functions }
  \item{sessnum}{ session number if \code{object$capthist} spans
  multiple sessions}
  \item{border}{ width of blank margin around the outermost detectors }
  \item{nx}{ dimension of interpolation grid in x-direction }
  \item{levels}{ numeric vector of confidence levels for Pr(X|wi)}
  \item{p}{ numeric vector of contour levels as probabilities }
  \item{plt}{logical to plot contours}
  \item{add}{logical to add contour(s) to an existing plot }
  \item{fitmode}{logical to refine estimate of mode of each pdf}
  \item{plotmode}{logical to plot mode of each pdf}
  \item{X}{2-column matrix of x- and y- coordinates}
  \item{normal}{ logical; should values of pdf be normalised?}
  \item{start}{ vector of x-y coordinates for maximization }
  \item{\dots}{additional arguments passed to \code{contour} or \code{nlm}}
}

\details{

  \code{fxi.contour} computes contours of probability density for one
  or more detection histories. Increase \code{nx} for smoother
  contours. If \code{levels} is not set, contour levels are set
  to approximate the confidence levels in n\code{p}.
  
  \code{fxi.secr} computes the probability density for a single
  detection history; \code{X} may contain coordinates for one or several
  points; a dataframe or vector (x then y) will be coerced to a matrix.

  \code{fxi.mode} finds the maximum of the pdf for a single detection history
  (i.e. \code{n} is of length 1). \code{fxi.mode} calls \code{\link[stats]{nlm}}.

  \code{fxi.contour} with \code{fitmode = TRUE} uses \code{fxi.mode} to
  find the maximum of each pdf. Otherwise the reported mode is an
  approximation (mean of coordinates of highest contour).

  If \code{i} is character it will be matched to row names of
  object$capthist (restricted to the relevant session in the case of a
  multi-session fit); otherwise it will be interpreted as a row number.
  
  Values of the pdf are optionally normalised by dividing by the
  integral of Pr(wi|X) over the habitat mask in \code{object}.
  
  If \code{start} is not provided then the first detector site is used,
  but this is not guaranteed to work.

  The \dots argument gives additional control over a contour plot; for
  example, set \code{drawlabels = FALSE} to suppress contour labels.
}
\value{
  
\code{fxi.contour} --
  
Coordinates of the plotted contours are returned as a list with one
component per polygon. The list is returned invisibly if plt = TRUE.

An additional component `mode' reports the x-y coordinates of the
highest point of each pdf (see Details).

\code{fxi.secr} --

Vector of probability densities

\code{fxi.mode} --

List with components `x' and `y'

}
\references{

  Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

}

\seealso{
\code{\link{pdot.contour}}, \code{\link[graphics]{contour}}
}

\note{
  These functions only work with homogeneous Poisson density models and
  do not allow incomplete usage (some detectors not used on some occasions).
}

\examples{

fxi.secr(secrdemo.0, i = 1, X = c(365,605))

plot(secrdemo.0$capthist)
## contour first 5 detection histories
fxi.contour (secrdemo.0, i = 1:5, add = TRUE,
    plotmode = TRUE, drawlabels = FALSE)

}

\keyword{ hplot }