\name{plot.traps}
\alias{plot.traps}

\title{ Plot traps Object }

\description{
  Map the locations of detectors (traps). 
}

\usage{
\method{plot}{traps}(x, border = 100, label = FALSE, offset = c(6,6), add = FALSE,
  hidetr = FALSE, detpar = list(), txtpar = list(), bg = "white",
  gridlines = TRUE, gridspace = 100, gridcol = "grey", 
  markvarying = FALSE, ...)
}

\arguments{
  \item{x}{ a \code{traps} object }
  \item{border}{ width of blank margin around the outermost detectors }
  \item{label}{ logical indicating whether a text label should appear by each detector }
  \item{offset}{ vector displacement of label from point on x and y axes }
  \item{add}{ logical to add detectors to an existing plot } 
  \item{hidetr}{ logical to suppress plotting of detectors }
  \item{detpar}{ list of named graphical parameters for detectors (passed to \code{par}) }
  \item{txtpar}{ list of named graphical parameters for labels (passed to \code{par}) }
  \item{bg}{ background colour }
  \item{gridlines}{ logical for plotting grid lines }
  \item{gridspace}{ spacing of gridlines }
  \item{gridcol}{ colour of gridlines }
  \item{markvarying}{ logical to distinguish detectors whose usage varies among ocasions }
  \item{\dots}{ arguments to pass to \code{eqscplot} }
}
\details{
\code{offset} may also be a scalar value for equal displacement on the x and y axes. The \code{hidetr} option is most likely to be used when plot.traps is called by plot.capthist. See \code{\link{par}} and \code{\link{colours}} for more information on setting graphical parameters. The initial values of graphical parameters are restored on exit.

Axes are not labeled. Use \code{\link{axis}} and \code{\link{mtext}} if necessary.
}

\value{
None
}

\seealso{ \code{\link{plot}}, \code{\link{traps}} }

\examples{
 temptrap <- make.grid()
 plot (temptrap, detpar = list(pch = 16, col = 'blue'), 
    label = TRUE, offset = 7)
}

\keyword{ hplot }

