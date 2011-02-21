\name{subset.traps}
\alias{subset.traps}
\alias{split.traps}

\title{ Subset traps Object }

\description{
  Retain selected rows of a traps object.
}
\usage{
\method{subset}{traps}(x, subset, \dots)
\method{split}{traps}(x, f, drop = FALSE, prefix = 'S', \dots)
}

\arguments{
  \item{x}{ \code{traps} object }
  \item{subset}{ vector to subscript the rows of \code{x} }
  \item{\dots}{ arguments passed to other functions }
  \item{f}{ factor or object that may be coerced to a factor }
  \item{drop}{ logical indicating if levels that do not occur should be dropped (if f is a factor) }
  \item{prefix}{ a character prefix to be used for component names when values of f are numeric }
}
\details{

The subscripts in \code{subset} may be of type integer, character or
logical as described in \code{\link{Extract}}.

In the case of `polygon' and `transect' detectors, subsetting is done at
the level of whole polygons or transects. \code{subset} should therefore
have the same length as \code{levels(polyID(x))} or
\code{levels(transectID(x))}.

\code{split} generates a list in which each component is a
\code{traps} object. Each component corresponds to a level of
\code{f}.  The argument 'x' of \code{split} cannot be a list.

}
\value{

An object of class \code{traps} with only the requested subset of rows.
Subsetting is applied to \code{usage} and \code{covariates} attributes
if these are present.

}

\seealso{ \code{\link{traps}}, \code{\link{rbind.traps}} }

\examples{
## odd-numbered traps only, using modulo operator
temptrap <- make.grid(nx = 7, ny = 7)
t2 <- subset(temptrap, as.logical(1:nrow(temptrap) \%\% 2))
plot(t2)
}

\keyword{ manip }

