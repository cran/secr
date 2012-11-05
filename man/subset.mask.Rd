\name{subset.mask}
\alias{subset.mask}
\alias{rbind.mask}

\title{ Subset Mask Object }

\description{
  Retain selected rows of a \code{mask} object.
}
\usage{

\method{subset}{mask}(x, subset, ...)

\method{rbind}{mask}(...)
}

\arguments{
  \item{x}{ \code{mask} object }
  \item{subset}{ numeric or logical vector to select rows of mask }
  \item{\dots}{ two or more \code{mask} objects (rbind only) }
}

\details{

The subscripts in \code{subset} may be of type integer, character or
logical as described in \code{\link{Extract}}.

Covariates are ignored by \code{rbind.mask}.

}

\value{

  For \code{subset}, an object of class `mask' with only the requested
  subset of rows and `type' attribute set to `subset'.

  For \code{rbind}, an object of class `mask' with all unique rows from
  the masks in \dots, and `type' attribute set to `rbind'.

}

\section{Warning}{

The spacing attribute is carried over from the input (it is not updated
automatically). In the case of very sparse masks (i.e. those with
isolated points) this may lead to an unexpected value for this
attribute. (Automatic updating requires excessive computation time
and/or memory for very large masks).
  
  }

\seealso{ \code{\link{mask}} }

\examples{
tempmask <- make.mask(make.grid())
OK <- (tempmask$x + tempmask$y) > 100
tempmask <- subset(tempmask, subset = OK)
plot(tempmask)

}

\keyword{ models }

