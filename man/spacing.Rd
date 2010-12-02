\name{spacing}
\alias{spacing}
\alias{spacing<-}

\title{ Detector or Mask Spacing }

\description{
  Extract or replace the spacing attribute of a detector array or mask. 
}

\usage{
spacing(object, \dots)
spacing(object) <- value
}

\arguments{
  \item{object}{ object with `spacing' attribute e.g. \code{traps} }
  \item{value}{ numeric value for spacing }
  \item{\dots}{ other arguments (not used) }
}

\details{
  
The `spacing' attribute of a detector array is the average distance from
one detector to the nearest other detector.

The attribute was not always set by \code{make.grid()} and
\code{read.traps()} in versions of \pkg{secr} before 1.5.0. If the
attribute is found to be NULL then \code{spacing} will compute it on the
fly.
  
}

\value{
  scalar numeric value of mean spacing, or a vector if \code{object} has multiple sessions
}

\seealso{ 
  \code{\link{traps}}
}

\examples{
temptrap <- make.grid(nx = 6, ny = 8)
spacing(temptrap)
}

\keyword{ models }