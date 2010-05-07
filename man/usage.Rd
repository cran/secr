\name{usage}
\alias{usage}
\alias{usage<-}

\title{ Detector Usage }

\description{
  Extract or replace usage information of a \code{traps} object.
}

\usage{
usage(object, \dots)
usage(object) <- value
}
\arguments{
  \item{object}{ a \code{traps} object }
  \item{value}{ a matrix of traps x occasions 1 if trap[i] used on occasion[j], zero otherwise. }
  \item{\dots}{ other arguments (not used) }
}
\details{
  For replacement, the number of rows of \code{value} must match exactly the number of traps in \code{object}.
}
\value{
  usage(object) returns the usage matrix of the \code{traps} object. \code{usage(object)} may be NULL.
}
\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{traps}} }

\examples{
demo.traps <- make.grid(nx = 6, ny = 8)
## random usage over 5 occasions
usage(demo.traps) <- matrix (sample(0:1, 48*5, replace = TRUE, 
    p = c(0.5,0.5)), nc = 5)
usage(demo.traps)
summary(demo.traps)
}
\keyword{ manip }

