\name{rbind.popn}
\alias{rbind.popn}

\title{ Combine popn Objects }

\description{
Form a single \code{popn} object from two or more existing \code{popn} objects, or a list.
}

\usage{
rbind.popn(..., renumber = TRUE)
}
\arguments{
  \item{\dots}{ one or more \code{popn} objects, or a single list of \code{popn} objects }
  \item{renumber}{ logical for whether row names in the new object should be set to the row indices }
}

\details{
  An attempt to combine objects will fail if they conflict in their \code{covariates} attributes.  This is not an S3 method.
}

\value{
An object of class \code{popn} with number of rows equal to the sum of the rows in the input objects.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{popn}} }

\examples{

## generate and combine two subpopulations
trapobj <- make.grid()
p1 <- sim.popn(D = 3, core = trapobj)
p2 <- sim.popn(D = 2, core = trapobj)
covariates(p1) <- data.frame(size = rep('small', nrow(p1)))
covariates(p2) <- data.frame(size = rep('large', nrow(p2)))
pop <- rbind.popn(p1,p2)

}

\keyword{ manip }
