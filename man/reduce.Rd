\name{reduce}
\alias{reduce}

\title{ Combine Columns }

\description{
Combine columns in a matrix-like object to create a new data set using the first non-zero value.
}

\usage{
reduce (object, columns, ...)
}

\arguments{
  \item{object}{ object that may be coerced to a matrix }
  \item{columns}{ list in which each component is a vector of subscripts for columns to be pooled }
  \item{\dots}{ other arguments (not used currently) }
}

\details{
The first element of \code{columns} defines the columns of \code{object} for the first new column, the second for the second new column etc.
This is a generic method. A method exists for objects of class \code{capthist}.
}

\value{
A matrix with number of columns equal to \code{length(columns)}.
}

\seealso{ \code{\link{capthist}}, \code{\link{reduce.capthist}} }

\examples{

## matrix with random zeros
temp <- matrix(runif(20), nc = 4)
temp[sample(20,10)] <- 0
temp

reduce(temp, list(1:2, 3:4))

}

\keyword{ manip }

