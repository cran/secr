\name{reduce.capthist}
\alias{reduce.capthist}

\title{ Combine Occasions }

\description{

Combine columns (occasions) in a \code{capthist} object to create a new
data set, possibly converting between detector types

}

\usage{

\method{reduce}{capthist}(object, columns = NULL, outputdetector = 
    detector(traps(object)), select='last', dropunused = TRUE,
    verify = TRUE, sessions = NULL, \dots) 

}

\arguments{
  \item{object}{ \code{capthist} object }
  \item{columns}{ list in which each component is a vector of subscripts for occasions to be pooled }
  \item{outputdetector}{ character value giving detector type for output }
  \item{select}{ character value for method to resolve conflicts }
  \item{dropunused}{ logical, if TRUE any never-used detectors are dropped }
  \item{verify}{ logical, if TRUE the \code{verify} function is applied to the output }
  \item{sessions}{ vector of session indices or names (optional) }
  \item{\dots}{ other arguments (not used currently) }

}

\details{

The first component of \code{columns} defines the columns of
\code{object} for new occasion 1, the second for new occasion 2, etc. If
\code{columns} is NULL then all occasions are output. Subscripts in a
component of \code{columns} that do not match an occasion in the input
are ignored. When the output detector is one of the trap types
('single', 'multi'), reducing capture occasions can result in locational
ambiguity for individuals caught on more than one occasion, and for
single-catch traps there may also be conflicts between individuals at
the same trap. The method for resolving conflicts among 'multi'
detectors is determined by \code{select} which should be one of 'first',
'last' or 'random'. With 'single' detectors \code{select} is ignored and
the method is: first, randomly select* one trap per animal per day;
second, randomly select* one animal per trap per day; third, when
collapsing multiple days use the first capture, if any, in each
trap. With signal detectors, \code{select} determines how multiple
signal measurements are combined; options 'min', 'max' or 'mean' are
also allowed.

Usage data in the \code{traps} attribute are also
pooled if present; a trap is considered 'used' on a pooled occasion if
it was used on any contributing occasion.

* i.e., in the case of a single capture, use that capture; in the case of
multiple 'competing' captures draw one at random.

}

\value{

An object of class \code{capthist} with number of columns equal to
\code{length(occasions)}. The detector type is inherited from \code{object}
unless a new type is specified with the argument \code{outputdetector}.

}

\seealso{ \code{\link{capthist}}, \code{\link{subset.capthist}} }

\examples{
tempcapt <- sim.capthist (make.grid(nx = 6, ny = 6), nocc = 6)
class(tempcapt)

pooled.tempcapt <- reduce(tempcapt, col = list(1,2:3,4:6))
summary (pooled.tempcapt)
}

\keyword{ manip }

