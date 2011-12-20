\name{reduce.capthist}
\alias{reduce.capthist}

\title{ Combine Occasions Or Convert Detector Types }

\description{

Combine columns (occasions) in a \code{capthist} object to create a new
data set, possibly converting between detector types

}

\usage{

\method{reduce}{capthist}(object, columns = NULL, outputdetector = 
    detector(traps(object)), select='last', dropunused = TRUE,
    verify = TRUE, sessions = NULL, by = 1, \dots) 

}

\arguments{
  \item{object}{ \code{capthist} object }
  \item{columns}{ list in which each component is a vector of subscripts for occasions to be pooled }
  \item{outputdetector}{ character value giving detector type for output }
  \item{select}{ character value for method to resolve conflicts }
  \item{dropunused}{ logical, if TRUE any never-used detectors are dropped }
  \item{verify}{ logical, if TRUE the \code{verify} function is applied to the output }
  \item{sessions}{ vector of session indices or names (optional) }
  \item{by}{ number of old occasions in each new occasion }
  \item{\dots}{ other arguments (not used currently) }

}

\details{

The first component of \code{columns} defines the columns of
\code{object} for new occasion 1, the second for new occasion 2, etc. If
\code{columns} is NULL then all occasions are output. Subscripts in a
component of \code{columns} that do not match an occasion in the input
are ignored. When the output detector is one of the trap types
(`single', `multi'), reducing capture occasions can result in locational
ambiguity for individuals caught on more than one occasion, and for
single-catch traps there may also be conflicts between individuals at
the same trap. The method for resolving conflicts among `multi'
detectors is determined by \code{select} which should be one of `first',
`last' or `random'. With `single' detectors \code{select} is ignored and
the method is: first, randomly select* one trap per animal per day;
second, randomly select* one animal per trap per day; third, when
collapsing multiple days use the first capture, if any, in each
trap.

Usage data in the \code{traps} attribute are also
pooled if present; a trap is considered `used' on a pooled occasion if
it was used on any contributing occasion.

* i.e., in the case of a single capture, use that capture; in the case of
multiple `competing' captures draw one at random.

If \code{columns} is not provided then old occasions are grouped into
new occasions as indicated by the \code{by} argument. For example, if
there are 15 old occasions and \code{by = 5} then new occasions will be
formed from occasions 1:5, 6:10, and 11:15. A warning is given when the
number of old occasions is not a multiple of \code{by} as then the final
new occasion will comprise fewer old occasions.

A special use of the \code{by} argument is to combine all occasions into
one for each session in a multi-session dataset. This is done by setting
\code{by = "all"}.

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

pooled.tempcapt2 <- reduce(tempcapt, by = 2)
summary (pooled.tempcapt2)

## collapse multi-session dataset to single-session 'open population'
onesess <- join(reduce(ovenCH, by = "all"))
summary(onesess)
}

\keyword{ manip }

