\name{rbind.capthist}
\alias{rbind.capthist}
\alias{MS.capthist}

\title{ Combine capthist Objects }

\description{
Form a single \code{capthist} object from two or more compatible \code{capthist} objects.
}
\usage{
rbind.capthist(..., renumber = TRUE, pool = NULL)
MS.capthist(...)
}

\arguments{
\item{\dots}{ one or more simple \code{capthist} objects (i.e., single-session) }
\item{renumber}{ assign new composite individual ID: sourceobject.oldID }
\item{pool}{ list of vectors of indices }
}

\details{

In its simplest usage, the source objects in \dots each provide
detection histories from a single session, and the result is a
single-session object. For this to work the objects must be compatible.
\code{capthist} objects are compatible if they use the same detectors
(traps) and have consistent \code{covariates and other attributes. }

If the \dots argument includes at least one multi-session capthist
object then the elements will be formed into a single multi-session
capthist object. If \dots is a single multi-session object then the
components of \code{pool} are used to define combinations of old
sessions (e.g. pool = list(1:3, 4:5) forms two new sessions from 5 old ones).

Although \code{rbind.capthist} looks like an S3 method, it isn't. The
full function name must be used.

\code{MS.capthist} treats each source object as the data for a separate
session. Compatibility is not required. The \dots argument may include
lists of single-session capthist objects regardless of whether the list
has the class ('list','capthist').

}

\value{

For \code{rbind.capthist}, an object of class \code{capthist} with
number of rows equal to the sum of the rows in the input objects.

For \code{MS.capthist}, a multi-session object of class \code{capthist}
with number sessions equal to the number of objects in \dots.

}

\seealso{ \code{\link{capthist}}, \code{\link{subset.capthist}} }

\examples{

## simulate 2-part mixture
temptrap <- make.grid(nx = 8, ny = 8)
temp1 <- sim.capthist(temptrap,
    detectpar = list(g0 = 0.1, sigma = 40))
temp2 <- sim.capthist(temptrap,
    detectpar = list(g0 = 0.2, sigma = 20))
temp3 <- rbind.capthist(temp1, temp2)

## compare mixture to sum of components
## note 'traps visited' is not additive for 'multi' detector
## nor is 'traps set'
(summary(temp1)$counts +  summary(temp2)$counts) -
    summary(temp3)$counts

## assemble a multi-session object
## we fake the 2010 data by copying from 2005
## note must provide name each session
fakeCH <- ovenCH[['2005']]
MS.capthist(ovenCH, '2010' = fakeCH)

}

\keyword{ manip }

