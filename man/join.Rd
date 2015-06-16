\name{join}
\alias{join}
\alias{unjoin}

\title{ Combine or Split Sessions of capthist Object }

\description{ Make a single-session capthist object from a list of
single-session objects, or a multi-session capthist object.  }

\usage{ 
join(object, remove.dupl.sites = TRUE, tol = 0.001) 
unjoin(object, interval, ...)
}

\arguments{
  
  \item{object}{list of single-session objects, or a multi-session
capthist object [\code{join}], or a single-session capthist object [\code{unjoin}]}

  \item{remove.dupl.sites}{ logical; if TRUE then a single record is
retained for each trap site used in multiple input sessions }

  \item{tol}{ absolute distance in metres within which sites are
considered identical }

  \item{interval}{ vector of times between occasions; zero indicates same session }

  \item{...}{ other arguments passed to \code{\link{subset.capthist}} }
}

\details{

  \describe{
\item{\code{join}}{

  The input sessions are assumed to be of the same detector type and to
  have the same attributes (e.g., covariates should be present for all
  or none).
  
  The number of occasions (columns) in the output is equal to the sum of
the number of occasions in each input.

A new dataframe of individual covariates is formed using the covariates
for the first occurrence of each animal.

Attributes xy and signal are handled appropriately, as is trap usage.
}
\item{\code{unjoin}}{

The input grouping of occasions (columns) into sessions is specified via
\code{interval}. This is a vector of length one less than the number of
occasions (columns) in \code{object}. Elements greater than zero
indicate a new session.

The \code{interval} argument may be omitted if \code{object} has a
valid `interval' attribute, as in the output from \code{join}.
}
}
}
\value{

For \code{join}, A single-session capthist object. The attribute `interval' records the
distinction between occasions that are adjacent in the input (interval =
0) and those that are in consecutive sessions (interval = 1); `interval'
has length one less than the number of occasions.

For \code{unjoin}, a multi-session capthist object. Sessions are named
with integers.

}

\note{

Do not confuse \code{unjoin} with \code{\link{split.capthist}} which
splits by row (animal) rather than by column (occasion).

Occasions survive intact; to pool occasions use
\code{\link{reduce.capthist}}.

\code{join} was modified in version 2.9.5 to check whether the
components of `object' all used the same detectors (`traps') (putting
aside differences in usage). If the traps are identical and
remove.dupl.sites = TRUE then the resulting `capthist' uses the common list of
detectors, with a usage attribute formed by concatenating the usage
columns of the input. This is faster than the previous filtering
algorithm using `tol'; the older algorithm is still used if the traps differ. 

}

\seealso{
  
  \code{\link{MS.capthist}}, \code{\link{rbind.capthist}}
  
}

\examples{

joined.ovenCH <- join (ovenCH)
summary(joined.ovenCH)
attr(joined.ovenCH, "interval")

summary(unjoin(joined.ovenCH))

}

\keyword{ manip }
