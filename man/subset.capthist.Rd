\name{subset.capthist}
\alias{subset.capthist}
\alias{split.capthist}

\title{ Subset or Split capthist Object }

\description{
  Create a new \code{capthist} object or list of objects by selecting rows (individuals), columns (occasions) and traps from an existing \code{capthist} object.
}
\usage{
\method{subset}{capthist}(x, subset = NULL, occasions = NULL, traps = NULL,
    sessions = NULL, cutval = NULL, dropnullCH = TRUE, dropnullocc = FALSE,
    dropunused = TRUE, droplowsignals = TRUE, renumber = FALSE, ...)
\method{split}{capthist}(x, f, drop = FALSE, prefix = "S", ...)
}
\arguments{
  \item{x}{ object of class \code{capthist} }
  \item{subset}{ vector of subscripts to select rows (individuals) }
  \item{occasions}{ vector of subscripts to select columns (occasions) }
  \item{traps}{ vector of subscripts to select detectors (traps) }
  \item{sessions}{ vector of subscripts to select sessions }
  \item{cutval}{ new threshold for signal strength }
  \item{dropnullCH}{ logical for whether null (all-zero) capture histories
  should be dropped }
  \item{dropnullocc}{ logical for whether occasions with no detections should be dropped }
  \item{dropunused}{ logical for whether never-used detectors should be
    dropped }
  \item{droplowsignals}{ logical for whether cutval should be applied at each
  microphone rather than to sound as a whole } 
  \item{renumber}{ logical for whether row.names should be replaced with sequence number in new \code{capthist} }
  \item{f}{ factor or object that may be coerced to a factor }
  \item{drop}{ logical indicating if levels that do not occur should be dropped (if f is a factor) }
  \item{prefix}{ a character prefix to be used for component names when values of f are numeric }
  \item{...}{ other arguments (not used currently) }
}
\details{

Subscript vectors may be either logical- (length equal to the relevant
dimension of \code{x}), character- or integer-valued. Subsetting is
applied to attributes (e.g. \code{covariates}, \code{traps}) as
appropriate. The default action is to include all animals, occasions,
and detectors if the relevant argument is omitted.

When \code{traps} is provided, detections at other detectors are set to
zero, as if the detector had not been used, and the corresponding rows
are removed from \code{traps}. If the detector type is `proximity' then
selecting traps also reduces the third dimension of the capthist array. 

\code{split} generates a list in which each component is a
\code{capthist} object. Each component corresponds to a level of
\code{f}.

To combine (pool) occasions use \code{\link{reduce.capthist}}. There is
no equivalent of \code{unlist} for lists of \code{capthist} objects.

The effect of \code{droplowsignals = FALSE} is to retain below-threshold
measurements of signal strength on all channels (microphones) as long as
the signal is above \code{cutval} on at least one. In this case all
retained sounds are treated as detected on all microphones. This fails
when signals are already missing on some channels.

}

\value{

\code{capthist} object with the requested subset of observations, or a
list of such objects (i.e., a multi-session \code{capthist} object).
List input results in list output, except when a single session is
selected.

}

\seealso{ \code{\link{capthist}}, \code{\link{rbind.capthist}}, \code{\link{reduce.capthist}} }

\examples{
tempcapt <- sim.capthist (make.grid(nx=6, ny=6), nocc=6)
summary(subset(tempcapt, occ=c(1,3,5)))

## Consider `proximity' detections at a random subset of detectors
## This would not make sense for `multi' detectors, as the 
## excluded detectors influence detection probabilities in 
## sim.capthist.

tempcapt2 <- sim.capthist (make.grid(nx = 6, ny = 6, 
    detector = "proximity"), nocc = 6)
tempcapt3 <- subset(tempcapt2, traps = sample(1:36, 18, 
    replace=FALSE))
summary(tempcapt3)
plot(tempcapt3)

split (tempcapt2, f = sample (c("A","B"), nrow(tempcapt2), 
    replace = TRUE))

}

\keyword{ manip }

