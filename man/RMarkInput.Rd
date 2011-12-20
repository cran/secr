\name{RMarkInput}
\alias{RMarkInput}

\title{ Convert Data to RMark Input Format }

\description{
  
A single-session capthist object is formed into a dataframe that may be
passed directly to RMark.

}

\usage{

RMarkInput(object, grouped = TRUE)

}

\arguments{

  \item{object}{secr capthist object}

  \item{grouped}{logical for whether to replace each group of identical
    capture histories with a single line}
  
}
\details{

  To convert a multi-session object first collapse the sessions with
\code{\link{join}}.

}

\value{
  
  Dataframe with fields ch and freq. `ch' is a character string of 0's
  and 1's. If \code{grouped = FALSE} the rownames are retained and the
  value of `freq' is 1 or -1. Negative values of `freq' indicate
  removal.

  The attribute `interval' is copied from `object', if present;
  otherwise it is set to a vector of zeros (indicating a
  closed-population sample).

}

\references{

  Laake, J. and Rexstad E. (2008) Appendix C. RMark - an alternative
approach to building linear models in MARK. In: Cooch, E. and White,
G. (eds) Program MARK: A Gentle Introduction. 6th edition. Available at
\url{http://www.phidot.org/software/mark/docs/book/}.

}

\seealso{
\code{\link{join}}

}

\examples{

## ovenCH is a 5-year mist-netting dataset
ovenRD <- RMarkInput (join(ovenCH))
ovenRD

\dontrun{
## fit robust-design model in RMark
library(RMark)
ovenRD.data <- process.data(ovenRD, model = 'Robust',
    time.interval = attr(ovenRD, 'interval'))
ovenRD.model <- mark(data = ovenRD.data, model = 'Robust',
    model.parameters = list(p = list(formula = ~1, share = TRUE),
    GammaDoublePrime = list(formula = ~1),
    GammaPrime = list(formula = ~1),
    N = list(formula = ~1)))
cleanup(ask = FALSE)
}


}

\keyword{ manip }
