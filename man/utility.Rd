\name{getMeanSD}
\alias{getMeanSD}
\alias{maskarea}
\alias{masklength}

\title{Utility Functions}

\description{
  
Minor function(s) exported for use in related packages.

}

\usage{

getMeanSD(xy)
maskarea(mask, sess = 1)
masklength(mask, sess = 1)

}
\arguments{
  
  \item{xy}{2-column matrix or dataframe}
  \item{mask}{mask or linearmask object}
  
}

\details{
  
  For \code{masklength} the input should be a linear mask from \pkg{secrlinear}.

}

\value{

  For \code{getMeanSD}, a dataframe with columns `x' and `y' and two
  rows, mean and SD.

  For \code{maskarea}, the summed area of mask cells in hectares (ha).

  For \code{masklength}, the summed length of mask cells in kilometers (km).

}

\examples{

getMeanSD(possummask)

}

\keyword{ internal }
