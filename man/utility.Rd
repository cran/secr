\name{utility}
\alias{getMeanSD}
\alias{maskarea}
\alias{masklength}
\alias{edist}

\title{Utility Functions}

\description{
  
Minor function(s) exported for use in related packages.

}

\usage{

getMeanSD(xy)
maskarea(mask, sess = 1)
masklength(mask, sess = 1)
edist(xy1, xy2) 

}
\arguments{
  
  \item{xy}{2-column matrix or dataframe}
  \item{xy1}{2-column matrix or dataframe}
  \item{xy2}{2-column matrix or dataframe}
  \item{mask}{mask or linearmask object}
  
}

\details{
  
  For \code{masklength} the input should be a linear mask from \pkg{secrlinear}.
  
  \code{edist} computes the Euclidean distance between each point in xy1 and each point in xy2. (This duplicates the functionality of `rdist' in package \pkg{fields}).

}

\value{

  For \code{getMeanSD}, a dataframe with columns `x' and `y' and two
  rows, mean and SD.

  For \code{maskarea}, the summed area of mask cells in hectares (ha).

  For \code{masklength}, the summed length of mask cells in kilometers (km).
  
  For \code{edist}, a matrix with dim = c(nrow(xy1), nrow(xy2)).

}

\examples{

getMeanSD(possummask)

}

\keyword{ internal }
