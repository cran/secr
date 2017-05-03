\name{shareFactorLevels}
\alias{shareFactorLevels}

\title{Fix Inconsistent Covariates}

\description{

Factor covariates can give trouble in multi-session models if the levels differ among sessions. A warning is provided by \code{verify.capthist}. This function forces factor covariates to use the same levels. 

}

\usage{

shareFactorLevels(object, columns = NULL)

}

\arguments{
  \item{object}{multi-session capthist object or list of traps or mask objects}
  \item{columns}{indices of columns to fix (default all)}
}

\details{

Factor-valued covariates are coerced to use the same set of levels for each session of a multi-session capthist object or each component of a list of traps or masks. The combined level set is the union of all levels in separate sessions. The order of levels follows the default in \code{\link{factor}} (alphabetical).

}

\value{

An object of the same class as input.

}

\seealso{

\code{\link{verify.capthist}}

}

\keyword{ manip }
