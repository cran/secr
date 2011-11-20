\name{sort.capthist}
\alias{sort.capthist}

\title{
Sort Rows of capthist Object
}
\description{
Rows are sorted by fields in covariates or by a provided sort key of
length equal to the number of rows.
}
\usage{
\method{sort}{capthist}(x, decreasing = FALSE, by = "",
    byrowname = TRUE,...)
}

\arguments{
  \item{x}{\code{capthist} object
}
  \item{decreasing}{logical. Should the sort be increasing or decreasing? }
  \item{by}{character vector (names of covariates) or data frame whose
  columns will be used as sort keys}
  \item{byrowname}{logical. Should row name be used as a final sort key?}
  \item{\dots}{other arguments (not used)}
}
\details{

  For multi-session \code{capthist} objects only the named covariate form
is suitable as the number of rows varies between sessions.

If requested, rows are sorted by rowname within \code{by}. The effect of
the defaultsis to sort by rowname.

}
\value{
\code{capthist} object with sorted rows; any relevant attributes are
  also sorted (covariates, signal, xy)}

\seealso{
\code{\link{capthist}}
}
\examples{
sort(ovenCH, by = "Sex")
covariates(ovenCH)[["2005"]]
covariates(sort(ovenCH, by = "Sex"))[["2005"]]

}
\keyword{ manip }
