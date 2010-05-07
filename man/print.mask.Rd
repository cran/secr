\name{print.mask}
\alias{print.mask}

\title{ Print Habitat Mask }

\description{
  Print the x-y coordinates of points in a \code{mask} object }

\usage{
\method{print}{mask}(x, ...)
}

\arguments{
\item{x}{ \code{mask} object }
\item{...}{ arguments passed to other functions }
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{mask}}, \code{\link{summary.mask}} }

\examples{
tempmask <- make.mask(make.grid())
print(tempmask)
}

\keyword{ print }

