\name{SPACECAP}
\alias{read.SPACECAP}
\alias{write.SPACECAP}

\title{
  Exchange data with SPACECAP package
}

\description{
  
  Data in a single-session \pkg{secr} \code{capthist} object may be written
  directly to the `csv' format used by \pkg{SPACECAP}, a package
  for Bayesian spatially explicit capture--recapture (Singh et
  al. 2010). Data in csv format may also be read to construct a
  \code{capthist} object for analysis in \pkg{secr}.
  
}

\usage{

write.SPACECAP(object, mask = NULL, buffer = 100, ndec = 2,
    filestem = "")
read.SPACECAP(AC, TD, detector = "proximity", session = "1")

}

\arguments{
  \item{object}{
    \code{capthist} object with the captures and trap locations to export
}
\item{mask}{
  \code{mask} object to use for state-space file
}
  \item{buffer}{ width of buffer in metres to use when creating a mask
  if none is specified
}
\item{ndec}{
  number of digits after decimal point for coordinates of mask on output
}
\item{filestem}{
  character value used to form names of output files
}
  \item{AC}{
  character value giving name of `animal capture' .csv file
}
  \item{TD}{
  character value giving name of `trap deployment' .csv file
}
  \item{detector}{
  detector type (`proximity' or `count')
}
  \item{session}{
  character value to use as session name
}
}
\details{
  If successful, \code{write.SPACECAP} creates three output files with names
  ending in `AC.csv',`TD.csv' and `SS.csv'. These are respectively
  the `Animal Capture', `Trap Deployment' and `State-Space' files
  required by \pkg{SPACECAP}.
}
\value{
 \code{write.SPACECAP} is used for its side effect of writing the
 required csv files.
 \code{read.SPACECAP} returns a \code{capthist} object.
}
\references{
  
Gopalaswamy, A.M., Royle, J.A., Hines, J.E., Singh, P., Jathanna, D.,
Kumar, N.S. and Karanth, K.U. (2012) Program SPACECAP: software for
estimating animal density using spatially explicit capture--recapture
models. \emph{Methods in Ecology and Evolution} \bold{3}, 1067--1072.

Singh, P., Gopalaswamy, A. M., Royle, A. J., Kumar, N. S. and Karanth,
K. U. (2010) SPACECAP: A program to estimate animal abundance and
density using Bayesian spatially explicit capture-recapture
models. Version 1.0. Wildlife Conservation Society - India Program, Centre for
Wildlife Studies, Bangalure, India.

}
\note{
  State-space csv files may be imported with \code{read.mask}.
}

\seealso{
  
  \code{\link{capthist}}, \code{\link{mask}}, \code{\link{read.mask}}
  
}

\examples{
\dontrun{

## coerce data to proximity detector type for export
demo <- reduce(captdata, output = "proximity") 
write.SPACECAP (demo, filestem = "demo")

## now read back the data just exported...
temp <- read.SPACECAP ("demoAC.csv", "demoTD.csv")
temp <- reduce(temp, output = "single")
summary (temp)
summary (captdata)

## should match exactly
identical(summary(temp), summary(captdata))

}
}

\keyword{ IO }

