\name{detector}
\alias{detector}
\alias{detector<-}

\title{ Detector Type }

\description{
  Extract or replace the detector type. 
}

\usage{
detector(object, \dots)
detector(object) <- value
}

\arguments{
  \item{object}{ object with 'detector' attribute e.g. \code{traps} }
  \item{value}{ character string for detector type }
  \item{\dots}{ other arguments (not used) }
}

\details{

Valid detector types in version 1.3 are 'single', 'multi', 'proximity', 
'count', 'quadratbinary', 'quadratcount' and 'signal'. The detector type is 
stored as an attribute of a \code{traps} object.  Detector types are 
mostly described by Efford et al. (2009a,b). 

}

\value{
  character string for detector type
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ 
  \code{\link{traps}}
}

\references{

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009a) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009b) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\examples{
## Default detector type is 'multi'
temptrap <- make.grid(nx = 6, ny = 8)
detector(temptrap) <- 'proximity'
summary(temptrap)
}

\keyword{ models }

