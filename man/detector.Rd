\name{detector}
\alias{detector}
\alias{detector<-}
\alias{polygon}
\alias{polygonX}
\alias{transect}
\alias{transectX}
\alias{single}
\alias{multi}
%\alias{signal}
\alias{proximity}
\alias{count}

\title{ Detector Type }

\description{
  Extract or replace the detector type. 
}

\usage{
detector(object, \dots)
detector(object) <- value
}

\arguments{
  \item{object}{ object with `detector' attribute e.g. \code{traps} }
  \item{value}{ character string for detector type }
  \item{\dots}{ other arguments (not used) }
}

\details{

Valid detector types are `single', `multi', `proximity', `count',
`signal', `polygon', `transect', `polygonX', and `transectX'. The
detector type is stored as an attribute of a \code{traps} object.
Detector types are mostly described by Efford et al. (2009a,b; see also
\url{../doc/secr-overview.pdf}). Polygon and transect detector types are
for area and linear searches as described in
\url{../doc/secr-polygondetectors.pdf} and Efford (2011). The `signal'
detector type is used for acoustic data as described in
\url{../doc/secr-sound.pdf}.

}

\value{
  character string for detector type
}

\seealso{ 
  \code{\link{traps}}, \code{\link{RShowDoc}}
}

\references{

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area or transect searches. Unpublished
manuscript.

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

