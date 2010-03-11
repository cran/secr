\name{make.mask}
\alias{make.mask}

\title{ Build Habitat Mask }

\description{
Construct a habitat mask object for spatially explicit capture-recapture. 
A mask object is a set of points with optional attributes. 
}

\usage{
make.mask(traps, buffer = 100, spacing = NULL, nx = 64, 
    type = "traprect", poly = NULL, pdotmin = 0.001, ...)
}

\arguments{
  \item{traps}{ object of class \code{traps} }
  \item{buffer}{ width of buffer in metres }
  \item{spacing}{ spacing between grid points (metres) }
  \item{nx}{ number of grid points in 'x' direction }
  \item{type}{ character string for method to use ('traprect', 'trapbuffer', 'pdot', 'polygon') }
  \item{poly}{ matrix of two columns interpreted as the x and y coordinates of a bounding polygon (optional) }
  \item{pdotmin}{ minimum detection probability for inclusion in mask when type = 'pdot' (optional) }
  \item{...}{ additional arguments passed to \code{pdot} when \code{type} = 'pdot' }
}

\details{

The 'traprect' method constructs a grid of points in the rectangle
formed by adding a buffer strip to the minimum and maximum x-y
coordinates of the detectors in \code{traps}. Both 'trapbuffer' and
'pdot' start with a 'traprect' mask and drop some points.

The 'trapbuffer' method restricts the grid to points within distance
\code{buffer} of any detector.

The 'pdot' method restricts the grid to points for which the net
detection probability \eqn{p.(\mathbf{X})}{p.(X)} (see
\code{\link{pdot}}) is at least \code{pdotmin}. Additional parameters
are used by \code{pdot} (detectpar, noccasions). Set these with the
\dots argument; otherwise \code{make.mask} will silently use the
arbitrary defaults.

The 'polygon' method places points on a rectangular grid clipped to the
polygon (\code{buffer} is not used).

If \code{spacing} is not specified then it is determined by dividing the
range of the x coordinates (including any buffer) by \code{nx}.

}

\value{
an object of class \code{mask}
}

\note{

A warning is displayed if type = 'pdot' and the buffer is too small to
include all points with \eqn{p.} > pdotmin.

A habitat mask is needed to fit an SECR model and for some related
computations. The default mask settings in \code{secr.fit} may be good
enough, but it is preferable to use \code{make.mask} to construct a mask
in advance and to pass that mask as an argument to \code{secr.fit}.


}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{mask}}, \code{\link{subset.mask}}, \code{\link{pdot}} }

\examples{

temptrap <- make.grid(nx = 10, ny = 10, spacing = 30)

## default method: traprect
tempmask <- make.mask(temptrap, spacing = 5)
plot(tempmask)
summary (tempmask)

## make irregular detector array by subsampling 
## form mask by 'trapbuffer' method
temptrap <- subset (temptrap, sample(nrow(temptrap), size = 30))
tempmask <- make.mask (temptrap, spacing = 5, type = 'trapbuffer')
plot (tempmask)
plot (temptrap, add = TRUE)

## form mask by 'pdot' method
temptrap <- make.grid(nx = 6, ny = 6)
tempmask <- make.mask (temptrap, buffer = 150, type = 'pdot', 
    pdotmin = 0.0001, detectpar = list(g0 = 0.1, sigma = 30),
    noccasions = 4)
plot (tempmask)
plot (temptrap, add = TRUE)

}

\keyword{ datagen }
