\name{plot.mask}
\alias{plot.mask}

\title{ Plot Habitat Mask }

\description{
  Plot a habitat mask either as points or as an \code{image} plot. Colours may be used to show the value of one mask covariate.
}

\usage{
\method{plot}{mask}(x, border = 20, add = F, covariate = NULL, axes = F,
    dots = T, col = "grey", breaks = 12, ppoly = T, polycol = "red", \dots)
}

\arguments{
  \item{x}{ \code{mask} object }
  \item{border}{ width of blank display border (metres) }
  \item{add}{ logical for adding mask points to an existing plot }
  \item{covariate}{ name (as character string in quotes) or column number of a covariate to use for colouring }
  \item{axes}{ logical for plotting axes }
  \item{dots}{ logical for plotting mask points as dots, rather than as square pixels }
  \item{breaks}{ number of levels to use when cutting continuous covariate for plotting }
  \item{col}{ colour(s) to use for plotting }
  \item{ppoly}{ logical for whether the bounding polygon should be plotted (applies only if mask type = 'polygon') }
  \item{polycol}{ colour for outline of polygon (\code{ppoly} = TRUE) }
  \item{\dots}{ other arguments passed to \code{eqscplot} }

}

\details{

The argument \code{dots} selects between two distinct types of plot. If
using a covariate to colour points, the \code{col} argument should be a
colour vector of length equal to the number of levels.

Border lines around pixels are drawn in the current foreground colour
(par('fg')). Set this to NA with \code{par(fg=NA)} to eliminate borders,
but remember to reset it when you've finished.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{colours}}, \code{\link{mask}} }

\examples{

# simple

temptrap <- make.grid()
tempmask <- make.mask(temptrap)
plot (tempmask)

## restrict to points over an arbitrary detection threshold,
## add covariate, plot image and overlay traps

tempmask <- subset(tempmask, pdot(tempmask,temptrap)>0.001)

covariates (tempmask) <- data.frame(circle = 
    exp(-(tempmask$x^2 + tempmask$y^2)/10000) )

par(fg='white')
plot (tempmask, covariate = 'circle', dots = FALSE, axes = TRUE, 
    add = TRUE, breaks = 8, col = terrain.colors(8))
par(fg='black')

plot (temptrap, add = TRUE)

## add a legend

par(cex = 0.9)
covrange <- range(covariates(tempmask)$circle)
step <- diff(covrange)/8
colourlev <- terrain.colors(9)
zlev <- formatC(seq(covrange[1],covrange[2],step), format='f', 
    digits=2, width=4)
legend (x = 'topright', fill = colourlev, legend = zlev, 
    y.intersp = 0.8, title = 'Covariate')
title('Colour mask points with p.(X) > 0.001')
mtext(side=3,line=-1, 'g0 = 0.2, sigma = 20, nocc = 5')

}

\keyword{ hplot }

