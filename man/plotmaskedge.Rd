\name{plotMaskEdge}
\alias{plotMaskEdge}

\title{Outline Around Mask Cells}

\description{
  
Plots the outer edge of a  mask.

}
\usage{

plotMaskEdge(mask, add = FALSE, ...)

}

\arguments{
  \item{mask}{secr habitat mask object}
  \item{add}{logical; if TRUE the line is added to an existing plot}
  \item{\dots}{other line plotting arguments passed to \code{\link{segments}}}
}

\details{
  
  May be slow. 
  
}

\value{
  None
}

\examples{

\dontrun{
plot(possummask)
plotMaskEdge (possummask, add = TRUE)
}

}

\keyword{ hplot }
