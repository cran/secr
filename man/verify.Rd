\name{verify}
\Rdversion{1.1}
\alias{verify}
\alias{verify.default}
\alias{verify.traps}
\alias{verify.capthist}
\alias{verify.mask}

\title{
Check SECR Data
}

\description{
Check that the data and attributes of an object are internally consistent to avoid crashing functions such as \code{secr.fit}
}

\usage{
\method{verify}{default}(object, report, ...)
\method{verify}{traps}(object, report = 2, ...)
\method{verify}{capthist}(object, report = 2, tol = 0.01, ...)
\method{verify}{mask}(object, report = 2, ...)
}

\arguments{
  \item{object}{
an object of class 'traps', 'capthist' or 'mask'
}
  \item{report}{
integer code for level of reporting to the console. 0 = no report, 1 = errors only, 2 = full.
}
  \item{tol}{
numeric tolerance for deviations from transect line (m)
}
  \item{...}{ 
other arguments (not used) 
}

}
\details{

Checks are performed specific to the class of 'object'. The default
method is called when no specific method is available (i.e. class not
'traps', capthist' or 'mask'), and does not perform any checks.

\code{verify.capthist}
\enumerate{ 
  \item No 'traps' component
  \item Invalid 'traps' component reported by verify.traps 
  \item No live detections
  \item Missing values not allowed in capthist
  \item Live detection(s) after reported dead
  \item More than one capture in single-catch trap(s)
  \item More than one detection per detector per occasion at proximity detector(s)
  \item Count detector counts less than zero
  \item Signal detector signal(s) less than threshold or invalid threshold
  \item Number of rows in 'traps' object not compatible with reported detections
  \item Number of rows in dataframe of individual covariates differs from capthist
  \item Number of occasions in usage matrix differs from capthist
  \item Detections at unused detectors
  \item Coordinates of detection(s) outside polygons ('polygon' detectors)
  \item Coordinates of detection(s) do not lie on any transect ('transect' detectors)
}

\code{verify.traps}
\enumerate{
  \item Missing detector coordinates not allowed
  \item Number of rows in dataframe of detector covariates differs from expected
  \item Number of detectors in usage matrix differs from expected
  \item Occasions with no used detectors
  \item Area detectors, but quadrats overlap or no area specified
}

\code{verify.mask}
\enumerate{
  \item Valid x and y coordinates
  \item Number of rows in covariates dataframe differs from expected
}

Earlier errors may mask later errors: fix & re-run.

}

\value{ 

A list with the component \code{errors}, a logical value indicating
whether any errors were found. If \code{object} contains multi-session
data then session-specific results are contained in a further list
component \code{bysession}.

Full reporting is the same as 'errors only' except that a message is
posted when no errors are found.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{
\code{\link{capthist}}, \code{\link{secr.fit}}
}

\examples{

data(secrdemo)
verify(captdata)

## create null (complete) usage matrix, and mess it up
temptraps <- make.grid()
usage(temptraps) <- matrix(1, nr = nrow(temptraps), nc = 5)
usage(temptraps)[,5] <- 0
verify (temptraps)

## create mask, and mess it up
tempmask <- make.mask(temptraps)
verify(tempmask)
tempmask[1,1] <- NA
verify(tempmask)

}

\keyword{ manip }
