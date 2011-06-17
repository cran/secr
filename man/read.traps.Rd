\name{read.traps}
\alias{read.traps}
\title{ Read Detector Data From File }

\description{
  Construct an object of class \code{traps} with detector locations from a text file or data frame. Usage per occasion and covariates may be included.
}

\usage{
read.traps(file = NULL, data = NULL, detector = "multi", covnames =
NULL, ...)
}

\arguments{
  \item{file}{ character string with name of text file }
  \item{data}{ data frame of detector coordinates }
  \item{detector}{ character string for detector type }
  \item{covnames}{ character vector of names for detector covariate
  fields }
  \item{\dots}{ other arguments to pass to \code{read.table} }
}

\details{

Reads a text file in which the first column is a character string
identifying a detector and the next two columns are its x- and
y-coordinates, separated by white space. The coordinates optionally may
be followed by a string of codes `0' or `1' indicating whether the
detector was operated on each occasion. Trap-specific covariates may be
added at the end of the line preceded by `/'. This format is compatible
with the Density software (Efford 2007), except that all detectors are
assumed to be of the same type (usage codes greater than 1 are treated
as 1), and more than one covariate may be specified.

If \code{file} is missing then x-y coordinates will be taken instead
from \code{data}. This option does not allow for \code{covariates} or
\code{usage}, but they maybe added later.

\code{detector} specifies the behaviour of the detector following Efford
et al. (2009). `single' refers to a trap that is able to catch at most
one animal at a time; `multi' refers to a trap that may catch more than
one animal at a time. For both `single' and `multi' detectors a trapped
animals can appear at only one detector per occasion. Detectors of type
`proximity', such as camera traps and hair snags for DNA sampling, allow
animals to be recorded at several detectors on one occasion. See
\code{\link{detector}} for further detector types.

For polygon and transect detector types, each line corresponds to a
vertex and starts with a code to identify the polygon or transect (hence
the same code appears on 2 or more lines). For input from a dataframe
the code column should be named `polyID'. Also, usage and covariates
are for the polygon or transect as a whole and not for each
vertex. Usage and covariates are appended to the end of the line, just
as for point detectors (traps etc.). The usage and covariates for each
polygon or transect are taken from its first vertex. Although the
end-of-line strings of other vertices are not used, they cannot be blank
and should use the same spacing as the first vertex.

}

\value{

An object of class \code{traps} comprising a data frame of x- and
y-coordinates, the detector type (`single', `multi', `proximity',
`count', `polygon' etc.), and possibly other attributes.

}

\references{ 

Efford, M. G. (2007) \emph{Density 4.1: software for spatially explicit
capture--recapture}. Department of Zoology, University of Otago,
Dunedin, New Zealand. \url{http://www.otago.ac.nz/density}

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\seealso{  \code{\link{traps}}, \code{\link{make.grid}}, \code{\link{detector}} }

\examples{
## Replace file name with a valid local name and remove '#' 
# read.traps ('c:\\\\myfolder\\\\mytraps.txt', detector='proximity')
## 'mytraps.txt' should have lines like this 
# 1     365     365
# 2     365     395
# 3     365     425
# etc.

}
\keyword{ manip }

