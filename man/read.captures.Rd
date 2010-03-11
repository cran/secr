\name{read.captures}
\alias{read.captures}

\title{ Read Capture Data From File }
\description{
  Input raw detection (capture) data for spatially explicit capture--recapture.
}
\usage{
read.captures(file, ...)
}

\arguments{
  \item{file}{ character string with name of text file }
  \item{\dots}{ other arguments to be passed to \code{read.table} }
}
\details{
  This brief function merely calls \code{read.table} and is therefore technically redundant.  It accepts whatever data are in \code{file}.  The result is used as input to \code{make.capthist}.
}
\value{
  A dataframe.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{make.capthist}} }
\examples{
\dontrun{
   ## assumes 'capt.txt' in working folder
   temp <- read.captures(file = 'capt.txt', col.names = 
               c('Session','ID','Occasion', 'X', 'Y'))
}
}

\keyword{ IO }

