\name{read.mask}
\alias{read.mask}

\title{ Read Habitat Mask From File }
\description{
  Read coordinates of points on a habitat mask from a text file.
}
\usage{
read.mask(file, spacing = NULL, ...)
}

\arguments{
  \item{file}{ character string with name of text file }
  \item{spacing}{ spacing of grid points in metres }
  \item{...}{ other arguments to pass to \code{read.table} }
}
\details{
Usually, the x and y coordinates are the first two values on each line, separated by white space. If the file starts with a line of column headers and 'header = TRUE' is passed to \code{read.table} in the \dots argument then 'x' and 'y' need not be the first two fields.

If the grid cell size \code{spacing} is not provided then an attempt is
made to infer it from the minimum spacing of points.  This can be slow
and may demand more memory than is available. In rare cases (highly
fragmented masks) it may also yield the wrong answer.

}

\value{
object of class \code{mask} with type 'user'
}

\seealso{ \code{\link{mask}} }

\note{
The package \pkg{SPACECAP} uses a 'state-space' file in 'csv' text
format with columns 'X_COORD', 'Y_COORD' and 'HABITAT'. Such a file may
be input directly to \code{read.mask}; rows with HABITAT != 1 are dropped.
}

\examples{
## Replace file name with a valid local name and remove '#'
# read.mask (file = 'c:\\\\myfolder\\\\mask.txt', spacing = 3, header = TRUE)
## 'mask.txt' should have lines like this
# x   y
# 265 265
# 268 265
# ...
}

\keyword{ IO }

