\name{read.capthist}
\alias{read.capthist}
\alias{write.capthist}

\title{ Import or export data }

\description{

Data in the DENSITY formats for capture data and trap layouts may be
imported as a \code{capthist} object for analysis in \pkg{secr}. Data in
a \code{capthist} object may also be exported in these formats for use
in DENSITY (Efford 2009). \code{read.capthist} inputs data from text
files and constructs a \code{capthist} object in one step using the
functions \code{read.traps} and \code{make.capthist}.

}

\usage{

read.capthist(captfile, trapfile, detector = "multi", fmt = "trapID",
    noccasions = NULL, covnames = NULL, verify = TRUE, ...)

write.capthist(object, filestem = deparse(substitute(object)),
     sess = "1", ndec = 2, ...)

}

\arguments{

  \item{captfile}{ name of capture data file }
  \item{trapfile}{ name of trap layout file }
  \item{detector}{ character value for detector type (`single', `multi' etc.) }
  \item{fmt}{ character value for capture format (`XY' or `trapID') }
  \item{noccasions}{ number of occasions on which detectors were
    operated }
  \item{covnames}{ character vector of names for individual covariate fields }
  \item{verify}{ logical if TRUE then the resulting capthist object is
  checked with \code{\link{verify}}}

  \item{\dots}{ other arguments passed to \code{read.table}, \code{write.table} and
      \code{count.fields} }

  \item{object}{ \code{capthist} object with the captures and trap locations to export }
  \item{filestem}{ character value used to form names of output files }
  \item{sess}{ character session identifier }
  \item{ndec}{ number of digits after decimal point for x,y coordinates }

}

\details{

\code{read.capthist}

\code{captfile} should record one detection on each line. A detection
comprises a session identifier, animal identifier, occasion number (1,
2,...,\emph{S} where \emph{S} is the number of occasions), and a
detector identifier (\code{fmt = 'trapID'}) or X- and Y-coordinates
(\code{fmt = 'XY'}). Each line of \code{trapfile} has a detector
identifier and its X- and Y-coordinates. In either file type the
identifiers (labels) may be numeric or alphanumeric values. Values
should be separated by blanks or tabs unless (i) the file name ends in
`.csv' or (ii) \code{sep = ','} is passed in \dots, in which case commas
are assumed. Blank lines and any text after `#' are ignored. For further
details see \url{../doc/secr-datainput.pdf}, \code{\link{make.capthist}}
and `Data formats' in the help for DENSITY.

The \code{noccasions} argument is needed only if there were no
detections on the final occasion; it may be a positive integer (constant
across all sessions) or a vector of positive integers, one for each
session. \code{covnames} is needed only when \code{captfile} includes
individual covariates. Values of \code{noccasions} and \code{covnames}
are passed directly to \code{make.capthist}.

A session identifier is required even for single-session
capture data. In the case of data from multiple sessions,
\code{trapfile} may be a vector of file names, one for each session.

Additional data may be coded as for DENSITY. Specifically,
\code{captfile} may include extra columns of individual covariates, and
\code{trapfile} may code varying usage of each detector over occasions
and detector covariates. [These features have yet to be thoroughly
tested in 1.4.0].

\code{write.capthist}

  For a single-session analysis, DENSITY requires one text file of
  capture data and one text file with detector coordinates (the `trap
  layout' file).  \code{write.capthist} constructs names for these files
  by appending `capt.txt' and `trap.txt' to \code{filestem} which
  defaults to the name of the capthist object. If \code{filestem} is
  empty then output goes to the console.

  If \code{object} contains multiple sessions with differing
  \code{traps} then a separate trap layout file is exported for each
  session and each file name includes the session name. All capture data
  are exported to one file regardless of the number of sessions. The
  DENSITY format used is `TrapID'.

  Existing text files will be replaced without warning. In the case of a
  multi-session capthist file, session names are taken from
  \code{object} rather than \code{sess}. Session names are truncated to
  17 characters with blanks and commas removed.

  To export data in comma-delimited (`.csv') format, pass \code{sep =
  ','} in \dots. The resulting files have extension `.csv' rather than
  `.txt' and may be opened with spreadsheet software.

}

\references{

Efford, M. G. (2009) \emph{DENSITY 4.4: software for spatially explicit
capture--recapture}. Department of Zoology, University of Otago, Dunedin,
New Zealand \url{http://www.otago.ac.nz/density}.

}

\note{

  The DENSITY formats accommodate `single', `multi' and `proximity' data but
  not the newer detector types `count', `signal', `polygon', `transect'
  etc. See \code{\link{detector}} for more.

  The \dots argument is useful for some special cases. For example, if
  your input uses `;' instead of `#' for comments (`;' is also valid in
  DENSITY) then set \code{comment.char = ';'} in \code{read.capthist}.

  In a similar fashion, write comma- or tab-separated values by
  setting \code{sep = ','} or \code{sep = '\\t'} respectively.

  The arguments of \code{count.fields} are a subset of those of
  \code{read.table} so \dots is limited to any of \{sep, quote,
  skip, blank.lines.skip, comment.char\}.

  If you fail to set \code{fmt} correctly in \code{read.capthist} then the
  error message from \code{verify} may be uninformative.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{

  \code{\link{read.traps}}, \code{\link{make.capthist}},
  \code{\link{write.captures}}, \code{\link{write.traps}}, \code{\link{read.table}}

}

\examples{

data(ovenbird)
write.capthist(ovenCH)

}

\keyword{ IO }
