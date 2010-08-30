\name{homerange}
\alias{dbar}
\alias{RPSV}
\alias{MMDM}
\alias{ARL}
\alias{moves}

\title{ Home Range Statistics }

\description{

Some ad hoc measures of home range size may be calculated in \pkg{secr}
from capture--recapture data:

\code{dbar} is the mean distance between consecutive capture locations,
pooled over individuals (e.g. Efford 2004). \code{moves} returns the
raw distances.

\code{RPSV} (for 'Root Pooled Spatial Variance') is a measure of the 2-D
dispersion of the locations at which individual animals are detected,
pooled over individuals.

\code{MMDM} (for 'Mean Maximum Distance Moved') is the average maximum
distance between detections of each individual i.e. the observed range
length averaged over individuals (Otis et al. 1978).

\code{ARL} or 'Asymptotic Range Length') is obtained by fitting an
exponential curve to the scatter of observed individual range length vs
the number of detections of each individual (Jett and Nichols 1987: 889).

}

\usage{
dbar(capthist)
RPSV(capthist)
MMDM(capthist, min.recapt = 1, full = FALSE)
ARL(capthist, min.recapt = 1, plt = FALSE, full = FALSE)
moves(capthist)
}

\arguments{
  \item{capthist}{ object of class \code{\link{capthist}} }
  \item{min.recapt}{ integer minimum number of recaptures for a
  detection history to be used }
  \item{plt}{ logical; if TRUE observed range length is plotted against number
  of recaptures }
  \item{full}{ logical; set to TRUE for detailed output }
}

\details{
 \code{dbar} is defined as
 \deqn{
\overline{d}=\frac{\sum\limits _{i=1}^{n}
                          \sum\limits _{j=1}^{n_i - 1}
                              \sqrt{(x_{i,j}-x_{i,j+1})^2 + (y_{i,j}-y_{i,j+1})^2}}
                         {\sum\limits _{i=1}^{n} (n_i-1)}}{see pdf manual}

 \code{RPSV} is defined as
 \deqn{
 RPSV = \sqrt{
    \frac {\sum\limits _{i=1}^{n} \sum\limits _{j=1}^{n_i} [
	    (x_{i,j} - \overline x_i)^2 + (y_{i,j} - \overline y_i)^2
		]}{\sum\limits _{i=1}^{n} (n_i-1) - 1}}
}{see pdf manual}

\code{dbar} and \code{RPSV} have a specific role as proxies for
detection scale in inverse-prediction estimation of density (Efford
2004; see \code{\link{ip.secr}}).

\code{RPSV} is used in \code{autoini} to obtain plausible starting
values for maximum likelihood estimation.

\code{MMDM} and \code{ARL} discard data from detection histories
containing fewer than \code{min.recapt}+1 detections.

}

\value{

Scalar distance in metres, or a list of such values if \code{capthist}
is a multi-session list.

The \code{full} argument may be used with \code{MMDM} and \code{ARL} to
return more extensive output, particularly the observed range length for
each detection history.

}

\references{

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Jett, D. A. and Nichols, J. D. (1987) A field comparison of nested grid
and trapping web density estimators. \emph{Journal of Mammalogy}
\bold{68}, 888--892.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\note{

All measures are affected by the arrangement of detectors. \code{dbar}
is also affected quite strongly by serial correlation in the sampled
locations. Using \code{dbar} with 'proximity' detectors raises a problem
of interpretation, as the original sequence of multiple detections
within an occasion is unknown. RPSV is a value analogous to the standard
deviation of locations about the home range centre.

The value returned by \code{dbar} for 'proximity' or 'count' detectors
is of little use because multiple detections of an individual within an
occasion are in arbitrary order.

Inclusion of these measures in the \pkg{secr} package does not mean they are
recommended for general use! It is usually better to use a spatial
parameter from a fitted model (e.g., \eqn{\sigma}{sigma} of the
half-normal detection function). Even then, be careful that
\eqn{\sigma}{sigma} is not 'contaminated' with behavioural effects (e.g.
attraction of animal to detector) or 'detection at a distance'.

}

\seealso{ \code{\link{autoini}} }

\examples{
data(secrdemo)
dbar(captdata)
RPSV(captdata)
}

\keyword{ models }

