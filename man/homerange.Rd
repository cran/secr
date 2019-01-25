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

\code{MMDM} (for `Mean Maximum Distance Moved') is the average maximum
distance between detections of each individual i.e. the observed range
length averaged over individuals (Otis et al. 1978).

\code{ARL} (or `Asymptotic Range Length') is obtained by fitting an
exponential curve to the scatter of observed individual range length vs
the number of detections of each individual (Jett and Nichols 1987: 889).

\code{RPSV} (for `Root Pooled Spatial Variance') is a measure of the 2-D
dispersion of the locations at which individual animals are detected,
pooled over individuals (cf Calhoun and Casby 1958, Slade and Swihart 1983).

\code{moves} reports the distance between successive detections of each animal.

}

\usage{
dbar(capthist, userdist = NULL, mask = NULL)
MMDM(capthist, min.recapt = 1, full = FALSE, userdist = NULL, mask = NULL)
ARL(capthist, min.recapt = 1, plt = FALSE, full = FALSE, userdist = NULL, mask = NULL)
RPSV(capthist, CC = FALSE)
moves(capthist, userdist = NULL, mask = NULL, names = FALSE)
}

\arguments{
  \item{capthist}{ object of class \code{\link{capthist}} }
  \item{userdist}{ function or matrix with user-defined distances}
  \item{mask}{ habitat mask passed to userdist function, if required}
  \item{names}{ logical; should results be ordered alphanumerically by row names? }
  \item{min.recapt}{ integer minimum number of recaptures for a
  detection history to be used }
  \item{plt}{ logical; if TRUE observed range length is plotted against number
  of recaptures }
\item{full}{ logical; set to TRUE for detailed output }
\item{CC}{ logical for whether to use Calhoun and Casby formula }
}

\details{
 \code{dbar} is defined as --
 \deqn{
\overline{d}=\frac{\sum\limits _{i=1}^{n}
                          \sum\limits _{j=1}^{n_i - 1}
                              \sqrt{(x_{i,j}-x_{i,j+1})^2 + (y_{i,j}-y_{i,j+1})^2}}
                         {\sum\limits _{i=1}^{n} (n_i-1)}}{see pdf manual}

When \code{CC = FALSE}, \code{RPSV} is defined as --
\deqn{
 RPSV = \sqrt{
    \frac {\sum\limits _{i=1}^{n} \sum\limits _{j=1}^{n_i} [
	    (x_{i,j} - \overline x_i)^2 + (y_{i,j} - \overline y_i)^2
		]}{\sum\limits _{i=1}^{n} (n_i-1) - 1}}
	  }{see pdf manual}.

Otherwise (\code{CC = TRUE}), \code{RPSV} uses the formula of Calhoun
and Casby (1958) with a different denominator -- 	  
\deqn{
 s = \sqrt{
    \frac {\sum\limits _{i=1}^{n} \sum\limits _{j=1}^{n_i} [
	    (x_{i,j} - \overline x_i)^2 + (y_{i,j} - \overline y_i)^2
		]}{2\sum\limits _{i=1}^{n} (n_i-1)}}
	  }{see pdf manual}.

	  The Calhoun and Casby formula (offered from 2.9.1) correctly estimates \eqn{\sigma}
	  when trapping is on an infinite, fine grid, and is preferred
	  for this reason. The original RPSV
	  (\code{CC = FALSE}) is retained as the default for compatibility with
	  previous versions of \pkg{secr}.
	  
\code{dbar} and \code{RPSV} have a specific role as proxies for
detection scale in inverse-prediction estimation of density (Efford
2004; see \code{\link{ip.secr}}).

\code{RPSV} is used in \code{autoini} to obtain plausible starting
values for maximum likelihood estimation.

\code{MMDM} and \code{ARL} discard data from detection histories
containing fewer than \code{min.recapt}+1 detections.

The \code{userdist} option is included for exotic non-Euclidean cases
(see e.g. \code{secr.fit} \link{details}). RPSV is not defined for
non-Euclidean distances.

If \code{capthist} comprises standalone telemetry data (all detector 'telemetry') 
then calculations are performed on the telemetry coordinates.

Movements are reliably reported by \code{moves} only if there is a maximum of one detection per animal per occasion. The sequence of detections within any occasion is not known; where these occur the sequence used by \code{moves} is arbitrary (sequence follows detector index).

}

\value{

For \code{dbar}, \code{MMDM}, \code{ARL} and \code{RPSV} --

Scalar distance in metres, or a list of such values if \code{capthist}
is a multi-session list.

The \code{full} argument may be used with \code{MMDM} and \code{ARL} to
return more extensive output, particularly the observed range length for
each detection history.

For \code{moves} --

List with one component for each animal, a vector of distances, or numeric(0) if the animal is detected only once. A list of such lists if \code{capthist} is a multi-session list.

}

\references{

Calhoun, J. B. and Casby, J. U. (1958) Calculation of home range and
  density of small mammals. Public Health
  Monograph. No. 55. U.S. Government Printing Office.
  
Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Jett, D. A. and Nichols, J. D. (1987) A field comparison of nested grid
and trapping web density estimators. \emph{Journal of Mammalogy}
\bold{68}, 888--892.

Otis, D. L., Burnham, K. P., White, G. C. and Anderson, D. R. (1978)
  Statistical inference from capture data on closed animal
  populations. \emph{Wildlife Monographs} \bold{62}, 1--135.

Slade, N. A. and Swihart, R. K. (1983) Home range indices for the hispid
cotton rat (Sigmodon hispidus) in Northeastern Kansas. \emph{Journal of
Mammalogy} \bold{64}, 580--590.

}

\note{

All measures are affected by the arrangement of detectors. \code{dbar}
is also affected quite strongly by serial correlation in the sampled
locations. Using \code{dbar} with `proximity' detectors raises a problem
of interpretation, as the original sequence of multiple detections
within an occasion is unknown. RPSV is a value analogous to the standard
deviation of locations about the home range centre.

The value returned by \code{dbar} for `proximity' or `count' detectors
is of little use because multiple detections of an individual within an
occasion are in arbitrary order.

Inclusion of these measures in the \pkg{secr} package does not mean they are
recommended for general use! It is usually better to use a spatial
parameter from a fitted model (e.g., \eqn{\sigma}{sigma} of the
half-normal detection function). Even then, be careful that
\eqn{\sigma}{sigma} is not `contaminated' with behavioural effects (e.g.
attraction of animal to detector) or `detection at a distance'.

The argument 'names' was added in 3.0.1. The default \code{names = FALSE} causes a change in behaviour from that version onwards.

}

\seealso{ \code{\link{autoini}} }

\examples{
dbar(captdata)
RPSV(captdata)
RPSV(captdata, CC = TRUE)
}

\keyword{ models }

