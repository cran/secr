\name{sim.capthist}
\alias{sim.capthist}
\alias{sim.resight}

\title{ Simulate Detection Histories }

\description{

Create a set of capture or marking-and-resighting histories by simulated sampling 
of a 2-D population using an array of detectors.

}

\usage{
sim.capthist(traps, popn = list(D = 5, buffer = 100,
    Ndist = 'poisson'), detectfn = 0, detectpar = list(), 
    noccasions = 5, renumber = TRUE, seed = NULL)
sim.resight(traps, ..., q = 1, pID = 1, unmarked = TRUE,
    nonID = TRUE)

}

\arguments{
  \item{traps}{ \code{traps} object with the locations and other attributes of detectors }
  \item{popn}{ locations of individuals in the population to be sampled, either as
 a \code{\link{popn}} object or a list with named elements 'D' (density) and 'buffer' }
  \item{detectfn}{ code for type of detection function }
  \item{detectpar}{ list of values for named parameters of detection function } 
  \item{noccasions}{ number of occasions to simulate }
  \item{renumber}{ logical for whether output rows should labeled sequentially (TRUE) or
 retain the numbering of the population from which they were drawn (FALSE) }
  \item{seed}{an object specifying if and how the random number
  generator should be initialized ('seeded') }
  \item{\dots}{ arguments to pass to \code{sim.capthist} }
  \item{q}{ number of marking occasions }
  \item{pID}{ probability of individual identification for marked animals }
  \item{unmarked}{ logical, if true unmarked individuals are not recorded during 'sighting' }
  \item{nonID}{ logical, if true then unidentified marked individuals are not recorded during 'sighting' }

}
\details{

If \code{popn} is not of class 'popn' then a homogeneous Poisson
population with the desired density (animals/ha) is first simulated over
the rectangular area of the bounding box of \code{traps} plus a buffer
of the requested width (metres). The detection algorithm depends on the
detector type of \code{traps}. For 'proximity' detectors, the actual
detection probability of animal \emph{i} at detector \emph{j} is the
naive probability given by the detection function. For 'single' and
'multi' detectors the naive probability is modified by competition
between detectors and, in the case of 'single' detectors, between animals. See
Efford (2004) and other papers below for details. 

Detection parameters are specific to the detection function, which is
indicated by a numeric code (\code{\link{detection functions}}).
Parameters may vary with time - for this provide a vector of length
\code{noccasion}.  The default detection parameters are \code{list(g0 = 0.2,
 sigma = 25, z = 1)}.

\code{detectpar} may also include 'binomN' and other arguments for
detectors that have yet to be documented. A zero value for \code{binomN}
indicates that counts should be modelled with a Poisson distribution.

If \code{popn} is specified by an object of class 'popn' then any
individual covariates will be passed on; the \code{covariates} attribute
of the output is otherwise set to NULL.

The random number seed is managed as in \code{\link[stats]{simulate}}.

\code{sim.resight} generates mark-resight data for 'q' marking occasions
followed by 'noccasion -- q' sighting occasions. \code{sim.capthist} is
first called with the arguments 'traps' and \dots. The detector type
must be 'proximity'. The 'usage' attribute of \code{traps} is ignored at
present, so the same detectors are operated on all occasions.  Any
detection-parameter vector of length 2 in \dots is interpreted as
providing differing constant values for the marking and sighting phases.

}
\value{

For \code{sim.capthist}, an object of class \code{capthist}, a matrix or
3-dimensional array with additional attributes. Rows represent
individuals and columns represent occasions; the third dimension, used
when detector type = 'proximity', codes presence or absence at each
detector. For trap detectors ('single', 'multi') each entry in
\code{capthist} is either zero (no detection) or the sequence number of
the trap.

The initial state of the R random number generator is stored in the
'seed' attribute.

For \code{sim.resight}, an object of class \code{capthist}, always a
3-dimensional array, with additional attributes Tu and Tm containing
counts of 'unmarked' and 'marked, not identified' sightings.

}

\references{ 

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\note{ 

External code is called to speed the simulations. The present version
assumes a 'null model' i.e. naive detection probability is constant
except for effects of distance and possibly time (using vector-valued
detection parameters from 1.2.10). You can, however, use
\code{\link{rbind.capthist}} to combine detections of population
subclasses (e.g. males and females) simulated with different parameter
values. This is not valid for detector type 'single' because it fails to
allow for competition for traps between subclasses. Future versions may
allow more complex models. 

}

\seealso{ 
    \code{\link{sim.popn}}, \code{\link{capthist}}, \code{\link{traps}}, 
    \code{\link{popn}}, \code{\link{detection functions}}, \code{\link[stats]{simulate}}
}

\examples{
## simple example
## detector = 'multi' (default)
temptrap <- make.grid(nx = 6, ny = 6, spacing = 20)
sim.capthist (temptrap, detectpar = list(g0 = 0.2, sigma = 20))

## with detector = 'proximity, there may be more than one 
## detection per individual per occasion
temptrap <- make.grid(nx = 6, ny = 6, spacing = 20, detector = 
    'proximity')
summary(sim.capthist (temptrap, detectpar = 
    list(g0 = 0.2, sigma = 20)))

}
\keyword{ datagen }

