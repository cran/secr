\name{sim.secr}
\Rdversion{1.1}
\alias{simulate.secr}
\alias{sim.secr}
\alias{simulate}

\title{
  Simulate From Fitted secr Model
}

\description{

Simulate a spatially distributed population, sample from that population
with an array of detectors, and optionally fit an SECR model to the
simulated data.

}

\usage{

\method{simulate}{secr}(object, nsim = 1, seed = NULL, chat = 1, ...)

sim.secr(object, nsim = 1, extractfn = function(x) c(deviance =
    deviance(x), df = df.residual(x)), seed = NULL, data = NULL,
    tracelevel = 1, hessian = "none", start = object$fit$par)

}

\arguments{
  \item{object}{
an secr object
}
  \item{nsim}{
number of replicates
}
  \item{seed}{
value for setting .Random.seed - either NULL or an integer
}
  \item{chat}{
real value for overdispersion parameter
}
  \item{extractfn}{
function to extract output values from fitted model
}
  \item{data}{
optional list of simulated data saved from previous call to \code{simulate.secr}
}
  \item{tracelevel}{
integer for level of detail in reporting (0,1,2)
}
  \item{hessian}{
character string controlling the computation of the Hessian matrix
}
  \item{start}{
vector of starting 'beta' values for \code{secr.fit}
}
  \item{...}{
other arguments (not used)
}
}
\details{

For each replicate, \code{simulate.secr} calls \code{\link{sim.popn}} to
generate session- and group-specific realizations of the (possibly
inhomogeneous) 2-D Poisson distribution fitted in \code{object}, across
the habitat mask(s) in \code{object}.  Group subpopulations are combined
using \code{\link{rbind.popn}} within each session; information to
reconstruct groups is retained in the individual-level factor
covariate(s) of the resulting \code{popn} object (corresponding to
\code{object\$groups}).  Each population is then sampled using the fitted
detection model and detector (trap) array(s) in \code{object}.

The random number seed is managed as in \code{simulate.lm}.

\code{simulate.secr} does not yet work with models fitted using
conditional likelihood (\code{object$CL = TRUE}).  Detector type is
determined by \code{detector(traps(object$capthist))}, which should be one of
'single', 'multi', 'proximity', 'areasearch' or 'count'.

\code{sim.secr} is a wrapper function. If \code{data = NULL} (the
default) then it calls \code{simulate.secr} to generate new datasets. If
\code{data} is provided then \code{nsim} is taken to be
\code{length(data)}. \code{secr.fit} is called to fit the original model
to each new dataset. Results are summarized according to the
user-provided function \code{extractfn}. The default \code{extractfn}
returns the deviance and its degrees of freedom; a NULL value for
\code{extractfn} returns the fitted secr objects after
\code{\link{trim}}ming to reduce bulk. Simulation uses the detector type
of the data, even when another likelihood is fitted (this is the case
with single-catch data, for which a multi-catch likelihood is fitted).
Warning messages from \code{secr.fit} are suppressed.

\code{extractfn} should be a function that takes an \code{secr} object
as its only argument.

\code{tracelevel=0} suppresses most messages; \code{tracelevel=1} gives a
terse message at the start of each fit; \code{tracelevel=2} also sets
'details\$trace = TRUE' for \code{secr.fit}, causing each likelihood
evaluation to be reported.

It is OK (and faster) to use \code{hessian='none'} unless
\code{extractfn} needs variances or covariances.

\code{sim.capthist} is a more direct way to simulate data from a null
model (i.e. one with constant parameters for density and detection). It
is limited to a single session.

}
\value{

For \code{simulate.secr}, a list of data sets ('capthist' objects). This
list has class=('list','secrdata'); the initial state of the random
number generator (roughly, the value of .Random.seed) is stored as the
attribute 'seed'.

The value from \code{sim.secr} depends on \code{extractfn}: if that
returns a numeric vector of length \code{n.extract} then the value is a
matrix with \code{dim = c(nsim, n.extract)} (i.e., the matrix has one
row per replicate and one column for each extracted value). Otherwise,
the value returned by \code{sim.secr} is a list with one component per
replicate (strictly, an object of class = c('list','secrlist')). Each
simulated fit may be retrieved \emph{in toto} by specifying
\code{extractfn = identity}, or slimmed down by specifying
\code{extractfn = NULL} or \code{extractfn = trim}, which are
equivalent.

For either form of output from \code{sim.secr} the initial state of the
random number generator is stored as the attribute 'seed'.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{
\code{\link{sim.capthist}}, \code{\link{secr.fit}}, \code{\link[stats]{simulate}}
}

\note{

The value returned by \code{simulate.secr} is a list of 'capthist'
objects; if there is more than one session, each 'capthist' is itself a
sort of list .

The classes 'secrdata' and 'secrlist' are used only to override
the ugly and usually unwanted printing of the seed attribute.

The default value for \code{start} in \code{sim.secr} is the previously
fitted parameter vector. Alternatives are NULL or \code{object$start}.

}

\examples{

data(secrdemo)
simulate(secrdemo.0, nsim = 2)

\dontrun{

## this would take a long time...
sims <- sim.secr(secrdemo.0, nsim = 99)
deviance(secrdemo.0)
devs <- c(deviance(secrdemo.0),sims$deviance)
quantile(devs, probs=c(0.95))
rank(devs)[1] / length(devs)

## to assess bias
extrfn <- function (object) unlist(predict(object)['D',-1])
sims <- sim.secr(secrdemo.0, nsim = 50, extractfn=extrfn)
sims

## with a larger sample, could get parametric bootstrap CI
quantile(sims[,1], c(0.025, 0.975))

}

}
\keyword{ models }

