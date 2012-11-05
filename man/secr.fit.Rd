\name{secr.fit}
\alias{secr.fit}
\alias{binomN}

\title{ Spatially Explicit Capture--Recapture }

\description{

Estimate animal population density with data from an array of passive
detectors (traps) by fitting a spatial detection model by maximizing the
likelihood. Data must have been assembled as an object of class
\code{capthist}. Integration is by summation over the grid of points in
\code{mask}.

}

\usage{
secr.fit (capthist, model = list(D~1, g0~1, sigma~1),
    mask = NULL, buffer = NULL, CL = FALSE, detectfn = NULL,
    binomN = NULL, start = NULL, link = list(), fixed = list(),
    timecov = NULL, sessioncov = NULL, groups = NULL,
    dframe = NULL, details = list(), method = "Newton-Raphson",
    verify = TRUE, biasLimit = 0.01, trace = NULL, ncores = 1, ...)
}

\arguments{
  \item{capthist}{ \code{\link{capthist}} object including capture data and detector (trap) layout }
  \item{mask}{ \code{\link{mask}} object}
  \item{buffer}{ scalar mask buffer radius if \code{mask} not specified
  (default 100 m)}
  \item{CL}{ logical, if true then the model is fitted by maximizing the conditional likelihood }
  \item{detectfn}{integer code or character string for shape of
  detection function 0 = halfnormal, 1 = hazard rate etc. -- see
  \link{detectfn}}
  \item{binomN}{ integer code for distribution of counts (see Details) }
  \item{start}{ vector of initial values for beta parameters, or \code{secr} object from which they may be derived }
  \item{link}{ list with optional components `D', `g0', `sigma' and `z', each a character string in \{"log", "logit", "identity", "sin"\} for the link function of the relevant real parameter }
  \item{fixed}{ list with optional components corresponding to each `real' parameter (e.g., `D', `g0', `sigma'), the scalar value to which parameter is to be fixed }
  \item{model}{ list with optional components `D', `g0', `sigma' and `z', each symbolically defining a linear predictor for the relevant real parameter using \code{formula} notation }
  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
  \item{groups}{ optional vector of one or more variables with which to form groups. Each element should be the name of a factor variable in the \code{covariates} attribute of \code{capthist}. }
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{details}{ list of additional settings, mostly model-specific (see Details) }
  \item{method}{ character string giving method for maximizing log likelihood }
  \item{verify}{ logical, if TRUE the input data are checked with \code{\link{verify}} }
  \item{biasLimit}{ numeric threshold for predicted relative bias due to
  buffer being too small }
  \item{trace}{ logical, if TRUE then output each evaluation of the likelihood, and other messages}
  \item{ncores}{ integer number of cores available for parallel processing}
  \item{\dots}{ other arguments passed to the maximization function }
}

\details{

\code{secr.fit} fits a SECR model by maximizing the likelihood. The
likelihood depends on the detector type ("multi", "proximity", "count",
"polygon" etc.) of the \code{traps} attribute of \code{capthist}
(Borchers and Efford 2008, Efford, Borchers and Byrom 2009, Efford,
Dawson and Borchers 2009, Efford 2011). The `multi' form of the
likelihood is also used, with a warning, when detector type = "single"
(see Efford et al. 2009 for justification). The default \code{model} is
null (constant density and detection probability). The set of variables
available for use in linear predictors includes some that are
constructed automatically (t, T, b, B, bk, Bk, k, K), group (g), and
others that appear in the \code{covariates} of the input data. See also
\code{\link{timevaryingcov}} to construct time-varying detector
covariates, and \link{secr models} and \url{../doc/secr-overview.pdf} for
more on defining models.

\code{buffer} and \code{mask} are alternative ways to define the region
of integration (see \link{mask}).

The length of \code{timecov} should equal the number of sampling
occasions (\code{ncol(capthist)}). Arguments \code{timecov},
\code{sessioncov} and \code{groups} are used only when needed for terms
in one of the model specifications. Default \code{link} is \code{list(D="log",
g0="logit", sigma="log")}.

If \code{start} is missing then \code{\link{autoini}} is used for D, g0
and sigma, and other beta parameters are set initially to arbitrary
values, mostly zero. \code{start} may be a previously fitted nested
model. In this case, a vector of starting beta values is constructed
from the nested model and additional betas are set to zero. Mapping of
parameters follows the default in \code{\link{score.test}}, but user
intervention is not allowed.

\code{binomN} (previously a component of \code{details}) determines the
distribution that is fitted for the number of detections of an individual
at a particular detector, on a particular occasion, when the detectors
are of type `count', `polygon' or `transect':

\itemize{
\item{binomN > 1} {binomial with size binomN}
\item{binomN = 1} {Bernoulli}
\item{binomN = 0} {Poisson}
\item{binomN < 0} {negative binomial with size abs(binomN) -- see 
  \code{\link[stats]{dnbinom}}}
}

The default with these detectors is to fit a Poisson distribution. The
`size' parameter of the negative binomial is not estimated: it must be
supplied. \code{binomN} should be an integer unless negative.

\code{details} is used for various specialized settings --

\code{details$distribution} specifies the distribution of the number of
individuals detected; this may be conditional on the number in the
masked area ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The default is "poisson". See also Note.

\code{details$hessian} is a character string controlling the computation
of the Hessian matrix from which variances and covariances are obtained.
Options are "none" (no variances), "auto" (the default) or "fdhess" (use
the function fdHess in \pkg{nlme}).  If "auto" then the Hessian from the
optimisation function is used.

\code{details$LLonly} = TRUE causes the function to returns a single
evaluation of the log likelihood at the `start' values.

\code{details$scalesigma} = TRUE causes sigma to be scaled by
\eqn{\mathrm{D}^{-0.5}}{1/sqrt(D)}.

\code{details$scaleg0} = TRUE causes g0 to be scaled by
\eqn{\mathrm{sigma}^{-2}}{sigma^-2}.

\code{details$centred} = TRUE causes coordinates of both traps and mask
to be centred on the centroid of the traps, computed separately for each
session in the case of multi-session data. This may be necessary to
overcome numerical problems when x- or y-coordinates are large
numbers. The default is not to centre coordinates.

\code{details$param} = 1 causes the Gardner & Royle parameterisation of
the detection model (\emph{p0}, \eqn{\sigma}; Gardner et al. 2009) to be
used for multi-catch detectors (default 0 for Borchers and Efford). This
parameterisation does not allow detector covariates.

\code{telemetrysigma} = TRUE uses coordinate information from telemetry
when capthist has attribute `xylist' (see \code{\link{addTelemetry}}).

If \code{method = "Newton-Raphson"} then \code{\link[stats]{nlm}} is
used to maximize the log likelihood (minimize the negative log
likelihood); otherwise \code{\link[stats]{optim}} is used with the
chosen method ("BFGS", "Nelder-Mead", etc.).  If maximization fails a
warning is given appropriate to the method.

If \code{verify} = TRUE then \code{\link{verify}} is called to check
capthist and mask; analysis is aborted if "errors" are found. Some
conditions that trigger an "error" are benign (e.g., no detections in
some sessions of a multi-session study of a sparse population); use
\code{verify = FALSE} to avoid the check. See also Note.

If \code{buffer} is used rather than \code{mask}, and \code{biasLimit}
    is valid, then the estimated density is checked for bias due to the
    choice of buffer. A warning is generated when \code{buffer} appears
    to be too small (predicted RB(D-hat) > \code{biasLimit}, default 1\%
    relative bias). The prediction uses \code{\link{bias.D}}. No check
    is performed when \code{mask} is specified, when \code{biasLimit} is 0,
    negative or NA, or when the detector type is "polygon", "transect",
    "polygonX" or "transectX".

If \code{ncores > 1} the \pkg{parallel} package will be used to create
processes on multiple cores (see \link{Parallel} for more). Specifying
extra cores may improve the speed of multi-session analyses (it may also
slow them down, as data must be copied back and forth). There is
presently no benefit for single-session analyses.

}

\value{
The function \code{secr.fit} returns an object of class secr.  This has
components

  \item{call }{function call (as character string prior to secr 1.5) }

  \item{capthist }{saved input}
  \item{mask }{saved input}
  \item{detectfn }{saved input}
  \item{CL }{saved input}
  \item{timecov }{saved input}
  \item{sessioncov }{saved input}
  \item{groups }{saved input}
  \item{dframe }{saved input}
  
  \item{design }{reduced design matrices, parameter table and parameter
    index array for actual animals (see \code{\link{secr.design.MS}})}
  
  \item{design0 }{reduced design matrices, parameter table and parameter
    index array for `naive' animal (see \code{\link{secr.design.MS}})}

  \item{start }{vector of starting values for beta parameters} 

  \item{link }{list with one component for each real parameter
  (typically `D', `g0', `sigma'),giving the name of the link function
  used for each real parameter.}

\item{fixed }{saved input}  

  \item{parindx }{list with one component for each real parameter giving
    the indices of the `beta' parameters associated with each real
    parameter}  

  \item{model }{saved input}

  \item{details }{saved input}

  \item{vars }{vector of unique variable names in \code{model} }

  \item{betanames }{names of beta parameters}

  \item{realnames }{names of fitted (real) parameters }

  \item{fit }{list describing the fit (output from \code{nlm} or
  \code{optim}) }

  \item{beta.vcv }{variance-covariance matrix of beta parameters }  

  \item{N }{ if CL = FALSE, array of predicted number in each group at
    in each session, summed across mask, dim(N) = c(ngroups, nsessions),
    otherwise NULL }

  \item{version }{secr version number }
  \item{starttime }{character string of date and time at start of fit }
  \item{proctime }{processor time for model fit, in seconds }
}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum
likelihood methods for capture--recapture studies. \emph{Biometrics}
\bold{64}, 377--385.

Efford, M. G. (2004) Density estimation in live-trapping studies.
\emph{Oikos} \bold{106}, 598--610.

Efford, M. G. (2011) Estimation of population density by spatially
explicit capture--recapture with area or transect searches. Unpublished
manuscript.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

Gardner, B., Royle, J. A. and Wegan, M. T. (2009) Hierarchical models
for estimating density from DNA mark-recapture studies. \emph{Ecology}
\bold{90}, 1106--1115.

}

\note{

One system of units is used throughout \pkg{secr}. Distances are in metres and
areas are in hectares (ha). The unit of density is animals per
hectare. 1 ha = 10000 m^2 = 0.01 km^2. To convert density to animals /
km^2, multiply by 100.

\code{print}, \code{AIC}, \code{vcov}, and \code{predict} methods
are provided. \code{derived} is used to compute the derived parameters
`esa' (effective sampling area) and `D' (density) for models fitted by
maximizing the conditional likelihood (CL = TRUE).

The component `distribution' of argument `details' may also take a
    numeric value larger than nrow(capthist), rather than "binomial" or
    "poisson". The likelihood then treats n as a binomial draw from a
    superpopulation of this size, with consequences for the variance of
    density estimates. This can help to reconcile MLE with Bayesian
    estimates using data augmentation.
    
Components `version' and `starttime' were introduced in version 1.2.7,
and recording of the completion time in `fitted' was discontinued.

 The Newton-Raphson algorithm is fast, but it sometimes fails to compute
the information matrix correctly, causing some or all standard errors to
be set to NA. This usually indicates a major problem in fitting the
model, and parameter estimates should not be trusted. The alternative
\code{method = "BFGS"} often works better in these cases, or use
\code{details = list(hessian = "fdhess")}.

The component D in output was replaced with N from version 2.3. Use
\code{\link{region.N}} to obtain SE or confidence intervals for N-hat,
or to infer N for a different region.

Prior to version 2.3.2 the buffer bias check could be switched off by
setting \code{verify = FALSE}. This is now done by setting
\code{biasLimit = 0}.

}

\seealso{

  \link{detection functions},
  \code{\link{AIC.secr}},
  \code{\link{capthist}},
  \code{\link{derived}},
  \code{\link{mask}},
  \code{\link{predict.secr}},
  \code{\link{print.secr}},
  \code{\link{region.N}},
  \code{\link{troubleshooting}}
  \code{\link{vcov.secr}},
  \code{\link{verify}},

}

\examples{

\dontrun{

## construct test data (array of 48 `multi-catch' traps)

detectors <- make.grid (nx = 6, ny = 8, detector = "multi")
detections <- sim.capthist (detectors, popn = list(D = 10,
    buffer = 100), detectpar = list(g0 = 0.2, sigma = 25))

## fit & print null (constant parameter) model
secr0 <- secr.fit (detections)
secr0   ## uses print method for secr

## compare fit of null model with learned-response model for g0

secrb <- secr.fit (detections, model = g0~b)
AIC (secr0, secrb)

## typical result

##                  model   detectfn npar    logLik     AIC    AICc dAICc  AICwt
## secr0 D~1 g0~1 sigma~1 halfnormal    3 -347.1210 700.242 700.928 0.000 0.7733
## secrb D~1 g0~b sigma~1 halfnormal    4 -347.1026 702.205 703.382 2.454 0.2267
}
}

\keyword{ models }

