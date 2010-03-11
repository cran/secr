\name{secr.fit}
\alias{secr.fit}

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
    mask = NULL, buffer = 100, CL = FALSE, detectfn = NULL,
    start = NULL, link = list(), fixed = list(),
    timecov = NULL, sessioncov = NULL, groups = NULL, dframe = NULL,
    details = list(), method = 'Newton-Raphson', verify = TRUE, trace = NULL, ...)
}

\arguments{
  \item{capthist}{ \code{\link{capthist}} object including capture data and detector (trap) layout }
  \item{mask}{ \code{\link{mask}} object}
  \item{buffer}{ scalar mask buffer radius if \code{mask} not specified}
  \item{CL}{ logical, if true then the model is fitted by maximizing the conditional likelihood }
  \item{detectfn}{ integer code for shape of detection function. 0 halfnormal, 1 hazard-rate, 2 exponential etc.}
  \item{start}{ vector of initial values for beta parameters, or \code{secr} object from which they may be derived }
  \item{link}{ list with optional components 'D', 'g0', 'sigma' and 'z', each a character string in \{'log', 'logit', 'identity', 'sin'\} for the link function of the relevant real parameter }
  \item{fixed}{ list with optional components corresponding to each 'real' parameter (e.g., 'D', 'g0', 'sigma'), the scalar value to which parameter is to be fixed }
  \item{model}{ list with optional components 'D', 'g0', 'sigma' and 'z', each symbolically defining a linear predictor for the relevant real parameter using \code{formula} notation }
  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
  \item{groups}{ optional vector of one or more variables with which to form groups. Each element should be the name of a factor variable in the \code{covariates} attribute of \code{capthist}. }
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{details}{ list of additional settings, mostly model-specific (see Details) }
  \item{method}{ character string giving method for maximizing log likelihood }
  \item{verify}{ logical, if TRUE the input data are checked with \code{\link{verify}} }
  \item{trace}{ logical, if TRUE then output each evaluation of the likelihood, and other messages}
  \item{\dots}{ other arguments passed to the maximization function }
}

\details{

\code{secr.fit} fits a SECR model by maximizing the likelihood. The
likelihood depends on the detector type ('multi' or 'proximity') of the
'\code{traps}' attribute of \code{capthist} (Borchers and Efford 2008,
Efford, Borchers and Byrom 2009, Efford, Dawson and Borchers 2009). The
'multi' form of the likelihood is also used, with a warning, when
detector type = 'single' (see Efford et al. 2009). Default \code{model}
is null (constant density and detection probability). The set of
variables available for use in linear predictors includes some that are
constructed automatically (t, T, b, B), group (g), and others that appear
in the \code{covariates} of the input data. See \link{secr models} for
more on defining models.

The length of \code{timecov} should equal the number of sampling
occasions (\code{ncol(capthist)}). Arguments \code{timecov},
\code{sessioncov} and \code{groups} are used only when needed for terms
in one of the model specifications. Default \code{link} is list(D='log',
g0='logit', sigma='log').

If \code{start} is missing then \code{\link{autoini}} is used for D, g0
and sigma, and other beta parameters are set initially to arbitrary
values, mostly zero. \code{start} may be a previously fitted nested
model. In this case, a vector of starting beta values is constructed
from the nested model and additional betas are set to zero. Mapping of
parameters follows the default in \code{\link{score.test}}, but user
intervention is not allowed.

\code{details} is used for various specialized settings --

\code{details$distribution} specifies the distribution of the number of
individuals detected; this may be conditional on the number in the
masked area ('binomial') or unconditional ('poisson').
\code{distribution} affects the sampling variance of the estimated
density. The default is 'poisson'.

\code{details$hessian} is a character string controlling the computation
of the Hessian matrix from which variances and covariances are obtained.
Options are 'none' (no variances), 'auto' (the default) or 'fdhess' (use
the function fdHess in \pkg{nlme}).  If 'auto' then the Hessian from the
optimisation function is used.

\code{details$LLonly} = TRUE causes the function to returns a single
evaluation of the log likelihood at the 'start' values.

\code{details$scalesigma} = TRUE causes sigma to be scaled by
\eqn{\mathrm{D}^{-0.5}}{1/sqrt(D)}.

\code{details$scaleg0} = TRUE causes g0 to be scaled by
\eqn{\mathrm{sigma}^{-2}}{sigma^-2}. The corresponding 'real' parameters
are marked with an asterisk in output (e.g. g0*).

If \code{method} = 'Newton-Raphson' then \code{\link[stats]{nlm}} is
used to maximize the log likelihood; otherwise
\code{\link[stats]{optim}} is used with the chosen method
('BFGS','Nelder-Mead', etc.).  A feature of \code{nlm} is that it takes a
large step early on in the maximization that may cause floating point
underflow or overflow in one or more real parameters.  This can be
controlled by passing the 'stepmax' argument of \code{nlm} in the \dots
argument of \code{secr.fit} (see first example).

If \code{verify} = TRUE then \code{\link{verify}} is called to check
capthist and mask; analysis is aborted if errors are found.

}

\value{
The function \code{secr.fit} returns an object of class secr.  This has components
  \item{call }{function call}
  \item{capthist }{saved input}
  \item{mask }{saved input}
  \item{detectfn }{saved input}
  \item{CL }{saved input}
  \item{timecov }{saved input}
  \item{sessioncov }{saved input}
  \item{groups }{saved input}
  \item{dframe }{saved input}
  \item{design }{reduced design matrices, parameter table and parameter index array for actual animals (see \code{\link{secr.design.MS}})}
  \item{design0 }{reduced design matrices, parameter table and parameter index array for 'naive' animal (see \code{\link{secr.design.MS}})}
  \item{start }{vector of starting values for beta parameters}
  \item{link }{list with components for each real parameter (e.g., 'D', 'g0'), the name of the link function used for each real parameter. Component 'z' is NULL unless detectfn = 1 (hazard-rate). }
  \item{fixed }{saved input}
  \item{parindx }{list with possible components 'D', 'g0', 'sigma' and 'z', for the indices of the 'beta' parameters associated with each real parameter ('z' NULL unless detectfn = 1).}
  \item{model }{saved input}
  \item{details }{saved input}
  \item{vars }{vector of unique variable names in \code{model} }
  \item{betanames }{names of beta parameters}
  \item{realnames }{names of fitted (real) parameters }
  \item{fit }{list describing the fit (output from \code{nlm} or \code{optim}) }
  \item{beta.vcv }{variance-covariance matrix of beta parameters }
  \item{D }{ array of predicted densities of each group at each mask point in each session, dim(D) = c(nrow(mask), ngroups, nsessions) }
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

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture--recapture: likelihood-based methods. In:
D. L. Thompson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer. Pp. 255--269.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\note{

\code{print}, \code{AIC}, \code{vcov}, and \code{predict} methods
are provided. \code{derived} is used to compute the derived parameters
'esa' (effective sampling area) and 'D' (density) for models fitted by
maximizing the conditional likelihood (CL = TRUE).

Components 'version' and 'starttime' were introduced in version 1.2.7,
and recording of the completion time in 'fitted' was discontinued.

}

\seealso{
\code{\link{capthist}}, \code{\link{mask}}, \link{detection functions},
\code{\link{print.secr}}, \code{\link{vcov.secr}}, \code{\link{AIC.secr}},
\code{\link{derived}}, \code{\link{predict.secr}}, \code{\link{verify}}
}

\examples{

## construct test data (array of 48 'multi-catch' traps)

detectors <- make.grid (nx = 6, ny = 8, detector = 'multi')
detections <- sim.capthist (detectors, popn = list(D = 10,
    buffer = 100), detectpar = list(g0 = 0.2, sigma = 25))

## fit & print null (constant parameter) model
## stepmax is passed to nlm (not usually needed)

secr0 <- secr.fit (detections, stepmax = 50)
secr0   ## uses print method for secr

## compare fit of null model with learned-response model for g0

secrb <- secr.fit (detections, model = g0~b)
AIC (secr0, secrb)

## typical result

##                  model   detectfn npar    logLik     AIC    AICc dAICc  AICwt
## secr0 D~1 g0~1 sigma~1 halfnormal    3 -347.1210 700.242 700.928 0.000 0.7733
## secrb D~1 g0~b sigma~1 halfnormal    4 -347.1026 702.205 703.382 2.454 0.2267

}

\keyword{ models }

