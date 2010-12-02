\name{secr.model}
\alias{secr.model}
\alias{secr models}

\title{ Spatially Explicit Capture--Recapture Models }

\description{

A family of capture--recapture models (e.g. SECR) may include submodels
that constrain variation in core parameters and include the
effects of covariates. The language of generalised linear models is
convenient for describing submodels (e.g. Huggins 1989,
Lebreton et al. 1992). Each parameter is treated as a linear combination
of effects on its transformed ('link') scale. This is useful for
combining effects because, given a suitable link function, any
combination maps to a feasible value of the parameter. The logit scale
has this property for probabilities in (0,1), and the natural log scale
works for positive parameters i.e. (0, +Inf).

Submodels for spatially explicit capture--recapture in \pkg{secr} are defined
symbolically using the \R formula notation. A separate linear predictor
is used for each core parameter. Core parameters are 'real' parameters
in the terminology of MARK, and \pkg{secr} uses that term to reduce
confusion. Four real parameters are commonly modelled in \pkg{secr}: D
(density), g0, sigma and z. Only the last three real parameters, the
ones jointly defining detection probability as a function of location,
can be estimated directly when the model is fitted by maximizing the
conditional likelihood. D is then a derived parameter. 'z' is a shape
parameter used only for a 'hazard-rate' detection function (Hayes and
Buckland 1983). Other real parameters are used for acoustic models
(beta0, beta1; \url{../doc/secr-sound.pdf}) and for the mixture
proportion (pmix) in finite mixture models (\url{../doc/secr-finitemixtures.pdf}).

Each real parameter has a linear predictor of the form 

y = X * beta,

where y is vector of parameter values on the link scale, X is a design
matrix of predictor values, beta is a vector of coefficients, and '*'
stands for matrix multiplication. The elements of beta are estimated
when we fit the model; in MARK these are called 'beta parameters' to
distinguish them from the 'real' parameter values in y. X has one column
for each element of beta. To repeat: there is an X and a beta for each
real parameter; elsewhere in the documentation we use 'beta' to refer to
the vector got by concatenating \emph{all} the parameter-specific beta's. We now
describe design matrices in more detail.

[Some variations on the basic SECR model do not fit easily into this
framework. An example is the choice of detection function (halfnormal vs
hazard-rate). These are treated as higher-level choices.]

\bold{Design matrices}

The design matrix contains a column of '1's (for the constant or
intercept term) and additional columns as needed to describe the effects
in the submodel. Depending on the model, these may be continuous
predictors (e.g. air temperature to predict occasion-to-occasion
variation in g0), indicator variables (e.g. 1 if animal i was caught
before occasion s, 0 otherwise), or coded factor levels.

Within \code{secr.fit}, a design matrix is constructed automatically
from the input data (\code{capthist}) and the model formula (e.g.
\code{model$g0}) in a 2-stage process. First, a data frame is built
containing 'design data' with one column for each variable in the
formula. Second, the R function \code{model.matrix()} is used to
construct the design matrix. This process is hidden from the user. The
design matrix will have at least one more column than the design data,
and more if the formula includes interactions or factors with more than
two levels. For a good description of the general approach see the
documentation for RMark (Laake and Rexstad 2008). The key point is that
the necessary design data can be either extracted from the inputs
(\code{capthist} and \code{mask}) or generated automatically (e.g.
indicator of previous capture, mentioned in the previous paragraph).

Real parameters fall into two groups: density (D) and detection (g0,
sigma and z). Density and detection parameters are subject to different
types of effect, so they use different design matrices and are described
separately here:

\code{\link{secr detection models}}, \code{\link{secr density models}}

}

\references{ 

Laake, J. and Rexstad E. (2008) Appendix C. RMark - an alternative approach to building linear models in MARK. In: Cooch, E. and White, G. (eds) \emph{Program MARK: A Gentle Introduction.} 6th edition. Available online at \url{www.phidot.org}.

}

\note{

The structure of \pkg{secr} precludes certain types of model. Unlike
density, detection parameters (g0, sigma etc.) cannot be modelled as
varying in space \emph{per se}, whether continuously or discretely (e.g.
as a function of habitat class). However, such variation may be modelled
between detectors or between sessions. As an example, consider a measure
of vegetation cover in a 50-m circle centred on each detector. This may
be used as a detector covariate in models for g0 or sigma. A
'detector-centred' view of habitat effects is almost as sensible as an
'animal-centred' view; the one reservation is that the spatial scale
(radius of the circle) is arbitrary rather than being driven by sigma as
you might like. Perhaps this could be fixed in future versions by
computing the trap covariate 'on the fly' from covariates in the habitat
mask, given the current magnitude of sigma.

}

\keyword{ models }

