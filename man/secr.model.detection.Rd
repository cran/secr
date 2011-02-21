\name{secr.model.detection}
\alias{secr.model.detection}
\alias{secr detection models}

\title{ Models for Detection Parameters }

\description{

For spatially explicit capture--recapture estimation of a closed
population, we model the detection of individual \eqn{i} on occasion
\eqn{s} at detector \eqn{k}. Given \eqn{n} observed individuals on
\eqn{S} occasions at \eqn{K} detectors there are therefore \eqn{n.S.K}
detection probabilities of interest. We can think of these as elements
of a 3-dimensional array. Strictly, we are also interested in the
detection probabilities of unobserved individuals, but these are
estimated only by extrapolation from those observed so we do not
consider them in the array.

In a null (constant) model, all \eqn{n.S.K} detection probabilities are the
same. The conventional sources of variation in capture probability (Otis
et al. 1978) appear as variation in the \eqn{n} dimension (`individual
heterogeneity' h), in the \eqn{S} dimension (`time variation' t) or as a
particular interaction in these two dimensions (`behavioural response to
capture' b). Combined effects are possible.

Spatially explicit capture--recapture introduces two sorts of additional
complexity. Firstly, detection probability is no longer a scalar (even
for a particular animal, occasion and detector combination); it is
described by the detection function, which may have two parameters (e.g.
g0, sigma for half-normal), three parameters (e.g. g0, sigma, z for the
hazard-rate function), or potentially more.

Secondly, many more types of variation are possible. Any of the
parameters of the detection function may vary with respect to individual
(\eqn{i}), occasion (\eqn{s}) or detector (\eqn{k}). For example, there
may be a covariate associated with trap location that influences
detection probability.

The full design matrix for each detection submodel has one row for each
combination of \eqn{i}, \eqn{s} and \eqn{k} (animal, occasion and trap).
Allowing a distinct probability for each animal (the `\eqn{n}'
dimension) may seem excessive, as continuous individual-specific
covariates are feasible only when a model is fitted by maximizing the
conditional likelihood (cf Huggins 1989). However, the full \eqn{n.S.K}
array is convenient for coding both group membership (Lebreton et al.
1992, Cooch and White 2008) and experience of capture, even when
individual-level heterogeneity cannot be modelled.

Variation between `sessions' and between latent classes in a finite
mixture adds two further dimensions: in principle there is an
\eqn{n.S.K} array for each latent class (classes are numbered
1..\eqn{M}), and an \eqn{n.S.K.M} array for each session (sessions are
numbered 1..\eqn{R}). The full design matrix has \eqn{n.S.K.M.R}
rows. We do not expand on this here.

}

\section{ Specifying effects on detection parameters}{

Effects on parameters of detection probability are specified with \R
formulae using standard variable names or named covariates supplied by
the user. The formula for each detection parameter (g0, sigma, z) may be
constant (\eqn{\sim{~}}{~}1, the default) or some combination of terms
in standard \R formula notation (see \code{\link{formula}}).

\tabular{llll}{
Variable \tab Description \tab Data source \tab Dim \cr
t        \tab time factor (one level for each occasion) \tab automatic \tab \eqn{S} \cr
T        \tab time trend (integer covariate 0:(S-1)) \tab automatic \tab \eqn{S} \cr
tcov     \tab default time covariate \tab timecov[,1] \tab \eqn{S} \cr
kcov     \tab default trap covariate \tab covariates (traps)[,1] \tab \eqn{K} \cr
b        \tab learned response \tab \code{capthist} \tab \eqn{n.S} \cr
B        \tab transient (Markovian) response \tab \code{capthist} \tab \eqn{n.S} \cr
g        \tab group \tab see below \tab \eqn{n} \cr
h2       \tab 2-class mixture \tab -- \tab 2 \cr
h3       \tab 3-class mixture \tab -- \tab 2 \cr
session  \tab session factor (one level for each session) \tab automatic \tab \eqn{R} \cr
Session  \tab session number 0:(R-1) \tab automatic \tab \eqn{R} \cr
[user]   \tab individual covariate \tab \code{covariates (capthist)} \tab \eqn{n} \cr
[user]   \tab session covariate \tab \code{sessioncov} \tab \eqn{R} \cr
[user]   \tab time covariate \tab \code{timecov} \tab \eqn{S} \cr
[user]   \tab detector covariate \tab \code{covariates (traps)} \tab \eqn{K} \cr
}

The classic `learned response' is a step change following first
detection; this is implemented with the predictor variable `b' which is
FALSE up to and including the time of first capture and TRUE afterwards.
An alternative is a response that depends only on detection at the last
opportunity (`B').

Groups (`g') are defined by the interaction of the \code{capthist}
categorical (factor) individual covariates identified in \code{secr.fit}
argument `groups'. Groups are redundant with conditional likelihood
because individual covariates of whatever sort (continuous or
categorical) may be included freely in the model.

Individual heterogeneity (`h' in the notation of Otis et al. 1978) may
modelled by treating any detection parameter as a 2-part or 3-part finite mixture e.g. g0 \eqn{\sim{~}}{~}
  h2. See \url{../doc/secr-finitemixtures.pdf}.

Any other variable name appearing in a formula is assumed to refer to a
user-defined predictor. These will be interpreted by searching for name
matches in the dataframes of individual, session, time and trap
covariates, in that order (remembering that individual covariates other
than groups are allowed only when the model is fitted by maximizing the
conditional likelihood). The type of the predictor is inferred from the
data frame in which it first occurs. Thus if the model included the
formula `g0 ~ wetness', and `wetness' was a column in the data frame of
time covariates (timecov), then 'wetness' would be interpreted as a time
covariate, and a column of the same name in covariates(traps) would be
ignored. In this case, renaming the column in timecov would expose the
traps covariate, and `wetness' would be interpreted as an attribute of
detectors, rather than sample intervals. This is a good reason to give
covariates distinctive names!

The design matrix for detection parameters may also
be provided manually in the argument \code{dframe}. This feature is
untested.

The submodels for `g0', `sigma' and `z' are named components of the
\code{model} argument of \code{secr.fit}. They are expressed in \R
formula notation by appending terms to \eqn{\sim{~}}{~}.
The name of the response may optionally appear on the left hand side of
the formula (e.g. g0\eqn{\sim{~}}{~}b).

}

\examples{

## constant (null) model
list(g0 = ~1, sigma = ~1)

## both detection parameters change after first capture
list(g0 = ~b, sigma = ~b)

## group-specific parameters; additive time effect on g0
## groups are defined via the 'groups' argument of secr.fit
list(g0 = ~ g + t, sigma = ~ g)

## g0 depends on trap-specific covariate
list(g0 = ~ kcov)

}

\note{

The parameter `z' was previously called `b'; it was renamed to
avoid confusion with the predictor b used in a formula for a learned
trap response.

}

\references{

Cooch, E. and White, G. (eds) (2008) \emph{Program MARK: A Gentle Introduction.} 6th edition. Available online at \url{www.phidot.org}.

Hayes, R. J. and Buckland, S. T. (1983) Radial-distance models for the line-transect method. \emph{Biometrics} \bold{39}, 29--42.

Huggins, R. M. (1989) On the statistical analysis of capture experiments. \emph{Biometrika} \bold{76}, 133--140.

Lebreton, J.-D., Burnham, K. P., Clobert, J., Anderson, D. R. (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} \bold{62}, 67--118.

}

\seealso{ \code{\link{secr models}}, \code{\link{secr density models}}, \code{\link{secr.fit}} }

\keyword{ models }

