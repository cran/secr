\name{derived}
\alias{derived}
\alias{esa}

\title{ Derived Parameters of Fitted SECR Model }

\description{
Compute derived parameters of spatially explicit capture-recapture model. 
Density is a derived parameter when a model is fitted by maximizing the conditional likelihood. So also is the effective sampling area (in the sense of Borchers and Efford 2008).
}
\usage{
derived(object, sessnum = NULL, groups = NULL, alpha = 0.05, 
    se.esa = FALSE, se.D = TRUE, loginterval = TRUE, 
    distribution = NULL)
esa(object, sessnum = 1, beta = NULL, real = NULL)
}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit} }
  \item{sessnum}{ index of session in object\$capthist for which output required }
  \item{groups}{ indices defining group(s) (see Details) }
  \item{alpha}{ alpha level for confidence intervals }
  \item{se.esa}{ logical for whether to calculate SE(mean(esa)) }
  \item{se.D}{ logical for whether to calculate SE(D-hat) }
  \item{loginterval}{ logical for whether to base interval on log(D)}
  \item{distribution}{ character string for distribution of the number of individuals detected }
  \item{beta}{ vector of fitted parameters on transformed (link) scale }
  \item{real}{ vector of 'real' parameters } 
}

\details{

The derived estimate of density is a Horvitz-Thompson-like estimate:
\deqn{\hat{D} = \sum\limits _{i=1}^{n} {a_i (\hat{\theta})^{-1}}}{ D-hat
= sum( a_i (theta-hat)^--1)} where \eqn{a_i (\hat{\theta})}{% a_i
(theta-hat)} is the estimate of effective sampling area for animal
\eqn{i} with detection parameter vector \eqn{\theta}.

A non-null value of the argument \code{distribution} overrides the value
in \code{object$details}. The sampling variance of \eqn{\hat{D}}{D-hat}
from \code{secr.fit} by default is spatially unconditional (\code{distribution}
= 'Poisson'). For sampling variance conditional on the population of the
habitat mask (and therefore dependent on the mask area), specify
\code{distribution} = 'binomial'. The equation for the conditional
variance includes a factor \eqn{(1 - a/A)} that disappears in the
unconditional (Poisson) variance (Borchers and Efford 2007). Thus the
conditional variance is always less than the unconditional variance. The
unconditional variance may in turn be an overestimate or (more likely)
an underestimate if the true spatial variance is non-Poisson.

Derived parameters may be estimated for population subclasses (groups)
defined by the user with the \code{groups} argument. Each named factor
in \code{groups} should appear in the covariates dataframe of
object\$capthist (or each of its components, in the case of a
multi-session dataset).

The effective sampling area 'esa' reported by \code{derived} is equal
to the mean of the \eqn{a_i (\hat{\theta})}{a_i (theta-hat)}.

A 100(1--alpha)\% asymptotic confidence interval is reported for
density. By default, this is asymmetric about the estimate because the
variance is computed by backtransforming from the log scale. You may 
also choose a symmetric interval (variance obtained on natural scale). 

\code{esa} is used by \code{derived} to compute individual-specific
effective sampling areas: \deqn{a_i (\hat{\theta}) = \int _A \:
p.(\mathbf{X};\mathbf{z}_i, \mathbf{\hat{\theta}}) \; \mathrm{d}
\mathbf{X}}{a_i = integral p.(X; z_i, theta_i) dX} where
\eqn{p.(\mathbf{X})}{p.(X)} is the probability an individual at X is
detected at least once and the \eqn{\mathbf{z}_i}{z_i} are optional
individual covariates. Integration is over the area \eqn{A} of the
habitat mask.

The vector of detection parameters for \code{esa} may be specified via
\code{beta} or \code{real}, with the former taking precedence. If
neither is provided then the fitted values in \code{object$fit$par} are
used. Specifying \code{real} parameter values bypasses the various
linear predictors. Strictly, the 'real' parameters are for a naive
capture (animal not detected previously).

The computation of sampling variances is relatively slow and may be
suppressed with \code{se.esa} and \code{se.D} as desired.

}

\value{

Dataframe with one row for each derived parameter ('esa', 'D') and
columns as below

\tabular{ll}{
estimate    \tab estimate of derived parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
varcomp1    \tab variance due to variation in \eqn{n} (Huggins' \eqn{s^2}) \cr
varcomp2    \tab variance due to uncertainty in estimates of detection parameters \cr
}

For a multi-session or multi-group analysis the value is a list
with one component for each session and group.

}

\references{
Borchers, D. L. and Efford, M. G. (2007) Supplements to Biometrics paper. Available online at \url{http://www.otago.ac.nz/density}.

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics}, \bold{64}, 377--385.

Huggins, R. M. (1989) On the statistical analysis of capture experiments. \emph{Biometrika} \bold{76}, 133--140.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{predict.secr}}, \code{\link{print.secr}}, \code{\link{secr.fit}} }

\examples{
## extract derived parameters from a model fitted previously
## by maximizing the conditional likelihood 
data(secrdemo)
derived (secrdemo.CL)

## what happens when sampling variance is conditional on mask N?
derived(secrdemo.CL, distribution = 'binomial')

## fitted g0, sigma
esa(secrdemo.CL)
## force different g0, sigma
esa(secrdemo.CL, real = c(0.2, 25))

}

\keyword{ models }

