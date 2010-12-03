\name{secr.model.density}
\alias{secr.model.density}
\alias{secr density models}

\title{ Density Models}

\description{

SECR can fit an inhomogeneous Poisson model to describe the distribution
of animals.  This may be viewed as a surface of expected density across
the study area.

The log likelihood is evaluated in \code{secr.fit} by summing values at
points on a 'habitat mask'.  Each point in a habitat mask represents a
grid cell of potentially occupied habitat (their combined area may be
almost any shape and may include disjunct patches).  The full design
matrix for density (D) has one row for each point in the mask.  The
design matrix has one column for the intercept (constant) term and one
for each predictor.

Predictors may be based on Cartesian coordinates (e.g. 'x' for an
east-west trend), a continuous habitat variable (e.g. vegetation cover)
or a categorical (factor) habitat variable.  Predictors must be known for
all points in the mask (non-habitat excluded).  The variables 'x', 'y',
'session' and 'g' are provided automatically.  Other covariates should be
named columns in the 'covariates' attribute of the habitat mask.

The fitted model for density is linear on the link scale (see the
\code{link} argument of \code{\link{secr.fit}}. The default link for
density is 'log'.

\tabular{lll}{
Variable \tab Description \tab Data source \cr
x        \tab x-coordinate   \tab automatic \cr
y        \tab y-coordinate   \tab automatic \cr
session  \tab session factor \tab automatic \cr
Session  \tab session number 0:(\eqn{R}-1) \tab automatic \cr
g        \tab group factor   \tab automatic \cr
[user]   \tab mask covariate \tab \code{covariates (mask)} as named in formula \cr
}

The submodel for density (D) is a named component of the list used in
the \code{model} argument of \code{secr.fit}.  It is expressed in \R
formula notation by appending terms to \eqn{\sim{~}}{\~}.

}

\note{

This implementation is still experimental.  Note that no density model is
fitted when \code{secr.fit} is called with CL = TRUE.

}

\references{

Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics} \bold{64}, 377--385.

}

\seealso{ \code{\link{secr models}}, \code{\link{secr detection models}}, \code{\link{secr.fit}} }

\examples{
list(D = ~ 1)     ## constant density (homogeneous Poisson)
list(D = ~ x)     ## east-west trend
list(D = ~ cover) ## requires 'cover' as a mask covariate
}

\keyword{ models }

