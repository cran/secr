\name{print.secr}
\alias{print.secr}

\title{ Print secr Object }

\description{
  Print results from fitting a spatially explicit capture--recapture model.
}

\usage{
\method{print}{secr} (x, newdata = NULL, alpha = 0.05, deriv = FALSE,
call = TRUE, ...)
}

\arguments{
  \item{x}{ \code{secr} object output from \code{secr.fit} }
  \item{newdata}{ optional dataframe of values at which to evaluate model }
  \item{alpha}{ alpha level }
  \item{deriv}{ logical for calculation of derived D and esa }
  \item{call}{ logical; if TRUE the call is printed }
  \item{...}{ other arguments (not used currently) }
}

\details{

Results are potentially complex and depend upon the analysis (see
below). Optional \code{newdata} should be a dataframe with a column for
each of the variables in the model. If \code{newdata} is missing then a
dataframe is constructed automatically.  Default \code{newdata} are for
a naive animal on the first occasion; numeric covariates are set to zero
and factor covariates to their base (first) level. Confidence intervals
are 100 (1 -- alpha) \% intervals.

\tabular{ll}{
call            \tab the function call (optional) \cr
version,time    \tab secr version, and date and time fitting started \cr

Detector type   \tab `single', `multi', `proximity' etc. \cr
Detector number \tab number of detectors \cr
Average spacing \tab \cr
x-range         \tab \cr
y-range         \tab \cr
New detector type \tab as fitted when details$newdetector specified \cr
N animals       \tab number of distinct animals detected \cr
N detections    \tab number of detections \cr
N occasions     \tab number of sampling occasions \cr
Mask area       \tab \cr

Model          \tab model formula for each `real' parameter \cr
Fixed (real)   \tab fixed real parameters \cr
Detection fn   \tab detection function type (halfnormal or hazard-rate) \cr
N parameters   \tab number of parameters estimated  \cr
Log likelihood \tab log likelihood \cr
AIC            \tab Akaike's information criterion \cr
AICc           \tab AIC with small sample adjustment (Burnham and Anderson 2002) \cr
Beta parameters \tab coef of the fitted model, SE and confidence intervals \cr
vcov            \tab variance-covariance matrix of beta parameters  \cr
Real parameters \tab fitted (real) parameters evaluated at base levels of covariates \cr
Derived parameters \tab derived estimates of density and mean effective sampling area (optional) \cr
}

Derived parameters (see \code{\link{derived}}) are computed only if
\code{deriv = TRUE}.

}
\references{ 
Burnham, K. P. and Anderson, D. R. (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. Second edition. New York: Springer-Verlag.
}

\seealso{ \code{\link{AIC.secr}}, \code{\link{secr.fit}} }

\examples{

## load & print previously fitted null (constant parameter) model
print(secrdemo.0)

\dontrun{
print(secrdemo.CL, deriv = TRUE)
}
}

\keyword{ print }
