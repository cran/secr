\name{predict.secr}
\alias{predict.secr}
\alias{predict.secrlist}
\alias{detectpar}

\title{ SECR Model Predictions }

\description{
Evaluate a spatially explicit capture--recapture model. That is, compute the `real' parameters corresponding to the `beta' parameters of a fitted model for arbitrary levels of any variables in the linear predictor.
}

\usage{
\method{predict}{secr}(object, newdata = NULL, type = c("response", "link"), se.fit = TRUE,
    alpha = 0.05, savenew = FALSE, scaled = FALSE, \dots)

\method{predict}{secrlist}(object, newdata = NULL, type = c("response", "link"), se.fit = TRUE,
    alpha = 0.05, savenew = FALSE, scaled = FALSE, \dots)

detectpar (object, \dots) 

}

\arguments{
  \item{object}{ \code{secr} object output from \code{secr.fit}, or list
  of \code{secr} objects (secrlist) }
  \item{newdata}{ optional dataframe of values at which to evaluate
    model }
  \item{type}{ character; type of prediction required. The default ("response") provides estimates of the 'real' parameters.} 
  \item{se.fit}{ logical for whether output should include SE and confidence intervals }
  \item{alpha}{ alpha level for confidence intervals }
  \item{savenew}{ logical for whether newdata should be saved }
  \item{scaled}{ logical for scaling of sigma and g0 (see Details) }
  \item{\dots}{ other arguments }
}
\details{
  
The variables in the various linear predictors are described in
\link{secr models} and listed for the particular model in the
\code{vars} component of \code{object}.

Optional \code{newdata} should be a dataframe with a column for each of
the variables in the model (see `vars' component of \code{object}).  If
\code{newdata} is missing then a dataframe is constructed automatically.
Default \code{newdata} are for a naive animal on the first occasion;
numeric covariates are set to zero and factor covariates to their base
(first) level.

Standard errors for parameters on the response (real) scale are by the
delta method (Lebreton et al. 1992), and confidence intervals are
backtransformed from the link scale.

The argument \code{scaled} applies only to the detection parameters g0
and sigma, and only to models fitted with \code{scalesigma} or
\code{scaleg0} switched on. If \code{scaled} is TRUE then each estimate is multiplied by
its scale factor (1/D^0.5 and 1/sigma^2 respectively).

The value of \code{newdata} is optionally saved as an attribute.

\code{detectpar} is used to extract the detection parameter estimates
from a simple model to pass to functions such as
\code{esa.plot}. \code{detectpar} calls \code{predict.secr}. Parameters
will be evaluated by default at base levels of the covariates, although
this may be overcome by passing a one-line \code{newdata} to
\code{predict} via the \dots argument. Groups and mixtures are a
headache for \code{detectpar}: it merely returns the estimated detection
parameters of the first group or mixture.

}

\value{
  When \code{se.fit} = FALSE, a dataframe identical to \code{newdata} except for the addition of one column for each `real' parameter. Otherwise, a list with one component for each row in \code{newdata}. Each component is a dataframe with one row for each `real' parameter (density, g0, sigma, b) and columns as below

\tabular{ll}{
link        \tab link function \cr 
estimate    \tab estimate of real parameter \cr
SE.estimate \tab standard error of the estimate \cr
lcl         \tab lower 100(1--alpha)\% confidence limit \cr
ucl         \tab upper 100(1--alpha)\% confidence limit \cr
}

When \code{newdata} has only one row, the structure of the list is
`dissolved' and the return value is one data frame.

For \code{detectpar}, a list with the estimated values of detection
parameters (e.g., g0 and sigma if detectfn = "halfnormal"). In the case
of multi-session data the result is a list of lists (one list per
session).

}

\note{

  \code{\link{predictDsurface}} should be used for predicting density at many
  points from a model with spatial variation. This deals automatically
  with scaling of x- and y-coordinates, and is much is faster than
  predict.secr. The resulting Dsurface object has its own plot method.
  
}

\references{ 
Lebreton, J.-D., Burnham, K. P., Clobert, J., Anderson, D. R. (1992) Modeling survival and testing biological hypotheses using marked animals: a unified approach with case studies. \emph{Ecological Monographs} \bold{62}, 67--118. }

\seealso{ \code{\link{secr.fit}}, \code{\link{predictDsurface}} }

\examples{

## load previously fitted secr model with trap response
## and extract estimates of `real' parameters for both
## naive (b = 0) and previously captured (b = 1) animals

predict (secrdemo.b, newdata = data.frame(b=0:1))

temp <- predict (secrdemo.b, newdata = data.frame(b=0:1), 
    save = TRUE)
attr(temp, "newdata")

detectpar(secrdemo.0)
}

\keyword{ models }

