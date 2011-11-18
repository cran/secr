\name{model.average}
\alias{model.average}
\alias{collate}

\title{ Model averaging for SECR Models }

\description{
AICc-weighted average of estimated `real' or `beta' parameters from multiple fitted secr models.
}

\usage{

model.average(..., realnames = NULL, betanames = NULL, newdata = NULL, 
    alpha = 0.05, dmax = 10, covar = FALSE, average = "link")

collate (..., realnames = NULL, betanames = NULL, newdata = NULL, 
    scaled = FALSE, alpha = 0.05, perm = 1:4, fields = 1:4)

}

\arguments{
  \item{\dots}{ \code{secr} objects }
  \item{realnames}{ character vector of real parameter names }
  \item{betanames}{ character vector of beta parameter names }
  \item{newdata}{ optional dataframe of values at which to evaluate models }
  \item{scaled}{ logical for scaling of sigma and g0 (see Details) }
  \item{alpha}{ alpha level for confidence intervals }
  \item{dmax}{ numeric, the maximum AIC difference for inclusion in confidence set }
  \item{covar}{ logical, if TRUE then return variance-covariance matrix }
  \item{average}{ character string for scale on which to average real parameters }
  \item{perm}{ permutation of dimensions in output from \code{collate} }
  \item{fields}{ vector to restrict summary fields in output }
}

\details{

Models to be compared must have been fitted to the same data and use the
same likelihood method (full vs conditional). If \code{realnames} ==
NULL and \code{betanames} == NULL then all real parameters will be
averaged; in this case all models must use the same real parameters. To
average beta parameters, specify \code{betanames} (this is ignored if a
value is provided for \code{realnames}). See \code{\link{predict.secr}}
for an explanation of the optional argument \code{newdata};
\code{newdata} is ignored when averaging beta parameters. 

Model-averaged estimates for parameter \eqn{\theta}{theta} are given by
\deqn{\hat{\theta} = \sum\limits _k w_k \hat{\theta}_k}{theta-hat = sum(
w_k * theta-hat_k)} where the subscript \eqn{k} refers to a specific
model and the \eqn{w_k} are AIC weights with small sample adjustment
(see \code{\link{AIC.secr}} for details). Averaging of real parameters
may be done on the link scale before back-transformation
(\code{average="link"}) or after back-transformation
(\code{average="real"}). 

Models for which dAICc > \code{dmax} are given a weight of zero and
effectively are excluded from averaging.

Also, \deqn{\mbox{var} (\hat{\theta}) = \sum\limits _{k} { w_{k}
( \mbox{var}(\hat{\theta}_{k} | \beta _k) + \beta _k ^2)}
}{var(theta-hat) = sum(w_k (var(theta-hat_k) + beta_k^2))}

where \eqn{\hat{\beta} _k = \hat{\theta}_k - \hat{\theta}}{beta-hat_k =
theta-hat_k -- theta-hat} and the variances are asymptotic estimates
from fitting each model \eqn{k}. This follows Burnham and Anderson
(2004) rather than Buckland et al. (1997).

\code{collate} extracts parameter estimates from a set of fitted secr
model objects. \code{fields} may be used to select a subset of summary
fields ("estimate","SE.estimate","lcl","ucl") by name or number. 

The argument \code{scaled} applies only to the detection parameters g0
and sigma, and only to models fitted with \code{scalesigma} or
\code{scaleg0} switched on. If \code{scaled} is TRUE then each estimate is multiplied by
its scale factor (1/D^0.5 and 1/sigma^2 respectively).

}

\value{

For \code{model.average}, an array of model-averaged estimates, their
standard errors, and a \eqn{100(1-\alpha)}{100(1-alpha)}\% confidence
interval. The interval for real parameters is backtransformed from the
link scale. If there is only one row in \code{newdata} or beta
parameters are averaged or averaging is requested for only one parameter
then the array is collapsed to a matrix. If \code{covar = TRUE} then a
list is returned with separate components for the estimates and the
variance-covariance matrices.

For \code{collate}, a 4-dimensional array of model-specific parameter
estimates. By default, the dimensions correspond respectively to rows in
\code{newdata} (usually sessions), models, statistic fields (estimate, SE.estimate, lcl,
ucl), and parameters ("D", "g0" etc.). For particular comparisons it often helps
to reorder the dimensions with the \code{perm} argument.

}
\references{ 

Buckland S. T., Burnham K. P. and Augustin, N. H. (1997) Model
selection: an integral part of inference. \emph{Biometrics} \bold{53},
603--618.

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model Selection and
Multimodel Inference: A Practical Information-Theoretic Approach}.
Second edition. New York: Springer-Verlag.

Burnham, K. P. and Anderson, D. R. (2004) Multimodel inference -
understanding AIC and BIC in model selection. \emph{Sociological Methods
& Research} \bold{33}, 261--304.

}

\seealso{ \code{\link{AIC.secr}}, \code{\link{secr.fit}}}

\examples{
## Compare two models fitted previously
## secrdemo.0 is a null model
## secrdemo.b has a learned trap response

model.average(secrdemo.0, secrdemo.b)
model.average(secrdemo.0, secrdemo.b, betanames = c("D","g0","sigma"))

## In this case we find the difference was actually trivial...
## (subscripting of output is equivalent to setting fields = 1)

collate (secrdemo.0, secrdemo.b, perm = c(4,2,3,1))[,,1,]

}

\keyword{ models }

