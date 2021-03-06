\name{AIC.secr}
\alias{AIC.secr}
\alias{logLik.secr}
\alias{AIC.secrlist}
\alias{secrlist}
\alias{[.secrlist}

\title{ Compare SECR Models }

\description{

Terse report on the fit of one or more spatially explicit
capture--recapture models. Models with smaller values of AIC (Akaike's
Information Criterion) are preferred. Extraction ([) and logLik methods
are included.

}

\usage{
\method{AIC}{secr}(object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c("AICc","AIC"))
\method{AIC}{secrlist}(object, ..., sort = TRUE, k = 2, dmax = 10, criterion = c("AICc","AIC"))
\method{logLik}{secr}(object, ...)
secrlist(...)
\method{[}{secrlist}(x, i)
}

\arguments{
  \item{object}{ \code{secr} object output from the function
      \code{\link{secr.fit}}, or a list of such objects with class c("secrlist", "list") }
  \item{\dots}{ other \code{secr} objects }
  \item{sort}{ logical for whether rows should be sorted by ascending AICc }
  \item{k}{ numeric, penalty per parameter to be used; always k = 2 in this method}
  \item{dmax}{ numeric, maximum AIC difference for inclusion in confidence set}
  \item{criterion}{ character, criterion to use for model comparison and
  weights}
  \item{x}{secrlist}
  \item{i}{indices}
}

\details{

Models to be compared must have been fitted to the same data and use the
same likelihood method (full vs conditional). From version 4.1 a warning is 
issued if \code{\link{AICcompatible}} reveals a problem.

AIC with small sample adjustment is given by

\deqn{ \mbox{AIC}_c = -2\log(L(\hat{\theta})) + 2K +
\frac{2K(K+1)}{n-K-1} }{AICc = -2log(L(theta-hat)) + 2K +
2K(K+1)/(n-K-1)}

where \eqn{K} is the number of "beta" parameters estimated. The sample size
\eqn{n} is the number of individuals observed at least once (i.e. the
number of rows in \code{capthist}).

Model weights are calculated as \deqn{w_i = \frac{\exp(-\Delta_i / 2)}{
\sum{\exp(-\Delta_i / 2)}}}{w_i = exp(-delta_i / 2) / sum{ exp(-delta_i
/ 2) }}, where \eqn{\Delta}{delta} refers to differences in AIC or AICc depending on the
argument `criterion'.

Models for which delta > \code{dmax} are given a weight of zero and are
excluded from the summation. Model weights may be used to form
model-averaged estimates of real or beta parameters with
\code{\link{model.average}} (see also Buckland et al. 1997, Burnham and
Anderson 2002).

The argument \code{k} is included for consistency with the generic method \code{AIC}.

\code{secrlist} forms a list of fitted models (an object of class
`secrlist') from the fitted models in \dots.  Arguments may include
secrlists. If secr components are named the model names will be retained
(see Examples).

}

\value{
  A data frame with one row per model.  By default, rows are sorted by ascending AICc.
  \item{model }{character string describing the fitted model}
  \item{detectfn }{shape of detection function fitted (halfnormal vs hazard-rate) }
  \item{npar }{number of parameters estimated}
  \item{logLik }{maximized log likelihood}
  \item{AIC }{Akaike's Information Criterion}
  \item{AICc }{AIC with small-sample adjustment of Hurvich & Tsai
    (1989)}
  And depending on \code{criterion}:
  \item{dAICc }{difference between AICc of this model and the one with smallest AICc}
  \item{AICcwt }{AICc model weight}
  or
  \item{dAIC }{difference between AIC of this model and the one with smallest AIC}
  \item{AICwt }{AIC model weight}

  \code{logLik.secr} returns an object of class `logLik' that has
  attribute \code{df} (degrees of freedom = number of estimated
  parameters).
  
}
\references{ 

Buckland S. T., Burnham K. P. and Augustin, N. H. (1997) Model selection: an integral part of inference. \emph{Biometrics} \bold{53}, 603--618.

Burnham, K. P. and Anderson, D. R. (2002) \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach}. Second edition. New York: Springer-Verlag.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model selection in small samples. \emph{Biometrika} \bold{76}, 297--307.

Turek, D. and Fletcher, D. (2012) Model-averaged Wald confidence
intervals. \emph{Computational statistics and data analysis} \bold{56},
2809--2815.

}

\note{

It is not be meaningful to compare models by AIC if they relate to
different data (see \code{\link{AICcompatible}}).

Specifically:

\itemize{

  \item an `secrlist' generated and saved to file by \code{mask.check}
may be supplied as the object argument of \code{AIC.secrlist}, but the
results are not informative

\item models fitted by the conditional likelihood (\code{CL = TRUE}) and
full likelihood (\code{CL = FALSE}) methods cannot be compared

\item hybrid mixture models (using hcov argument of secr.fit) should not
be compared with other models

\item grouped models (using groups argument of secr.fit) should not be
compared with other models

\item multi-session models should not be compared with single-session
models based on the same data.

}

A likelihood-ratio test (\code{\link{LR.test}}) is a more direct way to
compare two models.

The issue of goodness-of-fit and possible adjustment of AIC for
overdispersion has yet to be addressed (cf QAIC in MARK).

From version 2.6.0 the user may select between AIC and AICc for
comparing models, whereas previously only AICc was used and AICc weights
were reported as `AICwt'). There is evidence that AIC may be better for
model averaging even when samples are small sizes - Turek and Fletcher
(2012).

}

\seealso{ 

\code{\link{AICcompatible}}, 
\code{\link{model.average}}, 
\code{\link[stats]{AIC}}, 
\code{\link{secr.fit}}, 
\code{\link{print.secr}}, 
\code{\link{score.test}}, 
\code{\link{LR.test}}, 
\code{\link{deviance.secr}} }

\examples{
## Compare two models fitted previously
## secrdemo.0 is a null model
## secrdemo.b has a learned trap response

AIC(secrdemo.0, secrdemo.b)

## Form secrlist and pass to AIC.secr
temp <- secrlist(null = secrdemo.0, learnedresponse = secrdemo.b)
AIC(temp)

}

\keyword{ models }

