\name{score.test}
\alias{score.test}
\alias{score.table}

\title{ Score Test for SECR Models }
\description{
  Compute score tests comparing a fitted model and a more general alternative model.
}
\usage{
score.test(secr, ..., betaindex = NULL, trace = FALSE)

score.table(object, ..., sort = TRUE, dmax = 10)

}

\arguments{
  \item{secr}{ fitted secr model }
  \item{...}{ one or more alternative models OR a fitted secr model }
  \item{trace}{ logical. If TRUE then output one-line summary at each evaluation of the likelihood }
  \item{betaindex}{ vector of indices mapping fitted values to parameters in the alternative model }

  \item{object}{ score.test object or list of such objects }
  \item{sort}{ logical for whether output rows should be in descending order of AICc }
  \item{dmax}{ threshold of dAICc for inclusion in model set }

}
\details{

Score tests allow fast model selection (e.g. Catchpole & Morgan 1996).
Only the simpler model need be fitted. This implementation uses the
observed information matrix, which may sometimes mislead (Morgan et al.
2007). The gradient and second derivative of the likelihood function are
evaluated numerically at the point in the parameter space of the second
model corresponding to the fit of the first model. This operation uses
the function \code{fdHess} of the \pkg{nlme} package; the likelihood
must be evaluated several times, but many fewer times than would be
needed to fit the model. The score statistic is an approximation to the
likelihood ratio; this allows the difference in AIC to be estimated.

Mapping of parameters between the fitted and alternative models
sometimes requires user intervention via the \code{betaindex} argument.
For example \code{betaindex} = c(1,2,4) is the correct mapping when
comparing the null model (D\eqn{\sim{~}}{~}1, g0\eqn{\sim{~}}{~}1,
sigma\eqn{\sim{~}}{~}1) to one with a behavioural effect on g0
((D\eqn{\sim{~}}{~}1, g0\eqn{\sim{~}}{~}b, sigma\eqn{\sim{~}}{~}1).

\code{score.table} summarises one or more score tests in the form of a
model comparison table. The \dots argument here allows the inclusion of
additional score test objects (note the meaning differs from
\code{score.test}). Approximate AICc values are used to compute relative
AIC model weights for all models within dmax AICc units of the best
model.

}
\value{
  An object of class 'score.test' that inherits from 'htest', a list with components
  \item{statistic }{the value the chi-squared test statistic (score statistic)}
  \item{parameter }{degrees of freedom of the approximate chi-squared distribution of the test statistic (difference in number of parameters H0, H1)}
  \item{p.value }{probability of test statistic assuming chi-square distribution}
  \item{method }{a character string indicating the type of test performed }
  \item{data.name }{character string with null hypothesis, alternative hypothesis and arguments to function call from fit of H0 }
  \item{H0 }{simpler model }
  \item{np0 }{number of parameters in simpler model }
  \item{H1 }{alternative model }
  \item{H1.beta }{coefficients of alternative model }
  \item{AIC }{Akaike's information criterion, approximated from score statistic }
  \item{AICc }{AIC with small-sample adjustment of Hurvich & Tsai 1989 }

  If \dots defines several alternative models then a list of score.test objects is returned.

  The output from \code{score.table} is a dataframe with one row per model, including the reference model. 

}

\note{

This implementation is experimental. The AIC values, and values derived
from them, are approximations that may differ considerably from AIC
values obtained by fitting and comparing the respective models. Use of
the observed information matrix may not be optimal.

\code{score.test} cannot be used to compare models that differ in the arguments \code{scalesigma} or \code{scaleg0}.

}

\references{ 
Catchpole, E. A. and Morgan, B. J. T. (1996) Model selection of ring-recovery models using score tests. \emph{Biometrics} \bold{52}, 664--672.

Hurvich, C. M. and Tsai, C. L. (1989) Regression and time series model selection in small samples. \emph{Biometrika} \bold{76}, 297--307.

Morgan, B. J. T., Palmer, K. J. and Ridout, M. S. (2007) Negative score test statistic. \emph{American statistician} \bold{61}, 285--288.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz}  }
\seealso{ \code{\link{AIC}}, \code{\link{LR.test}} }
\examples{

\dontrun{
    data(secrdemo)
    AIC (secrdemo.0, secrdemo.b)
    st <- score.test (secrdemo.0, g0 ~ b)
    st
    score.table(st)
}

}
\keyword{ models }
