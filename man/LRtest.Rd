\name{LR.test}
\alias{LR.test}
\title{ Likelihood Ratio Test for SECR Models }

\description{
  Compute likelihood ratio test to compare two fitted models, one nested within the other.
}
\usage{
LR.test(secr1, secr2)
}
\arguments{
  \item{secr1}{ fitted secr model }
  \item{secr2}{ fitted secr model }
}
\details{
  The test statistic is twice the difference of the maximized likelihoods. It is compared to a chi-square distribution with df equal to the number of extra parameters in the more complex model.

  The models must be nested (no check is performed - this is up to the user), but either secr1 or secr2 may be the more general model.
}

\value{
  Object of class 'htest', a list with components
  \item{statistic }{value the test statistic }
  \item{parameter }{degrees of freedom of the approximate chi-squared distribution of the test statistic }
  \item{p.value }{probability of test statistic assuming chi-square distribution}
  \item{method }{character string indicating the type of test performed }
  \item{data.name }{character string with names of secr models compared }
}

\seealso{ \code{\link{AIC.secr}}, \code{\link{score.test}}  }
\examples{

## two pre-fitted models
AIC (secrdemo.0, secrdemo.b)
LR.test  (secrdemo.0, secrdemo.b)

}
\keyword{ htest }
