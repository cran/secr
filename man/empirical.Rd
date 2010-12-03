\name{empirical.varD}
\alias{empirical.varD}
\title{
Design-based Variance of Density
}
\description{

  Compute the variance of a density estimate from a fitted spatially
  explicit model, using the empirical spatial variance of the number
  observed in replicated sampling units.

}

\usage{
empirical.varD(object, esa = NULL, se.esa = NULL)
}

\arguments{
  \item{object}{
fitted secr model object or vector of number observed in each sampling unit
}
  \item{esa}{
scalar estimate of effective sampling area (\eqn{\hat{a}}{a-hat})
}
  \item{se.esa}{
estimated standard error of effective sampling area (\eqn{\widehat{SE}(\hat{a})}{SE-hat(a-hat)})
}
}

\details{

  The variance of a Horvitz-Thompson-like estimate of density may be
  estimated as the sum of two components, one due to uncertainty in the
  estimate of effective sampling area (\eqn{\hat{a}}{a-hat}) and the
  other due to spatial variance in the total number of animals \eqn{n}
  observed on \eqn{J} replicate sampling units (\eqn{n =
    \sum_{j=1}^{J}{n_j}}{sum(n_j)}). We use a delta-method approximation
  that assumes independence of the components:

  \deqn{ \widehat{\mbox{var}}(\hat{D}) = \hat{D}^2
    \{\frac{\widehat{\mbox{var}}(n)}{n^2} +
    \frac{\widehat{\mbox{var}}(\hat{a})}{\hat{a}}\} }{var(D-hat) =
    D-hat^2 (var(n)/n^2 + var(a-hat)/a-hat^2) }

where \eqn{\widehat{\mbox{var}}(n) = \frac{J}{J-1}
  \sum_{j=1}^{J}(n_j-n/J)^2}{var(n) = J/(J-1).sum((n_j - n/J)^2)}. The
  estimate of \eqn{\mbox{var}(\hat{a})}{var(a-hat)} is model-based while
  that of \eqn{\mbox{var}(n)}{var(n)} is design-based. This formulation follows
  that of Buckland et al. (2001, p. 78) for conventional distance
  sampling. Given sufficient independent replicates, it is a robust way
  to allow for unmodelled spatial overdispersion.

  There is a complication in SECR owing to the fact that
  \eqn{\hat{a}}{a-hat} is a derived quantity (actually an integral)
  rather than a model parameter. Its sampling variance
  \eqn{\mbox{var}(\hat{a})}{var(a-hat)} is estimated indirectly in
  \pkg{secr} by combining the asymptotic estimate of the covariance
  matrix of the fitted detection parameters \eqn{\theta}{theta} with a
  numerical estimate of the gradient of \eqn{a(\theta)}{a(theta)} with
  respect to \eqn{\theta}{theta}. This calculation is performed in
  \code{\link{derived}}.

  Input may be either a single fitted model or a vector of counts (the
  \eqn{n_j}) along with \eqn{\hat{a}}{a-hat} and
  \eqn{\widehat{SE}(\hat{a})}{SE(a-hat)}. If \code{object} is a single
  fitted model then it must span multiple sessions and each session is
  interpreted as a replicate sample. This form does not allow zero
  counts; if some units have zero counts use the second style of input.
  }

\value{
A vector of values --
\item{D }{Horvitz-Thompson-like estimate of population density}
\item{seD }{SE of density estimate}
\item{CVD }{relative SE of density estimate}
\item{CVn }{relative SE of number observed (across sampling units) }
\item{CVa }{relative SE of effective sampling area }
}
\references{

  Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers,
D. L. and Thomas, L. (2001) \emph{Introduction to Distance Sampling:
Estimating Abundance of Biological Populations}. Oxford University
Press, Oxford.

}

\seealso{
\code{\link{derived}}, \code{\link{esa}} 
}

\examples{

## The `ovensong' data are pooled from 75 replicate positions of a
## 4-microphone array. The array positions are coded as the first 4
## digits of each sound identifier. The sound data are initially in the
## object 'signalCH'. We first impose a 52.5 dB signal threshold as in
## Dawson & Efford (2009, J. Appl. Ecol. 46:1201--1209). The vector nj
## includes 33 positions at which no ovenbird was heard. The first and
## second columns of `temp' hold the estimated effective sampling area
## and its standard error.

signalCH.525 <- subset(signalCH, cutval = 52.5)
nonzero.counts <- table(substring(rownames(signalCH.525),1,4))
nj <- c(nonzero.counts, rep(0, 75 - length(nonzero.counts)))
temp <- derived(ovensong.model.1, se.esa = TRUE)
empirical.varD(nj, temp['esa',1], temp['esa',2])

## The result is very close to that reported by Dawson & Efford
## from a 2-D Poisson model fitted by maximizing the full likelihood.

}

\keyword{ models }