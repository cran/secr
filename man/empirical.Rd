\name{empirical.varD}
\alias{empirical.varD}
\alias{derived.nj}
\alias{derived.session}
\alias{derived.mash}
\alias{derived.cluster}
\alias{derived.external}
\title{
Empirical Variance of H-T Density Estimate
}
\description{

  Compute Horvitz-Thompson-like estimate of population density from a
  previously fitted spatial detection model, and estimate its sampling
  variance using the empirical spatial variance of the number observed
  in replicate sampling units. Wrapper functions are provided for
  several different scenarios, but all ultimately call
  \code{derived.nj}. The function \code{derived} also computes
  Horvitz-Thompson-like estimates, but it assumes a Poisson or binomial
  distribution of total number when computing the sampling variance.

}

\usage{

derived.nj ( nj, esa, se.esa, method = c("SRS", "local", "poisson",
    "binomial"), xy = NULL, alpha = 0.05, loginterval = TRUE, area =
    NULL )

derived.mash ( object, sessnum = NULL,  method = c("SRS", "local"),
    alpha = 0.05, loginterval = TRUE)

derived.cluster ( object, sessnum = NULL,  method = c("SRS", "local"),
    alpha = 0.05, loginterval = TRUE)

derived.session ( object,  method = c("SRS", "local"), xy = NULL,
    alpha = 0.05, loginterval = TRUE )

derived.external ( object, sessnum = NULL, nj, cluster, buffer = 100,
    mask = NULL, noccasions = NULL,  method = c("SRS", "local"), xy = NULL,
    alpha = 0.05, loginterval = TRUE)

}

\arguments{
  \item{object}{fitted secr model }
  \item{nj}{ vector of number observed in each sampling unit (cluster)}
  \item{esa}{scalar estimate of effective sampling area (\eqn{\hat{a}}{a-hat})}
  \item{se.esa}{ estimated standard error of effective sampling
    area (\eqn{\widehat{SE}(\hat{a})}{SE-hat(a-hat)})}
  \item{method}{ character string `SRS' or `local' }
  \item{xy}{ dataframe of x- and y- coordinates (\code{method = "local"} only)}
  \item{alpha}{ alpha level for confidence intervals }
  \item{loginterval}{ logical for whether to base interval on log(N)}
  \item{area}{ area of region for method = "binomial" (hectares) }
  \item{sessnum}{ index of session in object$capthist for which output required }
  \item{cluster}{ `traps' object for a single cluster}
  \item{buffer}{ width of buffer in metres (ignored if \code{mask}
  provided) } 
  \item{mask}{ mask object for a single cluster of detectors }
  \item{noccasions}{number of occasions (for \code{nj}) }
  
}

\details{

  \code{derived.nj} accepts a vector of counts (\code{nj}), along with
  \eqn{\hat{a}}{a-hat} and \eqn{\widehat{SE}(\hat{a})}{SE(a-hat)}.  The
  argument \code{esa} may include both \eqn{\hat{a}}{a-hat} and
  \eqn{\widehat{SE}(\hat{a})}{SE(a-hat)}) - any form will do if it can
  be coerced to a vector of length 2. In the special case that \code{nj}
  is of length 1, or \code{method} takes the values `poisson' or
  `binomial', the variance is computed using a theoretical variance
  rather than an empirical estimate. The value of \code{method}
  corresponds to `distribution' in \code{\link{derived}}, and defaults to
  `poisson'. For \code{method = 'binomial'} you must specify \code{area}
  (see Examples).

  \code{derived.cluster} accepts a model fitted to data from clustered
  detectors; each \link{cluster} is interpreted as a replicate
  sample. It is assumed that the sets of individuals sampled by
  different clusters do not intersect, and that all clusters have the
  same geometry (spacing, detector number etc.).
  
  \code{derived.mash} accepts a model fitted to clustered data that have
  been `mashed' for fast processing (see \code{\link{mash}}); each
  cluster is a replicate sample: the function uses the vector of cluster
  frequencies (\eqn{n_j}) stored as an attribute of the mashed
  \code{capthist} by \code{mash}.

  \code{derived.external} combines detection parameter estimates from a
  fitted model with a vector of frequencies \code{nj} from replicate
  sampling units configured as in \code{cluster}. Detectors in
  \code{cluster} are assumed to match those in the fitted model with
  respect to type and efficiency, but sampling duration
  (\code{noccasions}), spacing etc. may differ. The \code{mask} should
  match \code{cluster}; if \code{mask} is missing, one will be
  constructed using the \code{buffer} argument and defaults from
  \code{make.mask}.

  \code{derived.session} accepts a single fitted model that must span
  multiple sessions; each session is interpreted as a replicate sample.

  Spatial variance may be calculated assuming simple random sampling
  (\code{method = "SRS"}) or using the neighbourhood variance estimator
  recommended by Stevens and Olsen (2003) for generalized random
  tessellation stratified (GRTS) samples and implemented in package
  \pkg{spsurvey} (\code{method = "local"}). For `local' variance
  estimates, the centre of each replicate must be provided in \code{xy},
  except where centres may be inferred from the data. The options
  \code{method = "poisson"} and \code{method = "binomial"} use
  theoretical (model-based) variance rather than the empirical spatial
  variance.
   
  }

\value{

  A dataframe with one row and the columns --
\item{estimate }{Horvitz-Thompson-like estimate of population density}
\item{SE.estimate }{SE of density estimate}
\item{lcl}{lower 100(1--alpha)\% confidence limit}
\item{ucl}{upper 100(1--alpha)\% confidence limit}
\item{CVn }{relative SE of number observed (across sampling units) }
\item{CVa }{relative SE of effective sampling area }
\item{CVD }{relative SE of density estimate}

}
\references{

Buckland, S. T., Anderson, D. R., Burnham, K. P., Laake, J. L., Borchers,
D. L. and Thomas, L. (2001) \emph{Introduction to Distance Sampling:
Estimating Abundance of Biological Populations}. Oxford University
Press, Oxford.

Stevens, D. L. Jr and Olsen, A. R. (2003) Variance estimation for
spatially balanced samples of environmental resources.
\emph{Environmetrics} \bold{14}, 593--610.

}

\seealso{
\code{\link{derived}}, \code{\link{esa}} 
}

\note{

  In versions before 2.1, the functionality of \code{derived.nj} and
  \code{derived.session} was provided by \code{empirical.VarD}, which
  has been removed.
  
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
  
}

\examples{

## The `ovensong' data are pooled from 75 replicate positions of a
## 4-microphone array. The array positions are coded as the first 4
## digits of each sound identifier. The sound data are initially in the
## object `signalCH'. We first impose a 52.5 dB signal threshold as in
## Dawson & Efford (2009, J. Appl. Ecol. 46:1201--1209). The vector nj
## includes 33 positions at which no ovenbird was heard. The first and
## second columns of `temp' hold the estimated effective sampling area
## and its standard error.

signalCH.525 <- subset(signalCH, cutval = 52.5)
nonzero.counts <- table(substring(rownames(signalCH.525),1,4))
nj <- c(nonzero.counts, rep(0, 75 - length(nonzero.counts)))
temp <- derived(ovensong.model.1, se.esa = TRUE)
derived.nj(nj, temp["esa",1:2])

## The result is very close to that reported by Dawson & Efford
## from a 2-D Poisson model fitted by maximizing the full likelihood.

## If nj vector has length 1, a theoretical variance is used...
msk <- ovensong.model.1$mask
A <- nrow(msk) * attr(msk, "area")
derived.nj (sum(nj), temp["esa",1:2], method = "poisson")
derived.nj (sum(nj), temp["esa",1:2], method = "binomial", area = A)

\dontrun{

## Set up an array of small (4 x 4) grids,
## simulate a Poisson-distributed population,
## sample from it, plot, and fit a model.
## mash() condenses clusters to a single cluster

testregion <- data.frame(x = c(0,2000,2000,0),
    y = c(0,0,2000,2000))
t4 <- make.grid(nx = 4, ny = 4, spacing = 40)
t4.16 <- make.systematic (n = 16, cluster = t4,
    region = testregion)
popn1 <- sim.popn (D = 5, core = testregion,
    buffer = 0)
capt1 <- sim.capthist(t4.16, popn = popn1)
fit1 <- secr.fit(mash(capt1), CL = TRUE, trace = FALSE)

## Visualize sampling
tempmask <- make.mask(t4.16, spacing = 10, type =
    "clusterbuffer")
plot(tempmask)
plot(t4.16, add = TRUE)
plot(capt1, add = TRUE)

## Compare model-based and empirical variances.
## Here the answers are similar because the data
## were simulated from a Poisson distribution,
## as assumed by \code{derived}

derived(fit1)
derived.mash(fit1)

## Now simulate a patchy distribution; note the
## larger (and more credible) SE from derived.mash().

popn2 <- sim.popn (D = 5, core = testregion, buffer = 0,
    model2D = "hills", details = list(hills = c(-2,3)))
capt2 <- sim.capthist(t4.16, popn = popn2)
fit2 <- secr.fit(mash(capt2), CL = TRUE, trace = FALSE)
derived(fit2)
derived.mash(fit2)

## The detection model we have fitted may be extrapolated to
## a more fine-grained systematic sample of points, with
## detectors operated on a single occasion at each...
## Total effort 400 x 1 = 400 detector-occasions, compared
## to 256 x 5 = 1280 detector-occasions for initial survey.

t1 <- make.grid(nx = 1, ny = 1)
t1.100 <- make.systematic (cluster = t1, spacing = 100,
    region = testregion)
capt2a <- sim.capthist(t1.100, popn = popn2, noccasions = 1)
## one way to get number of animals per point
nj <- attr(mash(capt2a), "n.mash")
derived.external (fit2, nj = nj, cluster = t1, buffer = 100,
    noccasions = 1)

## Review plots
library(MASS)
base.plot <- function() {
    eqscplot( testregion, axes = FALSE, xlab = "",
        ylab = "", type = "n")
    polygon(testregion)
}
par(mfrow = c(1,3), xpd = T, xaxs = "i", yaxs = "i")
base.plot()
plot(popn2, add = TRUE, col = "blue")
mtext(side=3, line=0.5, "Population", cex=0.8, col="black")
base.plot()
plot (capt2a, add = TRUE,title = "Extensive survey")
base.plot()
plot(capt2, add = TRUE, title = "Intensive survey")

}

}

\keyword{ models }