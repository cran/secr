\name{stoatDNA}
\docType{data}

\alias{stoatDNA}
\alias{stoatCH}
\alias{stoat.model.EX}
\alias{stoat.model.HN}
\alias{stoat.model.HZ}

\title{ Stoat DNA Data }

\description{

Data of A. E. Byrom from a study of stoats (\emph{Mustela erminea}) in
New Zealand. Individuals were identified from DNA in hair samples.

}

\usage{data(stoatDNA)}

\details{

The data are from a pilot study of stoats in red beech (\emph{Nothofagus fusca})
  forest in the Matakitaki Valley, South Island, New Zealand. Sticky
  hair-sampling tubes (\emph{n} = 94) were placed on a 3-km x 3-km grid
  with 500-m spacing between lines and 250-m spacing along lines. Tubes
  were baited with rabbit meat and checked daily for 7 days, starting on
  15 December 2001. Stoat hair samples were identified to individual
  using DNA microsatellites amplified by PCR from follicular tissue
  (Gleeson et al. 2010). Six loci were amplified and the mean number of
  alleles was 7.3 per locus. Not all loci could be amplified in 27\% of
  samples. A total of 40 hair samples were collected (Gleeson et
  al. 2010), but only 30 appear in this dataset; the rest presumably
  did not yield sufficient DNA for genotyping.

The data are provided as a single-session \code{capthist} object
`stoatCH'. Hair tubes are treated as `proximity' detectors which allow
an individual to be detected at multiple detectors on one occasion
(day), although there are no multiple detections in this dataset. Three
pre-fitted models are included: \code{stoat.model.HN},
\code{stoat.model.HZ}, and \code{stoat.model.EX} (with halfnormal,
hazard-rate and negative exponential detection functions,
respectively). Some components of the fitted secr models have been
removed with \code{\link{trim}} to save space.

\tabular{ll}{
  Object \tab Description \cr
  stoatCH \tab capthist object \cr
  stoat.model.EX \tab fitted secr model -- null, exponential detection
  function \cr
  stoat.model.HN \tab fitted secr model -- null, halfnormal detection
  function \cr
  stoat.model.HZ \tab fitted secr model -- null, hazard-rate detection
  function \cr
}

}

\source{

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\seealso{

  \code{\link{capthist}}, \code{\link{detection functions}},
  \code{\link{secr.fit}}

}

\references{

Gleeson, D. M., Byrom, A. E. and Howitt, R. L. J. (2010) Non-invasive
methods for genotyping of stoats (\emph{Mustela erminea}) in New
Zealand: potential for field applications. \emph{New Zealand Journal of
Ecology} \bold{34}, 356--359. Available on-line at
\url{http://www.newzealandecology.org}.

}

\note{

  The log-likelihood values reported for these data by \code{secr.fit}
  differ by a constant from those published by Efford et al. (2009)
  because the earlier version of DENSITY used in that analysis did not
  include the multinomial coefficient, which in this case is log(20!) or
  about +42.336. The previous analysis also used a coarser habitat mask
  than the default in \pkg{secr} (32 x 32 rather than 64 x 64) and this
  slightly alters the log-likelihood and \eqn{\Delta}{delta}AIC
  values. Fitting the hazard-rate detection function previously required
  the shape parameter z (or b) to be fixed, but the model can be fitted
  in \pkg{secr} without fixing z.

  Gleeson et al. (2010) address the question of whether there is enough
  variability at the sampled microsatellite loci to distinguish
  individuals. The reference to 98 sampling sites in that paper is a minor
  error (A. E. Byrom pers. comm.).

}

\examples{
summary(stoatCH)

\dontrun{
stoat.model.HN <- secr.fit(stoatCH, buffer = 1000, detectfn = 0)
stoat.model.HZ <- secr.fit(stoatCH, buffer = 1000, detectfn = 1)
stoat.model.EX <- secr.fit(stoatCH, buffer = 1000, detectfn = 2)
confint(stoat.model.HN, 'D')
## Profile likelihood interval(s)...
##         lcl        ucl
## D 0.01381419 0.04563511
}

## plot fitted detection functions
xv <- seq(0,800,10)
plot(stoat.model.EX, xval = xv, ylim = c(0,0.12), limits = FALSE,
    lty = 2)
plot(stoat.model.HN, xval = xv, limits = FALSE, lty = 1, add = TRUE)
plot(stoat.model.HZ, xval = xv, limits = FALSE, lty = 3, add = TRUE)

## review density estimates
collate(stoat.model.HZ, stoat.model.HN, stoat.model.EX,
    realnames='D', perm=c(2,3,4,1))
model.average(stoat.model.HZ, stoat.model.HN, stoat.model.EX,
    realnames='D')

}
\keyword{datasets}
