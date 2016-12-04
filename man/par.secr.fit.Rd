\name{par.secr.fit}
\alias{par.secr.fit}
\alias{par.derived}
\alias{par.region.N}

\title{Fit Multiple SECR Models}

\description{
  
  These functions are wrappers for \code{\link{secr.fit}}, \code{\link{region.N}} and
  \code{\link{derived}} that allow multiple models to be fitted
  simultaneously on multi-core machines.

}

\usage{

par.secr.fit (arglist, ncores = 1, seed = 123, trace = TRUE, logfile = "logfile.txt")
par.derived (secrlist, ncores = 1, ...)
par.region.N (secrlist, ncores = 1, ...)

}

\arguments{
  
  \item{arglist}{list of argument lists for \code{secr.fit} or a character
    vector naming such lists}
  
  \item{ncores}{ integer number of cores to be used for parallel processing}

  \item{seed}{integer pseudorandom number seed}

  \item{trace}{logical; if TRUE intermediate output may be logged}
  
  \item{logfile}{character name of file to log progress reports}

  \item{\dots}{other arguments passed to \code{derived} or \code{region.N}}

  \item{secrlist}{secrlist object}
  
}

\details{

  Any attempt in \code{arglist} to set \code{ncores > 1} for a
  particular secr fit is ignored.

  \code{trace} overrides any settings in \code{arglist}. Reporting of
  intermediate results is unreliable on Windows when \code{ncores > 1}.

  It is convenient to provide the names of the capthist and mask arguments
  in each component of arglist as character values (i.e. in quotes); objects thus
  named are exported from the workspace to each worker process (see Examples).

}

\value{

For \code{par.secr.fit} - \code{\link{secrlist}} of model fits (see
\code{\link{secr.fit}}). Names are created by prefixing ``fit.'' to the
names of \code{argslist}. If \code{trace} is TRUE then the total
execution time and finish time are displayed.

For \code{par.derived} - a list of dataframes output from
\code{\link{derived}}, applied to each model in turn.

For \code{par.region.N} - a list of dataframes output from
\code{\link{region.N}}, applied to each model in turn.

}

\seealso{
  \code{\link{secr.fit}},
  \code{\link{region.N}},
  \code{\link{derived}},
  \code{\link{AIC.secr}},
  \code{\link{predict.secr}},
  \link{Parallel}
}

\examples{
\dontrun{

fit0 <- list(capthist = 'captdata', model = g0~1) 
fitb <- list(capthist = 'captdata', model = g0~b)
fits <- par.secr.fit (c('fit0','fitb'), ncores = 2)
AIC(fits)

par.derived(fits, ncores = 2, se.esa = FALSE)

par.region.N(fits, ncores = 2)

}
}
\keyword{ model }

