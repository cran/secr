\name{Parallel}
\alias{Parallel}
\alias{ncores}
\alias{Multi-core processing}

\title{ Multi-core Processing }

\description{

From version 2.4.0 \pkg{secr} makes limited use of multiple cores (CPUs)
through the package \pkg{\link{parallel}}. Only the few \pkg{secr}
functions listed below make any use of parallel processing. Increased
speed can be expected with \code{sim.secr} (e.g., x3 with 4 cores), but
gains in \code{secr.fit} are much smaller and may be negative.

The function \code{\link{par.secr.fit}} added in 2.7.1 is useful if you have a set of
different models you wish to fit: if you have as many cores as models
then all will complete in about the time required to run the longest. 

\tabular{llll}{
  Function \tab Unit \tab Benefit \tab Notes \cr
  \code{\link{secr.fit}} \tab session likelihood \tab small-moderate \tab multi-session models only \cr
  \code{\link{score.test}} \tab model \tab moderate \tab multi-model comparisons only \cr
  \code{\link{derived}} \tab session \tab moderate \tab SE by parameter if
  one session \cr
  \code{\link{mask.check}} \tab spacing x buffer \tab moderate-large \tab no file output, suppresses messages \cr
  \code{\link{sim.secr}} \tab replicate \tab large \tab all models, suppresses messages \cr
  \code{\link{ip.secr}} \tab replicate \tab large \tab  \cr
  \code{\link{LLsurface.secr}} \tab parameter combination \tab large \tab  \cr
  \code{\link{par.secr.fit}} \tab model \tab large \tab  \cr
  \code{\link{par.derived}} \tab fitted model \tab large \tab  \cr
  \code{\link{par.region.N}} \tab fitted model \tab large \tab  \cr
}
`Unit' refers to the unit of work sent to each worker process. As a guide, a `large' benefit means >60\% reduction in process time with
  4 CPUs.

\pkg{parallel} offers several different mechanisms, bringing together
the functionality of \pkg{multicore} and \pkg{snow}. The mechanism used
by \pkg{secr} is the simplest available, and is expected to work across all
operating systems. Technically, it relies on Rscript and communication
between the master and worker processes is \emph{via} sockets. As stated
in the \pkg{parallel} documentation "Users of Windows and Mac OS X may
expect pop-up dialog boxes from the firewall asking if an R process
should accept incoming connections".

To use multiple cores, install \pkg{parallel} from CRAN and set
\code{ncores > 1} in the function call. Use \code{detectCores()} to get
an idea of how many cores are available on your machine; this may (in
Windows) include virtual cores over and above the number of physical
cores. See RShowDoc("parallel", package = "parallel") in core R for
explanation.

You may possibly get warnings from R about closing unused
connections. These can safely be ignored.

In \code{sim.secr}, new datasets are generated in the master process, so
there is no need to manage the random number streams in the worker
processes.

In \code{secr.fit} the output component `proctime' misrepresents the
elapsed processing time when multiple cores are used.

Worker processes are created in \code{secr.fit} with \code{makeCluster}
and the options \code{methods = FALSE}, \code{useXDR = FALSE}. This has been
tested only on Windows systems.

The code used internally by \pkg{secr} is quite simple and could be
adapted as a wrapper for user-defined simulations. See Examples.

}

\examples{

\dontrun{

# R version 2.15.2 (2012-10-26)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# quad-core i7 CPU

library(parallel)
detectCores()
# [1] 8

# ovenCH is a 5-session dataset
# multiple cores help a little here

system.time(f5 <- secr.fit(ovenCH, buffer = 400, trace = FALSE, ncores = 1))
#   user  system elapsed 
#  61.21    0.95   62.28 
system.time(f5 <- secr.fit(ovenCH, buffer = 400, trace = FALSE, ncores = 5))
#   user  system elapsed 
#   8.51    8.66   35.81 

# however, there is substantial benefit when simulating

system.time(s1 <- sim.secr(secrdemo.0, nsim = 20))
#   user  system elapsed 
# 789.90    4.41  795.07 
system.time(s4 <- sim.secr(secrdemo.0, nsim = 20, ncores = 4))
#   user  system elapsed 
#  26.91    0.34  276.15 

system.time(ip.secr (captdata, ncores = 1))
#   user  system elapsed 
# 149.93    0.01  150.72 
system.time(ip.secr (captdata, ncores = 6))
#  user  system elapsed 
#   0.47    0.19   41.06 

system.time(score.test (secrdemo.0, g0 ~ b, g0~t, g0 ~ bk, ncores = 1))
#   user  system elapsed 
# 130.73    0.45  131.34 
system.time(score.test (secrdemo.0, g0 ~ b, g0~t, g0 ~ bk, ncores = 4))
#   user  system elapsed 
#   0.04    0.01  109.69

system.time(derived(ovenbird.model.D, ncores=1))
#   user  system elapsed 
#   7.99    0.02    8.00 
system.time(derived(ovenbird.model.D, ncores=5))
#   user  system elapsed 
#   0.05    0.04    4.06 

system.time( LLsurface.secr(secrdemo.0, ncores=1))
#   user  system elapsed 
#  40.97    0.64   41.66 
system.time( LLsurface.secr(secrdemo.0, ncores=4))
#   user  system elapsed 
#   0.05    0.06   13.82 
system.time( LLsurface.secr(secrdemo.0, ncores=8))
#   user  system elapsed 
#   0.03    0.11   11.14 

## the code used in LLsurface.secr() looks like this

if (ncores > 1) {
    require(parallel)
    clust <- makeCluster(ncores)
    # load 'secr' in each worker process
    clusterEvalQ(clust, library(secr))
    # send data to each worker process from master process
    clusterExport(clust, c("object", "details"), environment())
    # run function LL for each row of matrix 'grid'
    # LL accepts one argument, a vector of parameter values
    # (can also use parLapply, parSapply etc.)
    temp <- parRapply (clust, grid, LL)
    stopCluster(clust)
}

}

}
