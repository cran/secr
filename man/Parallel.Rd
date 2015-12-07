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
  \code{\link{fxi.contour}} \tab individual \tab moderate-large \tab multiple individuals only \cr
  \code{\link{fxi.secr}} \tab individual \tab moderate-large \tab multiple individuals only \cr
  \code{\link{fx.total}} \tab individual \tab moderate-large \tab \cr
  \code{\link{pmixProfileLL}} \tab level of pmix \tab large \tab \cr
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
and the options \code{methods = FALSE}, \code{useXDR = .Platform$endian=='big'} (consistently from 2.10.0).
\code{methods = TRUE} is required for S4 methods, but as these are not used directly 
in \pkg{secr} it seems worth saving the time needed to load \pkg{methods}. \code{useXDR = FALSE} saves time ``where large amounts of data are to be transferred and all the nodes are little-endian'', which applies on Windows systems and some linux systems. Please report any issues.

}

\examples{

\dontrun{

sessionInfo()
# R version 3.2.0 (2015-04-16)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# quad-core i7 CPU, 16 Gb RAM
# ...

## almost no benefit with secr.fit

for (i in 1:5) 
    print(system.time(secr.fit(ovenCH, buffer = 400, trace = FALSE, ncores = i)))

#   user  system elapsed 
#  52.26    0.76   53.02 
#   user  system elapsed 
#   3.38    2.77   48.36 
#   user  system elapsed 
#   5.69    5.36   47.18 
#   user  system elapsed 
#   6.50    6.90   54.51 
#   user  system elapsed 
#   7.74    8.25   50.23 

## better for simulation, up to a point
 
for (i in 1:8)
    print(system.time(sim.secr(secrdemo.0, nsim = 20, tracelevel = 0, ncores = i)))

#   user  system elapsed 
# 262.93    3.74  266.71 
#   user  system elapsed 
#   0.68    0.02  139.61 
#   user  system elapsed 
#   0.69    0.01  110.05 
#   user  system elapsed 
#   0.70    0.00   84.24 
#   user  system elapsed 
#   0.70    0.02   69.45 
#   user  system elapsed 
#   0.67    0.06   75.67 
#   user  system elapsed 
#   0.75    0.00   70.00 
#   user  system elapsed 
#   0.73    0.05   75.71 
 
for (i in 1:8) 
    print(system.time(ip.secr (captdata, trace = FALSE, ncores = i)))

#   user  system elapsed 
# 132.59    0.00  132.58 
#   user  system elapsed 
#   0.25    0.05   71.08 
#   user  system elapsed 
#   0.20    0.04   50.76 
#   user  system elapsed 
#   0.32    0.07   42.13 
#   user  system elapsed 
#   0.26    0.09   41.74 
#   user  system elapsed 
#   0.39    0.08   39.01 
#   user  system elapsed 
#   0.25    0.09   36.52 
#   user  system elapsed 
#   0.46    0.08   35.91 
 
for (i in 1:8)
    print(system.time(LLsurface.secr(secrdemo.0, ncores = i)))
    
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#  48.00    0.40   48.41 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.00    0.03   27.94 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.00    0.03   20.60 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.03    0.02   17.63 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.01    0.05   16.41 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.02    0.04   16.10 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.05    0.03   16.03 
#Evaluating log likelihood across grid of 121 points...
#   user  system elapsed 
#   0.02    0.07   15.99 

}

}
