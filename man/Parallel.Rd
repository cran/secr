\name{Parallel}
\alias{Parallel}
\alias{ncores}
\alias{Multi-core processing}

\title{ Multi-core Processing }

\description{

From version 4.0 \pkg{secr} uses multi-threading in C++ (package \pkg{RcppParallel}, Allaire et al. 2018) to speed likelihood evaluation and hence model fitting in \code{secr.fit}. Detection histories are distributed over threads; the number of threads is set by default to one less than the number of available cores.

Earlier versions of \pkg{secr} made more limited use of multiple cores (CPUs)
through the package \pkg{\link{parallel}}. This mechanism is still available in the functions listed here, but the speed gains are often small or even negative. Set \code{ncores > 1} in the function call to use multiple cores.

\tabular{llll}{
  Function \tab Unit \tab Benefit \tab Notes \cr
  \code{\link{ip.secr}} \tab replicate \tab large \tab  \cr
  \code{\link{par.secr.fit}} \tab model \tab none \tab  \cr
  \code{\link{par.derived}} \tab fitted model \tab none \tab  \cr
  \code{\link{par.region.N}} \tab fitted model \tab none \tab  \cr
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
should accept incoming connections". You may possibly get warnings from R 
about closing unused connections. These can safely be ignored.

Use \code{parallel::detectCores()} to get
an idea of how many cores are available on your machine; this may (in
Windows) include virtual cores over and above the number of physical
cores. See RShowDoc("parallel", package = "parallel") in core R for
explanation.

In \code{secr.fit} the output component `proctime' misrepresents the
elapsed processing time when multiple cores are used.

}

\section{Warning}{

It appears that multicore operations in \pkg{secr} using \pkg{parallel} may fail if the packages \pkg{snow} and \pkg{snowfall} are also loaded. The error message is obscure:

``Error in UseMethod("sendData") : 
  no applicable method for 'sendData' applied to an object of class "SOCK0node"''

}

\examples{

\dontrun{

sessionInfo()
# R version 3.5.1 (2018-07-02)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1
# quad-core i7 CPU, 16 Gb RAM
# ...

## benefit from multi-threading in secr.fit

for (i in 1:8) 
     print(system.time(secr.fit(ovenCH, buffer = 400, trace = FALSE, ncores = i)))
#   user  system elapsed 
#  39.41    0.41   40.18 
#   user  system elapsed 
#  39.76    0.23   22.20 
#   user  system elapsed 
#  35.71    0.33   14.83 
#   user  system elapsed 
#  39.37    0.34   13.50 
#   user  system elapsed 
#  30.04    0.28    8.89 
#   user  system elapsed 
#  42.31    0.49   11.16 
#   user  system elapsed 
#  42.77    0.31   10.59 
#   user  system elapsed 
#  43.01    0.43    9.97 

## and for simulation...
 
for (i in 1:8)
     print(system.time(sim.secr(secrdemo.0, nsim = 20, tracelevel = 0, ncores = i)))
#   user  system elapsed 
#  154.75    1.56  156.40 
#   user  system elapsed 
#   3.26    0.16   88.07 
#   user  system elapsed 
#   3.11    0.08   74.52 
#   user  system elapsed 
#   3.01    0.08   55.21 
#   user  system elapsed 
#   3.37    0.19   47.33 
#   user  system elapsed 
#   3.06    0.20   53.66 
#   user  system elapsed 
#   3.21    0.18   47.17 
#   user  system elapsed 
#   3.03    0.15   45.57  

for (i in 1:8) 
    print(system.time(ip.secr (captdata, trace = FALSE, ncores = i)))

#  user  system elapsed 
# 87.75    0.07   88.06 
#  user  system elapsed 
#  0.45    0.09   48.89 
#  user  system elapsed 
#  0.22    0.09   36.11 
#  user  system elapsed 
#  0.41    0.14   33.07 
#  user  system elapsed 
#  0.46    0.27   29.50 
#  user  system elapsed 
#  0.34    0.29   27.44 
#  user  system elapsed 
#  0.41    0.35   27.41 
#  user  system elapsed 
#  0.48    0.42   28.66  
 
for (i in 1:8)
    print(system.time(LLsurface(secrdemo.0, ncores = i)))
    
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   27.13    0.36   27.50 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   27.70    0.31   19.60 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.79    0.35   18.15 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.07    0.39   16.29 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   29.43    0.22   15.24 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   31.31    0.23   15.29 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   31.87    0.59   15.18 
# Evaluating log likelihood across grid of 121 points...
#    user  system elapsed 
#   32.62    0.49   15.49

}

}
