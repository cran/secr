\name{secrdemo}
\alias{secrdemo}
\alias{secrdemo.0}
\alias{secrdemo.b}
\alias{secrdemo.CL}

\docType{data}

\title{ SECR Models Fitted to Demonstration Data }

\description{
  Objects of class \code{secr} formed by 

      \code{secrdemo.0 <- secr.fit (captdata)}

      \code{secrdemo.b <- secr.fit (captdata, model = list(g0 = ~b))}

      \code{secrdemo.CL <- secr.fit (captdata, CL = TRUE) }

where \code{captdata} is a \code{capthist} object using the
demonstration data (files 'trap.txt','capt.txt') from program Density
4.1 (Efford 2007).

}

\usage{data(secrdemo)}

\details{

The raw data are 235 fictional captures of 76 animals in 100
single-catch traps over 5 occasions. The fitted models use a halfnormal
detection function and the likelihood for multi-catch traps (expect
estimates of g0 to be biased because of trap saturation Efford et al.
2009). The first is a null model (i.e. parameters constant) and the
second fits a learned trap response.

}

\source{ Efford, M.G. (2007) Density 4.1: software for spatially
explicit capture-recapture. Department of Zoology, University of Otago,
Dunedin, New Zealand. \url{http://www.otago.ac.nz/density}.

Efford, M. G., Borchers D. L. and Byrom, A. E. (2009) Density estimation
by spatially explicit capture-recapture: likelihood-based methods. In:
D. L. Thomson, E. G. Cooch and M. J. Conroy (eds) \emph{Modeling
Demographic Processes in Marked Populations}. Springer, New York. Pp.
255--269.

}

\seealso{ \code{\link{captdata}}, \code{\link{capthist}} }

\examples{
## display the default fit with the print method for secr
data (secrdemo)
secrdemo.0
}

\keyword{datasets}
