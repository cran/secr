\name{D.designdata}
\alias{D.designdata}

\title{ Construct Density Design Data }

\description{
  Internal function used by \code{\link{secr.fit}}, \code{\link{confint.secr}}, and \code{\link{score.test}}. 
}

\usage{
D.designdata (mask, Dmodel, grps, sessionlevels, sessioncov = NULL)
}

\arguments{
  \item{mask}{ \code{\link{mask}} object.}
  \item{Dmodel}{ formula for density model }
  \item{grps}{ vector of group names }
  \item{sessionlevels}{ vector of character values for session names }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
}
\details{

This is an internal \pkg{secr} function that you are unlikely ever to
use. Unlike \code{\link{secr.design.MS}}, this function does \emph{not}
call \code{model.matrix}.

}

\value{

Dataframe with one row for each combination of mask point, group and
session. The dataframe has an attribute 'dimD' that gives the relevant
dimensions: \code{ attr(dframe, 'dimD') = c(nmask, ngrp, R)}, where
\code{nmask} is the number of mask points, \code{ngrp} is the number of
groups, and \code{R} is the number of sessions. Columns correspond to
predictor variables in Dmodel.

}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\seealso{ \code{\link{secr.design.MS}} }

\keyword{ manip }
