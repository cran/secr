\name{secr.design.MS}
\alias{secr.design.MS}
\title{ Construct Detection Model Design Matrices and Lookups }
\description{
  Internal function used by \code{\link{secr.fit}}. 
}
\usage{
secr.design.MS(capthist, models, timecov = NULL, sessioncov = NULL, 
    groups = NULL, dframe = NULL, naive = FALSE, bygroup = FALSE,
    ...)
}
\arguments{
  \item{capthist}{ \code{capthist} object }
  \item{models}{ list of formulae for parameters of detection }
  \item{timecov}{ optional dataframe of values of time (occasion-specific) covariate(s). }
  \item{sessioncov}{ optional dataframe of values of session-specific covariate(s). }
  \item{groups}{ optional vector of one or more variables with which to form groups. Each element should be the name of a factor variable in the \code{covariates} attribute of \code{capthist}. }
  \item{dframe}{ optional data frame of design data for detection parameters }
  \item{naive}{ logical if TRUE then modelled detection probability is for a naive animal (not caught previously); if FALSE then detection probability is contingent on individual's history of detection }
  \item{bygroup}{ logical if TRUE then the individual dimension of the parameter matrix is automatically collapsed to one row per group; if FALSE then the full dimensionality is retained (one row per individual) }
  \item{\dots}{ other arguments passed to the \R function \code{\link{model.matrix}} }
}
\details{

This is an internal \pkg{secr} function that you are unlikely ever to
use. \dots may be used to pass \code{contrasts.arg} to
\code{model.matrix}. 

Each real parameter is notionally different for each unique combination
of session, individual, occasion and detector, i.e., for \eqn{R}
sessions, \eqn{n} individuals, \eqn{S} occasions and \eqn{K} detectors
there are \emph{potentially} \eqn{R \times n \times S \times K}{R x n x
S x K} different values. Actual models always predict a much reduced set
of distinct values, and the number of rows in the design matrix is
reduced correspondingly; a parameter index array allows these to
retrieved for any combination of session, individual, occasion and
detector.

}
\value{
  A list with the components 
  \item{designMatrices}{ list of reduced design matrices, one for each real detection parameter }
  \item{parameterTable}{ index to row of the reduced design matrix for each real detection parameter; dim(parameterTable) = c(uniquepar, np), where uniquepar is the number of unique combinations of paramater values (uniquepar < \eqn{RnSK}) and np is the number of parameters in the detection model. }
  \item{PIA}{ Parameter Index Array - index to row of parameterTable for a given session, animal, occasion and detector; dim(PIA) = c(R,n,S,K) }
}

\examples{
secr.design.MS (captdata, models = list(g0 = ~b))$designMatrices
secr.design.MS (captdata, models = list(g0 = ~b))$parameterTable
}

\seealso{ \code{\link{D.designdata}} }

\keyword{ manip }
