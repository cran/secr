\name{details}
\alias{details}

\title{ Detail Specification for secr.fit }

\description{

The function \code{secr.fit} allows many options. Some of these are used
infrequently and have been bundled as a single argument \code{details}
to simplify the documentation. They are described here.
}

\section{Detail components}{

\code{details$centred} = TRUE causes coordinates of both traps and mask
to be centred on the centroid of the traps, computed separately for each
session in the case of multi-session data. This may be necessary to
overcome numerical problems when x- or y-coordinates are large
numbers. The default is not to centre coordinates.

\code{details$distribution} specifies the distribution of the number of
individuals detected; this may be conditional on the number in the
masked area ("binomial") or unconditional ("poisson").
\code{distribution} affects the sampling variance of the estimated
density. The default is "poisson". The component `distribution' may also
take a numeric value larger than nrow(capthist), rather than "binomial"
or "poisson". The likelihood then treats n as a binomial draw from a
superpopulation of this size, with consequences for the variance of
density estimates. This can help to reconcile MLE with Bayesian
estimates using data augmentation.

\code{details$fixedbeta} may be used to fix values of beta
parameters. It should be a numeric vector of length equal to the total
number of beta parameters (coefficients) in the model. Parameters to be
estimated are indicated by NA. Other elements should be valid values on
the link scale and will be substituted during likelihood
maximisation. Check the order of beta parameters in a previously fitted
model.

\code{details$hessian} is a character string controlling the computation
of the Hessian matrix from which variances and covariances are obtained.
Options are "none" (no variances), "auto" (the default) or "fdhess" (use
the function fdHess in \pkg{nlme}).  If "auto" then the Hessian from the
optimisation function is used. See also method = "none" below.

\code{details$ignoreusage} = TRUE causes the function to ignore
usage (varying effort) information in the traps component. The default
(\code{details$ignoreusage} = FALSE) is to include usage in the model.

\code{details$intwidth2} controls the half-width of the interval
searched by optimise() for the maximum likelihood when there is a single
parameter. Default 0.8 sets the search interval to \eqn{(0.2s, 1.8s)} where \eqn{s}
is the `start' value.

\code{details$LLonly} = TRUE causes the function to returns a single
evaluation of the log likelihood at the `start' values.

\code{details$param} chooses between various parameterisations of the
SECR model. The default \code{details$param = 0} is the formulation in
orchers and Efford (2008) and later papers.

\code{details$param = 1} selects the Gardner & Royle parameterisation of
the detection model (\emph{p0}, \eqn{\sigma}; Gardner et al. 2009) when
the detector type is `multi'. This parameterisation does not allow
detector covariates.

\code{details$param = 2} selects a currently unpublished
parameterisation in terms of (\eqn{esa(g_0, \sigma)}, \eqn{\sigma}).

\code{details$param = 3} selects a currently unpublished
parameterisation in terms of (\eqn{a_0(g_0, \sigma)}, \eqn{\sigma}).

\code{details$scaleg0} = TRUE causes g0 to be scaled by
\eqn{\mathrm{sigma}^{-2}}{sigma^-2}. Deprecated from 2.6.2: use a0
parameterisation instead.

\code{details$scalesigma} = TRUE causes sigma to be scaled by
\eqn{\mathrm{D}^{-0.5}}{1/sqrt(D)}.

\code{details$telemetrysigma} = TRUE uses coordinate information from telemetry
when capthist has attribute `xylist' (see \code{\link{addTelemetry}}).

}

\references{

  Gardner, B., Royle, J. A. and Wegan, M. T. (2009) Hierarchical models
for estimating density from DNA mark-recapture studies. \emph{Ecology}
\bold{90}, 1106--1115.

}

\seealso{
  \code{\link{secr.fit} }
}

\keyword{ models }
