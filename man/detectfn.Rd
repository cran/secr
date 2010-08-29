\name{detectfn}
\alias{detection functions}
\alias{detectfn}

\title{ Detection Functions }

\description{

A detection function relates the probability of detection to the distance
of a detector from a point. The reference point is usually thought of as
an animal's home-range centre. In \pkg{secr} only simple 2- or
3-parameter functions are used. Each type of function is identified by a
numeric code (see below). In most cases the name may also be used (as a
quoted string).

Some functions are defined only for simulation: these either cannot be
fitted by maximum likelihood (uniform) or have yet to be implemented
(compound halfnormal).

\tabular{llll}{
  Code \tab Name \tab Parameters \tab Function \cr
  0 \tab halfnormal \tab g0, sigma \tab \eqn{ g(d) = g_0 \exp \left(\frac{-d^2} {2\sigma^2} \right) }{ g(d) = g0 * exp(--d^2 / (2 sigma^2)) } \cr
  1 \tab hazard rate \tab g0, sigma, z \tab \eqn{ g(d) = g_0 [1 - \exp\{ {-(^d/_\sigma)^{-z}} \}] }{ g(d) = g0 * (1 -- exp(-- (d / sigma) ^(--z) )) } \cr
  2 \tab exponential  \tab g0, sigma \tab \eqn{ g(d) = g_0  \exp \{ -(^d/_\sigma) \} }{ g(d) = g0 * exp(-- d / sigma) } \cr
  3 \tab compound halfnormal  \tab g0, sigma, z \tab \eqn{ g(d) = 1 - (1 - g_0 \exp \left(\frac{-d^2} {2\sigma^2} \right) ^ z) }{ g(d) = 1 -- (1 -- g0 * exp(--d^2 / (2 sigma^2)))^z } \cr
  4 \tab uniform  \tab g0, sigma \tab \eqn{ g(d) = g_0, d <= \sigma; g(d) = 0, otherwise \} }{ g(d) = g0, d<=sigma; g(d) = 0, otherwise } \cr
  5 \tab w exponential  \tab g0, sigma, w \tab \eqn{ g(d) = g_0, d < w; g(d) = g_0  \exp \{ -(^(d-w)/_\sigma), otherwise \} }{ g(d) = g0, d < w; g(d) = g0 * exp(-- (d -- w) / sigma), otherwise } \cr
  6 \tab annular normal \tab g0, sigma, w \tab \eqn{ g(d) = g_0 \exp \left(\frac{-(d-w)^2} {2\sigma^2} \right) }{ g(d) = g0 * exp(--(d-w)^2 / (2 sigma^2)) } \cr
  9 \tab binary signal strength  \tab b0, b1 \tab \eqn{ g(d) = F(b_0 + b_1  d)) }{ g(d) = F(b0 + b1 * d) } \cr
  10 \tab signal strength  \tab beta0, beta1, sdS \tab \eqn{ g(d) = F((c - (\beta_0 + \beta_1  d)) / sdS) }{ g(d) = F((c -- (beta0 + beta1 * d)) / sdS) } \cr
  11 \tab signal strength spherical \tab beta0, beta1, sdS \tab \eqn{ g(d) = F((c - (\beta_0 + \beta_1  (d-1) - 10 * log ( d^2 ) / log(10))) / sdS) }{ g(d) = F((c -- (beta0 + beta1 * (d--1) -- 10 * log ( d^2 ) / log(10))) / sdS) } \cr
}

For functions (9), (10) and (11), 'F' is the standard normal
distribution function and 'c' is an arbitrary signal threshold. The two
parameters of (9) are functions of the parameters of (10) and (11):
\eqn{b_0 = (\beta_0 - c) / sdS}{b0 = (beta0-c)/sdS} and \eqn{b_1 =
\beta_1 / sdS}{b1 = beta1/sdS} (see Efford et al. 2009).

Function (11) includes an additional 'hard-wired' term for sound
attenuation due to spherical spreading. Detection probability at
distances less than 1 m is given by \eqn{g(d) = F((c - \beta_0) / sdS)}{
g(d) = F((c -- beta0) / sdS)}

The hazard-rate detection function was described by Hayes and Buckland
(1983). The compound halfnormal detection function follows Efford and
Dawson (2009). The signal strength and binary signal strength functions
are from Efford et al. (2009).

}

\references{

Efford, M. G. and Dawson, D. K. (2009) Effect of distance-related
heterogeneity on population size estimates from point counts. \emph{Auk}
\bold{126}, 100--111.

Efford, M. G., Dawson, D. K. and Borchers, D. L. (2009) Population
density estimated from locations of individuals on a passive detector
array. \emph{Ecology} \bold{90}, 2676--2682.

Hayes, R. J. and Buckland, S. T. (1983) Radial-distance models for the
line-transect method. \emph{Biometrics} \bold{39}, 29--42.

}

\seealso{ \code{\link{detectfnplot}}, \code{\link{secr detection models}} }

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\keyword{ models }

