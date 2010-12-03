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

Some functions such as (4) uniform are defined only for simulation as
the pose problems for maximum likelihood.

For function (7), `F' is the standard normal distribution function and
\eqn{\mu}{\mu} and \eqn{s}{s} are the mean and standard deviation on the
log scale of a latent variable representing a threshold of detection
distance. See Notes for the relationship to the fitted parameters sigma
and z.

For function (8), `G' is the cumulative distribution function of the
gamma distribution with shape parameter \emph{k} ( = \code{z}) and scale
parameter \eqn{\theta}{\theta} ( = \code{sigma/z}). See R's
\code{\link{pgamma}}.

For functions (9), (10) and (11), `F' is the standard normal
distribution function and `c' is an arbitrary signal threshold. The two
parameters of (9) are functions of the parameters of (10) and (11):
\eqn{b_0 = (\beta_0 - c) / sdS}{b0 = (beta0--c)/sdS} and \eqn{b_1 =
\beta_1 / s}{b1 = beta1/s} (see Efford et al. 2009).

Function (11) includes an additional `hard-wired' term for sound
attenuation due to spherical spreading. Detection probability at
distances less than 1 m is given by \eqn{g(d) = 1 - F \lbrace(c -
  \beta_0) / sdS \rbrace}{ g(d) = 1 -- F\{(c -- beta0) / sdS \} }

The hazard-rate detection function was described by Hayes and Buckland
(1983). The compound halfnormal detection function follows Efford and
Dawson (2009). The signal strength and binary signal strength functions
are from Efford et al. (2009).

\tabular{llll}{
  Code \tab Name \tab Parameters \tab Function \cr
  
  0 \tab halfnormal \tab g0, sigma \tab \eqn{ g(d) = g_0 \exp
  \left(\frac{-d^2} {2\sigma^2} \right) }{ g(d) = g0 * exp{--d^2 / (2
  sigma^2) } } \cr

  1 \tab hazard rate \tab g0, sigma, z \tab \eqn{ g(d) = g_0 [1 - \exp\{
  {-(^d/_\sigma)^{-z}} \}] }{ g(d) = g0 * (1 -- exp(-- (d / sigma)
  ^(--z) )) } \cr

  2 \tab exponential  \tab g0, sigma \tab \eqn{ g(d) = g_0  \exp \{
    -(^d/_\sigma) \} }{ g(d) = g0 * exp(-- d / sigma) } \cr
  
  3 \tab compound halfnormal  \tab g0, sigma, z \tab \eqn{ g(d) = g_0 [1
  - \{1 - \exp \left(\frac{-d^2} {2\sigma^2} \right)\} ^ z] }{ g(d) = g0
  * [1 -- \{1 -- exp(--d^2 / (2 sigma^2))]^z\} } \cr

  4 \tab uniform  \tab g0, sigma \tab \eqn{ g(d) = g_0, d <= \sigma;
  g(d) = 0, \mbox{otherwise}  }{ g(d) = g0, d<=sigma; g(d) = 0,
  otherwise } \cr

  5 \tab w exponential  \tab g0, sigma, w \tab \eqn{ g(d) = g_0, d < w;
  g(d) = g_0  \exp \left( -\frac{d-w}{\sigma} \right), \mbox{otherwise}
  }{ g(d) = g0, d < w; g(d) = g0 * exp(-- (d -- w) / sigma), otherwise }
 \cr
 
  6 \tab annular normal \tab g0, sigma, w \tab \eqn{ g(d) = g_0 \exp
  \lbrace \frac{-(d-w)^2} {2\sigma^2} \rbrace }{ g(d) = g0 *
  exp(--(d-w)^2 / (2 sigma^2)) } \cr

  7 \tab cumulative lognormal \tab g0, sigma, z \tab \eqn{ g(d) = g_0
  [ 1 - F \lbrace(d-\mu)/s \rbrace ] }{ g(d) = g0 [1 -- F\{(d--\mu)/s)\}] }
 \cr

  8 \tab cumulative gamma \tab g0, sigma, z \tab \eqn{ g(d) = g_0
  \lbrace 1 - G (d; k, \theta)\rbrace }{ g(d) = g0 \{ 1 -- G (d; k, \theta) \} }
 \cr

  9 \tab binary signal strength  \tab b0, b1 \tab \eqn{ g(d) = 1 - F
    \lbrace - ( b_0 + b_1  d) \rbrace }{ g(d) = 1 -- F \{-- (b0 + b1 * d) \} } \cr
  
  10 \tab signal strength  \tab beta0, beta1, sdS \tab \eqn{ g(d) =1 -
  F[\lbrace c - (\beta_0 + \beta_1  d) \rbrace / s] }{ g(d) = 1 -- F[ \{c -- (beta0 + beta1
  * d)\} / sdS] } \cr

11 \tab signal strength spherical \tab beta0, beta1, sdS \tab
    \eqn{ g(d) = 1 - F [ \lbrace c - (\beta_0 + \beta_1 (d-1) -  10 \log
      _{10}  d^2 ) \rbrace / s ]}{ g(d) = 1 -- F[\{c -- (beta0 + beta1 * (d--1) -- 10 * log10
   ( d^2 ) ) \} / sdS ] } \cr

}

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

\note{
    
  The parameters of function (7) are potentially confusing.  The fitted
parameters describe a latent threshold variable on the natural scale:
sigma (mean) = \eqn{\exp(\mu + s^2 / 2)}{exp(\mu + s^2/2)} and z
(standard deviation) = \eqn{\sqrt{\exp(s^2 +
2\mu)(\exp(s^2)-1)}}{sqrt(exp(s^2+2*\mu)*(exp(s^2)--1)))}. As with other
detection functions, sigma is a spatial scale parameter, although in
this case it corresponds to the mean of the threshold variable; the
standard deviation of the threshold variable (z) determines the shape
(roughly 1/max(slope)) of the detection function.  }

\keyword{ models }

