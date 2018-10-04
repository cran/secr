\name{turnover}
\alias{turnover}

\title{ Specifying a Dynamic Population }

\description{

\code{sim.popn} can simulate a multi-session population with known between-session survival, recruitment and movement probabilities. The parameter settings to achieve this are passed to \code{sim.popn} in its `details' argument. Components of `details' that are relevant to turnover are described below; see \code{\link{sim.popn}} for others.

Multi-session populations are generated in \code{sim.popn} whenever its argument `nsessions' is greater than 1. If \code{details$lambda} remains NULL (the default) then the population for each successive session is generated de novo from the given density model (model2D, D etc.). If a value is specified for \code{details$lambda} then only the first population is generated de novo; remaining populations are generated iteratively with probabilistic mortality, recruitment and movement as described here.

}

\section{Turnover components of \code{sim.popn} details argument}{

\tabular{lll}{
Component \tab Description \tab Default \cr
\code{phi} \tab per capita survival rate \eqn{\phi} \tab 0.7 \cr
\code{survmodel} \tab probability model for number of survivors \tab ``binomial" \cr
\code{lambda} \tab finite rate of increase \eqn{\lambda = N_{t+1} / N_t} \tab none \cr
\code{recrmodel} \tab probability model for number of recruits \tab ``poisson" \cr
\code{superN} \tab optional superpopulation size for multinomial recruitment model\tab see below \cr
\code{sigma.m} \tab spatial scale of movement kernel \tab 0 \cr
\code{edgemethod} \tab treatment of animals that cross the boundary \tab ``wrap" \cr
\code{wrap} \tab deprecated in 3.1.6; use edgemethod \tab TRUE i.e. edgemethod = ``wrap" \cr
}

}

\section{Survival}{

Survival is usually thought of as a Bernoulli process (outcome 0 or 1 for each individual) so the number of survivors \eqn{S} is a binomial variable (survmodel = "binomial"). Another approach is to fix the proportion surviving, but this can be done exactly only when \eqn{\phi N} is an integer. A (slightly ad hoc) solution is to randomly choose between the two nearest integers with probability designed in the long term (over many sessions) to give the required \eqn{\phi} (survmodel = "discrete").

}

\section{Population growth and recruitment}{

Per capita recruitment (\eqn{f}) is the difference between lambda and phi (\eqn{f = \lambda - \phi}), which must be non-negative (phi > lambda causes an error). The number of recruits B is a random variable whose probability distribution is controlled by details$recrmodel:
\tabular{ll}{
Value \tab Probability model \cr
"constantN" \tab Exact replacement of animals that die (B = \eqn{N_t - S})\cr
"binomial" \tab Binomial number of recruits (B ~ bin(\eqn{N_t, f}) \cr
"poisson" \tab Poisson number of recruits (B ~ pois(\eqn{f N_t})) \cr
"discrete" \tab Minimum-variance number of recruits (see Survival) \cr
"multinomial" \tab The POPAN model, conditioning on superpopulation size (e.g., Schwarz and Arnason 1996)) \cr
}

In the case of binomial recruitment there is a maximum of one recruit per existing individual, so lambda <= (phi+1). Multinomial recruitment requires a value for the superpopulation size. This may be provided as the details component "superN". If not specified directly, a value is inferred by projecting a trial initial (simulated) population using the specified phi and lambda.

}

\section{Movement}{

Individuals may shift their home range centre between sessions. Movement probability is governed by a circular bivariate normal (Gaussian) kernel with scale sigma.m. 

sigma.m is also used to indicate two special cases:

sigma.m = 0 indicates no change in home range centres between sessions.

sigma.m < 0 causes a new home range centre for each animal to be drawn at random from within the buffered area.

If movement takes an animal across the boundary of the arena in \code{sim.popn} the component "edgemethod" comes into play. By default, locations are toroidally wrapped i.e. the animal re-joins the population on the opposing edge. Other options are ``clip'' (discard), ``stop'' (stop just inside the boundary), ``reflect'' (bounce off edges to the limit of the dispersal) and ``none" (allow centres outside the buffered area).

}

\seealso{
  \code{\link{sim.popn} }
}

\examples{

par (mfrow = c(2,3), mar = c(1,1,1,1))

## birth and death only
grid <- make.grid(nx = 7, ny = 4, detector = 'proximity', spacing = 10)
pop <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 0.8, phi = 0.6, sigma.m = 0))
sapply(pop, nrow)  ## how many individuals?
plot(pop)

## movement only
pop2 <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 1, phi = 1, sigma.m = 10, edgemethod = "wrap"))
pop3 <- sim.popn (Nbuffer = 100, core = grid, nsessions = 6,    
    details = list(lambda = 1, phi = 1, sigma.m = 10, edgemethod = "clip"))
sapply(pop2, nrow)  ## how many individuals?
plot(pop2)

## show effect of edgemethod --
tracks <- function (pop, ...) {
    rn <- unlist(lapply(pop, rownames))
    plot(pop[[1]], pch = 16)
    for (i in unique(rn))
       lines(t(sapply(pop, '[', i, 1:2)), ...)
}
par (mfrow=c(1,2))
tracks(pop2, type = 'o')
tracks(pop3, type = 'o')

}
\keyword{ datagen }
