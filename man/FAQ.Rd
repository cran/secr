\name{FAQ}
\alias{FAQ}
\alias{secr FAQ}

\title{ Frequently Asked Questions, And Others }

\description{
A place for hints and miscellaneous advice.
}

\section{How do I install secr?}{

Under Windows, it is simplest to install the package binary from the
Rgui. Save the file secr_1.3.x.zip to a local folder and use the menu
option "Packages | Install packages(s) from local zip files...".

For other systems (see below) you will need to install the source
package secr_1.3.x.tar.gz. 

Whatever your system, you also need to get the package \pkg{abind} (use
Packages | Install package(s)... to download from CRAN). Other required
packages (\pkg{MASS}, \pkg{nlme}, \pkg{stats}) should be available as
part of your \R installation.

}

\section{ Does secr work on operating systems other than Windows?}{

Yes, at least in theory. From version 1.2.8 you should be able to
install the source package secr\_1.x.x.tar.gz if you have the necessary
tools (check the 'R Installation and Administration' manual or the 'R
for Mac OS X FAQ' on CRAN for help on installing source packages).
Thanks to Ray Brownrigg and Joanne Potts for vital assistance. I
would be glad to hear of successful installations on non-Windows systems
(particularly if anyone produces a binary for Mac OS).

}

\section{ How can I get help?}{

There are three general ways of displaying documentation from within \R.
Firstly, you can bring up help pages for particular functions from the
command prompt.  For example:

\code{? secr.fit}
 
Secondly, help.search() lets you ask for a list of the help pages on a
vague topic. From \R version 2.8.0 you can use just ??.  For example:

\code{?? 'linear models'}

Thirdly, you can display various documents:

\code{RShowDoc ('secr-manual', package='secr')}

\code{RShowDoc ('secr-overview', package='secr')}

See below for more \R tips.

}


\section{ How should I report a problem? }{

If you get really stuck or find something you think is a bug then please
report the problem. For the moment, the correct address is
\email{density.software@otago.ac.nz}. 

When reporting problems please send an actual dataset (ideally, the
simplest one that exhibits the problem). Use \code{\link{save}} to wrap
several \R objects together in one .RData file, e.g., \code{save('captdata',
'secrdemo.0', 'secrdemo.b', file = 'mydata.RData')}. Also, paste the
output from \code{packageDescription("secr")} into the text of your
message.

Needless to say, we cannot promise to solve all problems.

}

\section{ Why do I get different answers from secr and Density? }{

Strictly speaking, this should not happen if you have specified the same
model and likelihood, although you may see a little variation due to the
different maximization algorithms. Likelihoods (and estimates) may
differ if you use different integration meshes (habitat masks), which
can easily happen because the programs differ in how they set up the
mesh.  If you want to make a precise comparison, save the Density mesh to a file and
read it into \pkg{secr}, or vice versa.

Extreme data, especially rare long-distance movements, may be handled
differently by the two programs. The 'minprob' component of the
'details' argument of \code{secr.fit} sets a threshold of probability
for capture histories (smaller values are all set to minprob), whereas
Density has no explicit limit. In the current version the default
minprob has been reduced from 1e-20 to 1e-50. If you find a discrepancy
with Density it may be worth lowering minprob even further.

}

\section{ How can I speed up model fitting and model selection?}{

If you don't need to model variation in density over space or time then
consider maximizing the conditional likelihood in \code{secr.fit} (CL =
TRUE). This reduces the complexity of the optimization problem,
especially where there are several sessions and you want
session-specific density estimates (by default, \code{derived} returns a
separate estimate for each session even if the detection parameters are
constant across sessions).

Check the extent and spacing of the habitat mask that you are using.
Execution time is roughly proportional to the number of mask points.
Default settings can lead to very large masks for detector arrays that
are elongated 'north-south' because the number of points in the
east-west direction is fixed. Compare results with a much sparser mask
(e.g., nx=32 instead of nx=64).

Do you really need to fit that complex model? Chasing down small
decrements in AIC is so last-century. Remember that detection parameters
are mostly nuisance parameters, and models with big differences in AIC
may barely differ in their density estimates. This is a good topic for
further research - we seem to need a 'focussed information criterion'
(Claeskens and Hjort 2008) to discern the differences that matter. 

Use \code{score.test} to compare nested models. At each stage this
requires only the more simple model to have been fitted in full; further
processing is required to obtain a numerical estimate of the gradient of
the likelihood surface for the more complex model, but this is much
faster than maximizing the likelihood.

}

\section{ Will secr soon include feature X?}{

\describe{

\item{Finite mixture models}{Possibly. Mixture models do seem
desirable as an option to allow for unmodelled heterogeneity. However,
mixtures add an extra level of complexity across the entire detection
model structure, and an additional 'real' parameter (for the mixing
proportion) so the decision is not trivial.  A decision was made to
release and evaluate \pkg{secr} 1.3 before adding even more complexity.}

\item{Inverse prediction}{Yes, from version 1.2.7.  Fitting SECR models by
simulation and inverse prediction is still the only rigorous method for
dealing with data from single-catch traps.  However, maximizing the
likelihood as if the data came from multi-catch traps yields fairly good
density estimates even with single-catch data (estimates of detection g0
are biased).  Limited capability for simulation-based fitting is provided
by the function \code{ip.secr}.  It is not intended to extend
\code{ip.secr} to include a wider range of models. }

}

}

\section{Things You Might Need To Know About \R}{

The function \code{findFn} in package \pkg{sos} lets you search CRAN for
R functions by matching text in their documentation.

There is now a vast amount of \R advice available on the web. For the
terminally frustrated, 'R inferno' by Patrick Burns is recommended
(\url{www.burns-stat.com/pages/Tutor/R_inferno.pdf}). "If you are using
R and you think you're in hell, this is a map for you".

Method functions for S3 classes cannot be listed in the usual way by
typing the function name at the \R prompt because they are 'hidden' in a
namespace. Get around this with getAnywhere(). For example: 

\code{getAnywhere(print.secr)}

R objects have 'attributes' that usually are kept out of sight.
Important attributes are 'class' (all objects), 'dim' (matrices and
arrays) and 'names' (lists). \pkg{secr} hides quite a lot of useful data
as named 'attributes'. Usually you will use summary and extraction
methods (\code{traps}, \code{covariates}, \code{usage} etc.) to view and change
the attributes of the various classes of object in \pkg{secr}. If you're
curious, you can reveal the lot with 'attributes'.  For example:

\code{data(captdata)}

\code{traps(captdata)       ## extraction method for 'traps'}

\code{attributes(captdata)  ## all attributes}

}

\references{
Claeskens, G. and Hjort N. L. (2008) \emph{Model Selection and Model Averaging}. Cambridge: Cambridge University Press.
}

\author{ Murray Efford \email{murray.efford@otago.ac.nz} }

\keyword{ manip }

