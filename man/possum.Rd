\name{possum}
\docType{data}

\alias{possum}
\alias{possumCH}
\alias{possumarea}
\alias{possumregion}
\alias{possummask}
\alias{possum.model.1}
\alias{possum.model.1b}

\title{ Brushtail Possum Trapping Dataset }

\description{

Data from a trapping study of brushtail possums at Waitarere, North
Island, New Zealand.

}

\usage{data(possum)}

\details{

Brushtail possums (\emph{Trichosurus vulpecula}) are an unwanted
invasive species in New Zealand.  Although most abundant in forests,
where they occasionally exceed densities of 15 / ha, possums live
wherever there are palatable food plants and shelter.

Efford et al. (2005) reported a live-trapping study of possums in
\emph{Pinus radiata} plantation on coastal sand dunes.  The 294-ha site
at Waitarere in the North Island of New Zealand was a peninsula, bounded
on one side by the sea and on two other sides by the Manawatu river. 
Cage traps were set in groups of 36 at 20-m spacing around the perimeter
of five squares, each 180 m on a side.  The squares ('hollow grids') were
centred at random points within the 294-ha area.  Animals were tagged and
released daily for 5 days in April 2002.  Subsequently, leg-hold trapping
was conducted on a trapping web centred on each square (data not
reported here), and strenuous efforts were made to remove all possums by
cyanide poisoning and further leghold trapping across the entire area. 
This yielded a density estimate of 2.26 possums / ha.

Traps could catch at most one animal per day.  The live-trapped animals
comprised 46 adult females, 33 adult males, 10 immature females and 11
immature males; sex and/or age were not recorded for 4 individuals (M.
Coleman unpubl. data).  One female possum was twice captured at two sites
on one day, having entered a second trap after being released; one
record in each pair was selected arbitrarily and discarded.

The data are provided as a single-session \code{capthist} object
'possumCH'.  'possummask' is a matching mask object - see Examples.  Two
fitted models ('possum.model.1' & 'possum.model.1b') are provided for
illustration.

From secr 1.5 a dataframe \code{possumarea} is included to indicate the
boundary of potential habitat, and this is used to clip
\code{possummask}. \code{possumarea} comprises a single polygon
representing the extent of terrestrial vegetation to the west, north and
east, and an arbitrary straight southern boundary. The boundary is also
included as a shapefile and as a text file (`possumarea.shp' etc. and
`possumarea.txt' in the package `extdata' folder). See Examples in
\code{\link{make.mask}}.

% From secr 2.1 a further dataframe \code{possumregion} is included that
% equates to the nominal removal study area (Efford et al. 2005 Fig. 1).

Clipping of the mask, introduced in secr 1.5, slightly changes the
fitted models (and hence the density estimates).

\tabular{ll}{
  Object \tab Description \cr
  possumCH \tab capthist object \cr
  possummask \tab mask object \cr
  possumarea \tab habitat perimeter \cr
%  possumregion \tab nominal boundary of removal region \cr
  possum.model.1 \tab fitted secr model -- null \cr
  possum.model.1b \tab fitted secr model -- trap response g0 \cr
}

Some components of the fitted secr models have been removed with
\code{\link{trim}} to save space.

}

\source{

Landcare Research, New Zealand.

}

\references{

Borchers, D.L. and Efford, M.G. (2008) Spatially explicit maximum
likelihood methods for capture-recapture studies. \emph{Biometrics} \bold{64},
377--385.

Efford, M. G., Dawson, D. K. and Robbins C. S. (2004) DENSITY: software
for analysing capture-recapture data from passive detector arrays.
\emph{Animal Biodiversity and Conservation} \bold{27}, 217--228. 

Efford, M. G., Warburton, B., Coleman, M. C. and Barker, R. J. (2005) A
field test of two methods for density estimation. \emph{Wildlife Society
Bulletin} \bold{33}, 731--738.

}

\seealso{

\code{\link{capthist}}

}

\examples{

\dontrun{

## code used to construct dataset
setwd('d:\\density secr 2.1\\package data\\possum')
possumCH <- read.capthist('foxton.txt', 'foxtraps.txt', detector = 'single')
## drop within-day duplicates of tag 5861 in traps 11,12, to conform DENSITY
possumCH['5861',] <- c(134, 133, 134, 133, 135)
trps <- traps(possumCH)
covariates(traps(possumCH)) <- data.frame(grid =
    factor(substring(rownames(trps),1,2)))
possumarea <- read.table ('possumarea.txt', header = TRUE)
possummask <- make.mask(trps, buffer = 300, type='pdot', pdotmin =
    0.001, detectpar = list(g0=0.2, sigma=60), spacing = 10, poly =
    possumarea)
possum.model.1 <- secr.fit(possumCH, mask = possummask, trace=F)
possum.model.1b <- secr.fit(possumCH, mask = possummask, model =
    list(g0~b), trace=F)
save(possumCH, possummask, possumarea, possum.model.1, possum.model.1b,
    file='d:\\density secr 2.0\\secr\\data\\possum.RData')
}

plot(possummask)
plot(possumCH, tracks = TRUE, add = TRUE)
plot(traps(possumCH), add = TRUE)
lines(possumarea)

summary(possumCH)

## compare & average pre-fitted models
AIC(possum.model.1, possum.model.1b)
model.average(possum.model.1, possum.model.1b)

}

\keyword{datasets}
