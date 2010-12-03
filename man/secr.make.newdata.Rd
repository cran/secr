\name{secr.make.newdata}
\Rdversion{1.1}
\alias{secr.make.newdata}

\title{
Create Default Design Data
}

\description{
Generate a dataframe containing design data for the base levels of all predictors in an secr object.
}

\usage{
secr.make.newdata(object)
}

\arguments{
  \item{object}{fitted secr model object}
}

\details{

\code{secr.make.newdata} is used by \code{predict} in lieu of
user-specified 'newdata'. There is seldom any need to call
\code{secr.make.newdata} directly.

}

\value{

A dataframe with one row for each session and group, and columns for the
predictors used by \code{object$model}.

}

\seealso{

\code{\link{predict.secr}}, \code{\link{secr.fit}}

}

\examples{

## from previously fitted model
secr.make.newdata(secrdemo.b)

}

\keyword{ models }

