\name{dowjones}
\alias{dowjones}
\title{Daily Closing Prices of The Dow Jones Index}
\usage{data(dowjones)}
\description{
  The \code{dowjones} data frame has 1304 rows and 2 columns.
  The second column contains daily closing prices of the Dow
  Jones Index over the period 1996 to 2000. The first column
  contains a \code{\link{POSIXct}} object giving the dates of
  each observation.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Date}{A \code{\link{POSIXct}} object containing dates.}
    \item{Index}{A numeric vector containing daily closing
      prices of the Dow Jones Index.}
  }
}
\source{
  Coles, S. G. (2001)
  \emph{An Introduction to Statistical Modelling of
    Extreme Values.} London: Springer.
}
\keyword{datasets}
