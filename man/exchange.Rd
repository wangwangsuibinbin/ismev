\name{exchange}
\alias{exchange}
\title{UK/US and UK/Canada Exchange Rates}
\usage{data(exchange)}
\description{
  The \code{exchange} data frame has 975 rows and 2 columns.
  The columns contain daily exchange rates; UK sterling
  against the US dollar (first column) and UK sterling
  against the Canadian dollar (second column).
  The rownames contain the corresponding dates in a character
  string with the format \code{"2000/05/26"}. This can be
  converted into a \code{POSIXct} or \code{POSIXlt} object
  using \code{\link{as.POSIXct}} or \code{\link{as.POSIXlt}}.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{USD.GBP}{US against UK exchange rate.}
    \item{CAD.GBP}{Canada against UK exchange rate.}
  }
}
\source{
  Coles, S. G. (2001)
  \emph{An Introduction to Statistical Modelling of
    Extreme Values.} London: Springer.
}

\keyword{datasets}
