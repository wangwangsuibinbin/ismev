\name{venice}
\alias{venice}
\title{Venice Sea Levels}
\usage{data(venice)}
\description{
  The \code{venice} data frame has 51 rows and 11 columns.
  The final ten columns contain the 10 largest sea levels
  observed within the year given by the first column.
  The ten largest sea levels are given for every year in
  the period 1931 to 1981, excluding 1935 in which only
  the six largest measurements are available.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Year}{A numeric vector of years.}
    \item{r1}{Annual sea level maxima.}
    \item{r2}{The second largest sea level.}
    \item{r3}{The third largest sea level.}
    \item{r4}{The forth largest sea level.}
    \item{r5}{The fifth largest sea level.}
    \item{r6}{The sixth largest sea level.}
    \item{r7}{The seventh largest sea level.}
    \item{r8}{The eigth largest sea level.}
    \item{r9}{The ninth largest sea level.}
    \item{r10}{The tenth largest sea level.}
  }
}
\source{
  Smith, R. L. (1986)
  Extreme value theory based on the \emph{r} largest annual
  events.
  \emph{Journal of Hydrology} \bold{86}, 27--43.
}
\references{
  Coles, S. G. (2001)
  \emph{An Introduction to Statistical Modelling of
    Extreme Values.} London: Springer.
}

\keyword{datasets}
