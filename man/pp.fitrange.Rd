\name{pp.fitrange}
\alias{pp.fitrange}

\title{Fitting the Point Process Model Over a Range of Thresholds}
\description{
  Maximum-likelihood fitting for a stationary point process model,
  over a range of thresholds. Graphs of parameter estimates
  which aid the selection of a threshold are produced.
}
\usage{
pp.fitrange(data, umin, umax, npy = 365, nint = 10, show = FALSE)
}
\arguments{
  \item{data}{A numeric vector of data to be fitted.}
  \item{umin, umax}{The minimum and maximum thresholds at which
    the model is fitted.}
  \item{npy}{The number of observations per year/block.}
  \item{nint}{The number of fitted models.}
  \item{show}{Logical; if \code{TRUE}, print details of each
    fit.}
}
\value{
  Three graphs showing maximum likelihood estimates and confidence
  intervals of the location, scale and shape parameters over a
  range of thresholds are produced.
}

\seealso{\code{\link{pp.fit}}, \code{\link{mrl.plot}},
  \code{\link{gpd.fit}}, \code{\link{gpd.fitrange}}}

\examples{
\dontrun{data(rain)}
\dontrun{pp.fitrange(rain, 10, 40)}
}
\keyword{hplot}
