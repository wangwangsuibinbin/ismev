\name{gpd.fitrange}
\alias{gpd.fitrange}

\title{Fitting the GPD Model Over a Range of Thresholds}
\description{
  Maximum-likelihood fitting for a stationary GPD model,
  over a range of thresholds. Graphs of parameter estimates
  which aid the selection of a threshold are produced.
}
\usage{
gpd.fitrange(data, umin, umax, nint = 10, show = FALSE)
}
\arguments{
  \item{data}{A numeric vector of data to be fitted.}
  \item{umin, umax}{The minimum and maximum thresholds at which
    the model is fitted.}
  \item{nint}{The number of fitted models.}
  \item{show}{Logical; if \code{TRUE}, print details of each
    fit.}
}
\value{
  Two graphs showing maximum likelihood estimates and confidence
  intervals of the shape and modified scale parameters over a
  range of thresholds are produced.
}

\seealso{\code{\link{gpd.fit}}, \code{\link{mrl.plot}},
  \code{\link{pp.fit}}, \code{\link{pp.fitrange}}}

\examples{
\dontrun{data(rain)}
\dontrun{gpd.fitrange(rain, 10, 40)}
}
\keyword{hplot}
