\name{mrl.plot}
\alias{mrl.plot}

\title{Mean Residual Life Plot}
\description{
  An empirical mean residual life plot, including confidence
  intervals, is produced. The mean residual life plot aids
  the selection of a threshold for the GPD or point process
  models.
}
\usage{
mrl.plot(data, umin = min(data), umax = max(data) - 0.1,
    conf = 0.95, nint = 100)
}
\arguments{
  \item{data}{A numeric vector of data to be fitted.}
  \item{umin, umax}{The minimum and maximum thresholds at which
    the mean residual life function is calculated.}
  \item{conf}{The confidence coefficient for the confidence
    intervals depicted in the plot.}
  \item{nint}{The number of points at which the mean residual
    life function is calculated.}
}

\seealso{\code{\link{gpd.fit}}, \code{\link{gpd.fitrange}},
  \code{\link{pp.fit}}}

\examples{
data(rain)
mrl.plot(rain)
}
\keyword{hplot}
