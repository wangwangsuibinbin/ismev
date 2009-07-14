\name{pp.fitrange}
\alias{pp.fitrange}

\title{Fitting the Point Process Model Over a Range of Thresholds}
\description{
  Maximum-likelihood fitting for a stationary point process model,
  over a range of thresholds. Graphs of parameter estimates
  which aid the selection of a threshold are produced.
}
\usage{
pp.fitrange(data, umin, umax, npy = 365, nint = 10, show = FALSE, ...)
}
\arguments{
  \item{data}{A numeric vector of data to be fitted.}
  \item{umin, umax}{The minimum and maximum thresholds at which
    the model is fitted.}
  \item{npy}{The number of observations per year/block.}
  \item{nint}{The number of fitted models.}
  \item{show}{Logical; if \code{TRUE}, print details of each
    fit.}
   \item{...}{Optional arguments to \code{pp.fit}.}
}
\value{
  Three graphs showing maximum likelihood estimates and confidence
  intervals of the location, scale and shape parameters over a
  range of thresholds are produced.  A list object is returned invisibly
  with components:
   'threshold' numeric vector of length 'nint' giving the thresholds used,
  'mle' an 'nint X 3' matrix giving the maximum likelihood parameter estimates
        (columns are location, scale and shape respectively),
  'se' an 'nint X 3' matrix giving the estimated standard errors for the 
        parameter estimates (columns are location, scale and shape, resp.),
  'ci.low', 'ci.up' 'nint X 3' matrices giving the lower and upper 95% confidence
        intervals, resp. (columns same as for 'mle' and 'se').
}

\seealso{\code{\link{pp.fit}}, \code{\link{mrl.plot}},
  \code{\link{gpd.fit}}, \code{\link{gpd.fitrange}}}

\examples{
\dontrun{data(rain)}
\dontrun{pp.fitrange(rain, 10, 40)}
}
\keyword{hplot}
