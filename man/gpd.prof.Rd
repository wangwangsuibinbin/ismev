\name{gpd.prof}
\alias{gpd.prof}
\alias{gpd.profxi}

\title{Profile Log-likelihoods for Stationary GPD Models}
\description{
  Produce profile log-likelihoods for shape parameters and
  m year/block return levels for stationary GPD models using
  the output of the function \code{gpd.fit}.
}
\usage{
gpd.prof(z, m, xlow, xup, npy = 365, conf = 0.95, nint = 100)
gpd.profxi(z, xlow, xup, conf = 0.95, nint = 100)
}
\arguments{
  \item{z}{An object returned by \code{gpd.fit}. The object
    should represent a stationary model.}
  \item{m}{The return level (i.e.\ the profile likelihood is for
    the value that is exceeded with probability 1/\code{m}).}
  \item{xlow, xup}{The least and greatest value at which to
    evaluate the profile likelihood.}
  \item{npy}{The number of observations per year.}
  \item{conf}{The confidence coefficient of the plotted profile
    confidence interval.}
  \item{nint}{The number of points at which the profile
    likelihood is evaluated.}
}
\value{
  A plot of the profile likelihood is produced, with a horizontal
  line representing a profile confidence interval with confidence
  coefficient \code{conf}.
}

\seealso{\code{\link{gpd.fit}}, \code{\link{gpd.diag}}}

\examples{
data(rain)
rnfit <- gpd.fit(rain, 10)
\dontrun{gpd.prof(rnfit, m = 10, 55, 75)}
\dontrun{gpd.profxi(rnfit, -0.02, 0.15)}
}
\keyword{hplot}
