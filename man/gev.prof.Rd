\name{gev.prof}
\alias{gev.prof}
\alias{gev.profxi}

\title{Profile Log-likelihoods for Stationary GEV Models}
\description{
  Produce profile log-likelihoods for shape parameters and
  m year/block return levels for stationary GEV models using
  the output of the function \code{gev.fit}.
}
\usage{
gev.prof(z, m, xlow, xup, conf = 0.95, nint = 100)
gev.profxi(z, xlow, xup, conf = 0.95, nint = 100)
}
\arguments{
  \item{z}{An object returned by \code{gev.fit}. The object
    should represent a stationary model.}
  \item{m}{The return level (i.e.\ the profile likelihood is for
    the value that is exceeded with probability 1/\code{m}).}
  \item{xlow, xup}{The least and greatest value at which to
    evaluate the profile likelihood.}
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

\seealso{\code{\link{gev.fit}}, \code{\link{gev.diag}}}

\examples{
data(portpirie)
ppfit <- gev.fit(portpirie[,2])
\dontrun{gev.prof(ppfit, m = 10, 4.1, 5)}
\dontrun{gev.profxi(ppfit, -0.3, 0.3)}
}
\keyword{hplot}
