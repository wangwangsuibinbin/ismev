\name{pp.diag}
\alias{pp.diag}

\title{Diagnostic Plots for Point Process Models}
\description{
  Produces diagnostic plots for point process models using the
  output of the function \code{pp.fit}.
}
\usage{
pp.diag(z)
}
\arguments{
  \item{z}{An object returned by \code{pp.fit}.} 
}
\value{
  For stationary models two plots are produced; a probability
  plot and a quantile plot.

  For non-stationary models two plots are produced; a residual
  probability plot and a residual quantile plot.
}

\seealso{\code{\link{pp.fit}}, \code{\link{gpd.fit}}}

\examples{
data(rain)
rnfit <- pp.fit(rain, 10)
pp.diag(rnfit)
}
\keyword{hplot}
