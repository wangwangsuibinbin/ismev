\name{gpd.diag}
\alias{gpd.diag}

\title{Diagnostic Plots for GPD Models}
\description{
  Produces diagnostic plots for GPD models using the output
  of the function \code{gpd.fit}.
}
\usage{
gpd.diag(z)
}
\arguments{
  \item{z}{An object returned by \code{gpd.fit}.} 
}
\value{
  For stationary models four plots are produced; a probability
  plot, a quantile plot, a return level plot and a histogram
  of data with fitted density.

  For non-stationary models two plots are produced; a residual
  probability plot and a residual quantile plot.
}

\seealso{\code{\link{gpd.fit}}, \code{\link{gpd.prof}},
  \code{\link{pp.fit}}}

\examples{
data(rain)
rnfit <- gpd.fit(rain, 10)
gpd.diag(rnfit)
}
\keyword{hplot}
