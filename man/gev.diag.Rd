\name{gev.diag}
\alias{gev.diag}

\title{Diagnostic Plots for GEV Models}
\description{
  Produces diagnostic plots for GEV models using the output
  of the function \code{gev.fit}.
}
\usage{
gev.diag(z)
}
\arguments{
  \item{z}{An object returned by \code{gev.fit}.} 
}
\value{
  For stationary models four plots are produced; a probability
  plot, a quantile plot, a return level plot and a histogram
  of data with fitted density.

  For non-stationary models two plots are produced; a residual
  probability plot and a residual quantile plot.
}

\seealso{\code{\link{gev.fit}}, \code{\link{gev.prof}}}

\examples{
data(portpirie)
ppfit <- gev.fit(portpirie[,2])
gev.diag(ppfit)
}
\keyword{hplot}
