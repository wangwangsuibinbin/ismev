\name{gum.diag}
\alias{gum.diag}

\title{Diagnostic Plots for Gumbel Models}
\description{
  Produces diagnostic plots for Gumbel models using the output
  of the function \code{gum.fit}.
}
\usage{
gum.diag(z)
}
\arguments{
  \item{z}{An object returned by \code{gum.fit}.} 
}
\value{
  For stationary models four plots are produced; a probability
  plot, a quantile plot, a return level plot and a histogram
  of data with fitted density.

  For non-stationary models two plots are produced; a residual
  probability plot and a residual quantile plot.
}

\seealso{\code{\link{gev.fit}}, \code{\link{gum.fit}}}

\examples{
data(portpirie)
ppfit <- gum.fit(portpirie[,2])
gum.diag(ppfit)
}
\keyword{hplot}
