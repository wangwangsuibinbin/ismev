\name{rlarg.diag}
\alias{rlarg.diag}

\title{Diagnostic Plots for Order Statistics Models}
\description{
  Produces diagnostic plots for order statistics models using
  the output of the function \code{rlarg.fit}.
}
\usage{
rlarg.diag(z, n = z$r)
}
\arguments{
  \item{z}{An object returned by \code{rlarg.fit}.}
  \item{n}{Probability and quantile plots are produced for the
    largest \code{n} order statistics.}
}
\value{
  For stationary models four plots are initially produced;
  a probability plot, a quantile plot, a return level plot
  and a histogram of data with fitted density.
  Then probability and quantile plots are produced for the
  largest \code{n} order statistics.

  For non-stationary models residual probability plots and
  residual quantile plots are produced for the largest
  \code{n} order statistics.
}

\seealso{\code{\link{rlarg.fit}}}

\examples{
\dontrun{data(venice)}
\dontrun{venfit <- rlarg.fit(venice[,-1])}
\dontrun{rlarg.diag(venfit)}
}
\keyword{hplot}
