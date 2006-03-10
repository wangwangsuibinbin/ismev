\name{pp.fit}
\alias{pp.fit}

\title{Maximum-likelihood Fitting for the Point Process Model}
\description{
  Maximum-likelihood fitting for the point process model,
  including generalized linear modelling of each parameter.
}
\usage{
pp.fit(xdat, threshold, npy = 365, ydat = NULL, mul = NULL, sigl =
    NULL, shl = NULL, mulink = identity, siglink = identity, shlink =
    identity, show = TRUE, method = "Nelder-Mead", maxit = 10000, \dots)
}
\arguments{
  \item{xdat}{A numeric vector of data to be fitted.}
  \item{threshold}{The threshold; a single number or a numeric
    vector of the same length as \code{xdat}.}
  \item{npy}{The number of observations per year/block.}
  \item{ydat}{A matrix of covariates for generalized linear modelling
    of the parameters (or \code{NULL} (the default) for stationary
    fitting). The number of rows should be the same as the length
    of \code{xdat}.}
  \item{mul, sigl, shl}{Numeric vectors of integers, giving the columns
    of \code{ydat} that contain covariates for generalized linear
    modelling of the location, scale and shape parameters repectively
    (or \code{NULL} (the default) if the corresponding parameter is
    stationary).}
  \item{mulink, siglink, shlink}{Inverse link functions for generalized
    linear modelling of the location, scale and shape parameters
    repectively.}
  \item{show}{Logical; if \code{TRUE} (the default), print details of
    the fit.}
  \item{method}{The optimization method (see \code{\link{optim}} for
    details).}
  \item{maxit}{The maximum number of iterations.}
  \item{\dots}{Other control parameters for the optimization. These
    are passed to components of the \code{control} argument of
    \code{optim}.}
}
\details{
  For non-stationary fitting it is recommended that the covariates
  within the generalized linear models are (at least approximately)
  centered and scaled (i.e.\ the columns of \code{ydat} should be
  approximately centered and scaled).  Otherwise, the numerics may
  become unstable.

  As of version 1.32, a more accurate estimate of the exceedance
  rate, in the face of covariates, is used (at the expense of
  computational efficiency).  In particular, when including covariates,
  parameter estimates may differ from those in Coles (2001).
}
\value{
  A list containing the following components. A subset of these
  components are printed after the fit. If \code{show} is
  \code{TRUE}, then assuming that successful convergence is
  indicated, the components \code{nexc}, \code{nllh}, \code{mle}
  and \code{se} are always printed.
  
  \item{trans}{An logical indicator for a non-stationary fit.}
  \item{model}{A list with components \code{mul}, \code{sigl}
    and \code{shl}.}
  \item{link}{A character vector giving inverse link functions.}
  \item{threshold}{The threshold, or vector of thresholds.}
  \item{npy}{The number of observations per year/block.}
  \item{nexc}{The number of data points above the threshold.}
  \item{data}{The data that lie above the threshold. For
    non-stationary models, the data is standardized.}
  \item{conv}{The convergence code, taken from the list returned by
    \code{\link{optim}}. A zero indicates successful convergence.}
  \item{nllh}{The negative logarithm of the likelihood evaluated at
    the maximum likelihood estimates.}
  \item{vals}{A matrix with four columns containing the maximum
    likelihood estimates of the location, scale and shape
    parameters, and the threshold, at each data point.}
  \item{gpd}{A matrix with three rows containing the maximum
    likelihood estimates of corresponding GPD location, scale
    and shape parameters at each data point.}
  \item{mle}{A vector containing the maximum likelihood estimates.}
  \item{cov}{The covariance matrix.}
  \item{se}{A vector containing the standard errors.}
}

\references{  
Beirlant J, Goegebeur Y, Segers J and Teugels J. (2004). Statistics of Extremes, Wiley,
  Chichester, England.

Coles, Stuart (2001).  An Introduction to Statistical Modeling of
 Extreme Values.  Springer-Verlag, London. }
\note{ This code is adapted by Eric Gilleland from code originally written for S-Plus by Stuart Coles, and ported to R by Alec Stephenson. See
details section above.}

\section{Warning }{Different optimization methods may result in wildly different parameter estimates.}

\seealso{\code{\link{pp.diag}}, \code{\link{optim}},
  \code{\link{pp.fitrange}}, \code{\link{mrl.plot}},
  \code{\link{gpd.fit}}}

\examples{
data(rain)
pp.fit(rain, 10)
}
\keyword{models}
