\name{gpd.fit}
\alias{gpd.fit}

\title{Maximum-likelihood Fitting for the GPD Model}
\description{
  Maximum-likelihood fitting for the GPD model,
  including generalized linear modelling of each parameter.
}
\usage{
gpd.fit(xdat, threshold, npy = 365, ydat = NULL, sigl = NULL,
    shl = NULL, siglink = identity, shlink = identity, siginit = NULL,
    shinit = NULL, show = TRUE,
    method = "Nelder-Mead", maxit = 10000, \dots)
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
  \item{sigl, shl}{Numeric vectors of integers, giving the columns
    of \code{ydat} that contain covariates for generalized linear
    modelling of the scale and shape parameters repectively
    (or \code{NULL} (the default) if the corresponding parameter is
    stationary).}
  \item{siglink, shlink}{Inverse link functions for generalized
    linear modelling of the scale and shape parameters repectively.}
  \item{siginit, shinit}{numeric giving initial value(s) for parameter
	estimates.  If NULL, default is sqrt(6 * var(xdat))/pi and 0.1
	for the scale and shape parameters, resp.  If using parameter
	covariates, then these values are used for the constant term,
	and zeros for all other terms.}
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
  approximately centered and scaled).
}
\value{
  A list containing the following components. A subset of these
  components are printed after the fit. If \code{show} is
  \code{TRUE}, then assuming that successful convergence is
  indicated, the components \code{nexc}, \code{nllh},
  \code{mle}, \code{rate} and \code{se} are always printed.
  
  \item{trans}{An logical indicator for a non-stationary fit.}
  \item{model}{A list with components \code{sigl} and \code{shl}.}
  \item{link}{A character vector giving inverse link functions.}
  \item{threshold}{The threshold, or vector of thresholds.}
  \item{nexc}{The number of data points above the threshold.}
  \item{data}{The data that lie above the threshold. For
    non-stationary models, the data is standardized.}
  \item{conv}{The convergence code, taken from the list returned by
    \code{\link{optim}}. A zero indicates successful convergence.}
  \item{nllh}{The negative logarithm of the likelihood evaluated at
    the maximum likelihood estimates.}
  \item{vals}{A matrix with three columns containing the maximum
    likelihood estimates of the scale and shape parameters, and
    the threshold, at each data point.} 
  \item{mle}{A vector containing the maximum likelihood estimates.}
  \item{rate}{The proportion of data points that lie above the
    threshold.}
  \item{cov}{The covariance matrix.}
  \item{se}{A vector containing the standard errors.}
  \item{n}{The number of data points (i.e.\ the length of
    \code{xdat}).}
  \item{npy}{The number of observations per year/block.}
  \item{xdata}{The data that has been fitted.} 
}

\seealso{\code{\link{gpd.diag}}, \code{\link{optim}},
  \code{\link{gpd.prof}}, \code{\link{gpd.fitrange}},
  \code{\link{mrl.plot}}, \code{\link{pp.fit}}}

\examples{
data(rain)
gpd.fit(rain, 10)
}
\keyword{models}
