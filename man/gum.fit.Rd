\name{gum.fit}
\alias{gum.fit}

\title{Maximum-likelihood Fitting of the Gumbel Distribution}
\description{
  Maximum-likelihood fitting for the gumbel distribution,
  including generalized linear modelling of each parameter.
}
\usage{
gum.fit(xdat, ydat = NULL, mul = NULL, sigl = NULL, mulink = identity,
    siglink = identity, show = TRUE, method = "Nelder-Mead",
    maxit = 10000, \dots)
}
\arguments{
  \item{xdat}{A numeric vector of data to be fitted.}
  \item{ydat}{A matrix of covariates for generalized linear modelling
    of the parameters (or \code{NULL} (the default) for stationary
    fitting). The number of rows should be the same as the length
    of \code{xdat}.}
  \item{mul, sigl}{Numeric vectors of integers, giving the columns
    of \code{ydat} that contain covariates for generalized linear
    modelling of the location and scale parameters repectively
    (or \code{NULL} (the default) if the corresponding parameter is
    stationary).}
  \item{mulink, siglink}{Inverse link functions for generalized
    linear modelling of the location and scale parameters
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
  approximately centered and scaled).
}
\value{
  A list containing the following components. A subset of these
  components are printed after the fit. If \code{show} is
  \code{TRUE}, then assuming that successful convergence is
  indicated, the components \code{nllh}, \code{mle} and \code{se}
  are always printed.
  
  \item{trans}{An logical indicator for a non-stationary fit.}
  \item{model}{A list with components \code{mul} and \code{sigl}.}
  \item{link}{A character vector giving inverse link functions.}
  \item{conv}{The convergence code, taken from the list returned by
    \code{\link{optim}}. A zero indicates successful convergence.}
  \item{nllh}{The negative logarithm of the likelihood evaluated at
    the maximum likelihood estimates.}
  \item{data}{The data that has been fitted. For non-stationary
    models, the data is standardized.}
  \item{mle}{A vector containing the maximum likelihood estimates.}
  \item{cov}{The covariance matrix.}
  \item{se}{A vector containing the standard errors.}
  \item{vals}{A matrix with two columns containing the maximum
    likelihood estimates of the location and scale parameters
    at each data point.}  
}

\seealso{\code{\link{gum.diag}}, \code{\link{optim}},
  \code{\link{gev.fit}}}

\examples{
data(portpirie)
gum.fit(portpirie[,2])
}
\keyword{models}
