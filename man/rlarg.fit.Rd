\name{rlarg.fit}
\alias{rlarg.fit}

\title{Maximum-likelihood Fitting of Order Statistics Model}
\description{
  Maximum-likelihood fitting for the order statistic model,
  including generalized linear modelling of each parameter.
}
\usage{
rlarg.fit(xdat, r = dim(xdat)[2], ydat = NULL, mul = NULL, sigl = NULL,
  shl = NULL, mulink = identity, siglink = identity, shlink = identity,
  show = TRUE, method = "Nelder-Mead", maxit = 10000, \dots)
}
\arguments{
  \item{xdat}{A numeric matrix of data to be fitted. Each row
    should be a vector of decreasing order, containing the
    largest order statistics for each year (or time period).
    The first column therefore contains annual (or period)
    maxima.
    Only the first \code{r} columns are used for the fitted
    model. By default, all columns are used.
    If one year (or time period) contains fewer order statistics
    than another, missing values can be appended to the end
    of the corresponding row.}
  \item{r}{The largest \code{r} order statistics are used for
    the fitted model.}
  \item{ydat}{A matrix of covariates for generalized linear modelling
    of the parameters (or \code{NULL} (the default) for stationary
    fitting). The number of rows should be the same as the number
    of rows of \code{xdat}.}
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
  approximately centered and scaled).
}
\value{
  A list containing the following components. A subset of these
  components are printed after the fit. If \code{show} is
  \code{TRUE}, then assuming that successful convergence is
  indicated, the components \code{nllh}, \code{mle} and \code{se}
  are always printed.
  
  \item{trans}{An logical indicator for a non-stationary fit.}
  \item{model}{A list with components \code{mul}, \code{sigl}
    and \code{shl}.}
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
  \item{vals}{A matrix with three columns containing the maximum
    likelihood estimates of the location, scale and shape parameters
    at each data point.}
  \item{r}{The number of order statistics used.}
}

\seealso{\code{\link{rlarg.diag}}, \code{\link{optim}}}

\examples{
\dontrun{data(venice)}
\dontrun{rlarg.fit(venice[,-1])}
}
\keyword{models}
