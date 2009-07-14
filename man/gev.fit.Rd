\name{gev.fit}
\alias{gev.fit}

\title{Maximum-likelihood Fitting of the GEV Distribution}
\description{
  Maximum-likelihood fitting for the generalized extreme value distribution,
  including generalized linear modelling of each parameter.
}
\usage{
gev.fit(xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, 
    mulink = identity, siglink = identity, shlink = identity, 
    muinit = NULL, siginit = NULL, shinit = NULL,
    show = TRUE, method = "Nelder-Mead", maxit = 10000, \dots)
}
\arguments{
  \item{xdat}{A numeric vector of data to be fitted.}
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
  \item{muinit, siginit, shinit}{numeric of length equal to total number
    of parameters used to model the location, scale or shape parameter(s),
    resp.  See Details section for default (NULL) initial values.}
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
  The form of the GEV used is that of Coles (2001) Eq (3.2).  Specifically,
  positive values of the shape parameter imply a heavy tail, and negative values
  imply a bounded upper tail.

  For non-stationary fitting it is recommended that the covariates
  within the generalized linear models are (at least approximately)
  centered and scaled (i.e.\ the columns of \code{ydat} should be
  approximately centered and scaled).

  Let m=mean(xdat) and s=sqrt(6*var(xdat))/pi.  Then, initial values
  assigend when 'muinit' is NULL are m - 0.57722 * s (stationary case).
  When 'siginit' is NULL, the initial value is taken to be s, and when
  'shinit' is NULL, the initial value is taken to be 0.1.  When
  covariates are introduced (non-stationary case), these same initial
  values are used by default for the constant term, and zeros for all
  other terms.  For example, if a GEV( mu(t)=mu0+mu1*t, sigma, xi) is
  being fitted, then the initial value for mu0 is m - 0.57722 * s, and
  0 for mu1.
}
\value{
  A list containing the following components. A subset of these
  components are printed after the fit. If \code{show} is
  \code{TRUE}, then assuming that successful convergence is
  indicated, the components \code{nllh}, \code{mle} and \code{se}
  are always printed.

  \item{nllh}{single numeric giving the negative log-likelihood value.}  
  \item{mle}{numeric vector giving the MLE's for the location, scale and shape parameters, resp.}
  \item{se}{numeric vector giving the standard errors for the MLE's for the location, scale and shape parameters, resp.}
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
}

\references{
Coles, S., 2001.  An Introduction to Statistical Modeling of Extreme Values.  Springer-Verlag, London, U.K., 208pp.
}

\seealso{\code{\link{gev.diag}}, \code{\link{optim}},
  \code{\link{gev.prof}}}

\examples{
data(portpirie)
gev.fit(portpirie[,2])
}
\keyword{models}
