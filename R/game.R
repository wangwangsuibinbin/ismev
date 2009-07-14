## Copyright (C) 2012 Marius Hofert and Valerie Chavez-Demoulin
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


## *G*eneralized *A*dditive *M*odel for *E*xtremes
##
## Remark:
## 1) Notation:
##    paper | Valerie's thesis | Valerie's former code | Coles' book + ismev
##      xi        -kappa                kappa                    xi
##    beta         sigma                sigma                 sigma
##      nu            nu                  eta                    --
## 2) GPD(xi in IR, beta > 0) distribution function (same as in ismev; up to notation):
##    G_{xi,beta}(x) = 1-(1+xi*x/beta)^(-1/xi) if xi!=0
##                   = 1-exp(-x/beta)          if xi =0
##    x>=0 when xi>=0 and x in [0,-beta/xi] when xi<0
## 3) GPD(xi, beta) density for xi>0:
##    g_{xi,beta}(x) = (1+xi*x/beta)^(-(1+1/xi))/beta if xi!=0, x>0


## loading required packages
require(mgcv) # comes with R natively (one of the recommended packages)
require(ismev) # for computing initial values xi.init, beta.init (not required to load within package ismev)


### auxiliary functions ########################################################

##' Reparameterized log-likelihood l^r(xi, nu; y_1,..,y_n) = l(xi, exp(nu)/(1+xi);
##' y_1,..,y_n), where l(xi, beta; y_1,..,y_n) = -n*log(beta)-(1+1/xi)*sum(
##' log(1+xi*y_i/beta)) denotes the log-likelihood of a GPD(xi, beta) distribution
##'
##' @title reparameterized log-likelihood l^r
##' @param y vector of exceedances (over a high threshold u)
##' @param xi GPD(xi,beta) (single) parameter xi
##' @param nu GPD(xi, beta) (single) parameter nu (orthogonal in the Fisher
##'        information metric to xi)
##' @param multivariate logical indicating whether xi and nu may be vectors of the
##'        same length as y
##' @return reparametrized log-likelihood l^r
##' @author Marius Hofert
rLogL <- function(y, xi, nu, multivariate=TRUE)
{
    stopifnot((n <- length(y)) > 0)
    if(multivariate){ # new setup
        stopifnot(length(xi)==n, length(nu)==n)
        sum(mapply(rLogL, y=y, xi=xi, nu=nu, multivariate=FALSE)) # note: this leads to n=length(y)==1 in the next call of rLogL (so n==1 and sum() goes over 1 obs. => correct)
    } else { # classical setup
        stopifnot(length(xi)==1, length(nu)==1) # classical setup
        -n*(nu-log1p(xi))-(1+1/xi)*sum(log1p(xi*(1+xi)*exp(-nu)*y))
    }
}

##' Compute derivatives of the reparameterized log-likelihood
##'
##' @title Compute derivatives of the reparameterized log-likelihood
##' @param y vector of exceedances (over a high threshold u)
##' @param xi vector of GPD(xi,beta) parameters xi
##' @param nu vector of (orthogonal in the Fisher information metric)
##'        GPD(xi, beta) parameters nu
##' @param verbose logical indicating whether wrong arguments to log1p are printed
##' @return (n x 4) matrix containing the partial derivatives of the
##'         reparameterized log-likelihood l^r where
##'         column 1: derivative of the reparameterized log-likelihood w.r.t. xi
##'         column 2: derivative of the reparameterized log-likelihood w.r.t. nu
##'         column 3: 2nd derivative of the reparameterized log-likelihood w.r.t. xi
##'         column 4: 2nd derivative of the reparameterized log-likelihood w.r.t. nu
##' @author Marius Hofert
##' Note: Column 3 and 4 have different sign than in the old code
DrLogL <- function(y, xi, nu, verbose=TRUE)
{
    stopifnot((n <- length(y)) > 0, length(xi)==n, length(nu)==n)

    ## auxiliary results
    enu <- exp(nu)
    xi1 <- 1+xi
    b <- enu/xi1 # beta
    bxiy <- b+xi*y # beta+xi*y
    ybxiy <- y/bxiy # y/(beta+xi*y)
    bbxiy <- b*bxiy # beta*(beta+xi*y)

    ## check (since this is sometimes a problematic case)
    if(verbose && any(ind <- ((value <- xi*y/b) <= -1))){
        warning("Some 1 + xi * y / b are <= 0, which, theoretically, should not happen:")
        tab <- cbind(xi=xi[ind], b=b[ind], y=y[ind], value=value[ind])
        nr <- nrow(tab)
        if (nr > 10) {
            print(tab[1:10,])
            cat("... ... ... ...\n")
        } else {
            print(tab)
        }
    }
    lxiyb <- log1p(xi*y/b) # log(1+xi*y/beta)

    ## components of the score function l(\bm{xi},\bm{beta};y_i, i=1,..,n')
    ## corresponding to l(xi(_i),beta(_i);y(_i))=log g_{xi(_i),beta(_i)}(y(_i)),
    ## => sum to get l(\bm{xi},\bm{beta};y_i, i=1,..,n')
    l.b <- (xi1*ybxiy-1) / b # derivative of rLogL w.r.t. beta
    l.xi <- 1/xi^2 * lxiyb - (1+1/xi) * ybxiy # derivative of rLogL w.r.t. xi

    ## components of the observed Fisher information
    l.bb <- -(y/b + (y-b)/bxiy) / bbxiy # 2nd derivative of rLogL w.r.t. beta
    l.xixi <- -2/xi^3 * lxiyb + ybxiy * (2/xi^2 + (1+1/xi) * ybxiy) # 2nd derivative of rLogL w.r.t. xi
    l.xib <- ybxiy * ((1+1/xi)/bxiy - 1/(xi*b)) # mixed partial derivative of rLogL w.r.t. xi and beta

    ## differentiate beta = beta(xi,nu) = exp(nu)/(1+xi) w.r.t. xi and nu
    b.xi <- -enu/xi1^2 # derivative of beta w.r.t. xi
    b.nu <- b # derivative of beta w.r.t. nu
    b.xixi <- 2*enu/xi1^3 # 2nd derivative of beta w.r.t. xi
    b.nunu <- b # 2nd derivative of beta w.r.t. nu

    ## Now let's derive the form of the reparameterized log-likelihood l^r (=rl)
    rl.xi <- l.xi + l.b * b.xi
    rl.xixi <- l.xixi + 2 * l.xib * b.xi + l.bb * b.xi^2 + l.b * b.xixi
    rl.nu <- l.b * b.nu
    rl.nunu <- l.bb * b.nu^2 + l.b * b.nunu

    ## return partial derivatives of the reparameterized log-likelihood l^r
    ## w.r.t. xi and nu
    cbind(rl.xi=rl.xi, rl.nu=rl.nu, rl.xixi=rl.xixi, rl.nunu=rl.nunu)
}

##' Function to check and adjust second-order derivatives
##'
##' @title Function to check and adjust second-order derivatives
##' @param der vector second-order derivatives
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives are printed
##' @return adjusted derivatives
##' @author Marius Hofert
##' Note: That's a helper function of gamGPDfitUp()
adjustD2 <- function(der, verbose=TRUE)
{
    ## setup
    isFinite <- is.finite(der)
    isNonFinite <- !isFinite
    isNonNeg <- isFinite & der >= 0
    ## adjust the derivatives
    isOK <- isFinite & !isNonNeg # logical indicating if der is okay (-Inf < . < 0)
    if(any(!isOK)){ # if not all derivativess are okay...
        if(sum(isOK)==0) stop("Can't adjust the derivatives, there are no finite, negative derivatives")
        der[!isOK] <- mean(der[isOK]) # Note: there is no justification for this (just a quick-and-dirty solution)
    }
    ## throw warning
    n <- length(der)
    if(verbose && any(isNonFinite)){
        numNA <- sum(isNonFinite)
        perNA <- round(100*numNA/n, 2)
        warning(numNA," (",perNA,"%) non-finite derivatives found and adjusted")
    }
    if(verbose && any(isNonNeg)){
        numNonNeg <- sum(isNonNeg)
        perNonNeg <- round(100*numNonNeg/n, 2)
        warning(numNonNeg," (",perNonNeg,"%) non-negative second-order derivatives found and adjusted")
    }
    ## return
    der
}

##' Compute update (one iteration) in gamGPDfit()
##'
##' @title Compute update (one iteration) in gamGPDfit()
##' @param y data.frame containing the exceedances over the threshold in a column
##'        labeled yname
##' @param xi.nu 2-column matrix of GPD parameters (xi,nu) to be updated where
##'        nu is orthogonal to xi in the Fisher information metric
##' @param xiFrhs right-hand side of the formula for xi in the gam() call
##'        for fitting xi
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param eyname string containing the name of the column of y which contains
##'        the exceedances
##' @param verbose logical indicating whether warnings about adjustments of
##'        the derivatives and wrong arguments in DrLogL() are printed
##' @param ... additional arguments passed to gam()
##' @return a list of length four containing
##'         element 1 (xi): object of class gamObject for xi as returned by mgcv::gam()
##'         element 2 (nu): object of class gamObject for nu as returned by mgcv::gam()
##'         element 3 (xi.weights): weights associated with xi
##'         element 4 (nu.weights): weights associated with nu
##' @author Marius Hofert
##' Note: That's a helper function of gamGPDfit()
gamGPDfitUp <- function(y, xi.nu, xiFrhs, nuFrhs, yname, verbose=TRUE, ...)
{
    stopifnot(is.data.frame(y), (dim. <- dim(y))[2] >= 1,
              length(which(colnames(y)==yname))==1,
              (n <- dim.[1]) > 0, dim(xi.nu)==c(n, 2),
              require(mgcv))

    ## pick out xi and nu (for readability)
    xi <- xi.nu[,1]
    nu <- xi.nu[,2]

    ## compute one Newton step in xi
    DrLL <- DrLogL(y[,yname], xi=xi, nu=nu, verbose=verbose) # (n1,4) matrix
    rl.xi <- DrLL[,"rl.xi"] # score in xi
    rl.xixi. <- adjustD2(DrLL[,"rl.xixi"], verbose=verbose) # -weight
    Newton.xi <- xi - rl.xi / rl.xixi. # Newton step

    ## concatenate Newton.xi and rl.xixi. to y, build formula, and estimate xi
    if("Newton.xi" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.xi'")
    y. <- cbind(y, Newton.xi=Newton.xi, rl.xixi.=rl.xixi.)
    xi.formula <- update(xiFrhs, Newton.xi~.) # build formula Newton.xi ~ xiFrhs
    xi.obj <- gam(xi.formula, data=y., weights=-rl.xixi., ...) # updated xi object of type gamObject

    ## build fitted (xi) object and check
    xi.fit <- fitted(xi.obj)
    if((n. <- length(xi.fit)) != n) stop(paste("length(xi.fit) = ",n.," != ",n,". This most likely comes from non-finite weights in the call to adjustD2()",sep=""))

    ## compute one Newton step in nu (for given new xi)
    DrLL <- DrLogL(y[,yname], xi=xi.fit, nu=nu, verbose=verbose) # (n1,4) matrix
    rl.nu <- DrLL[,"rl.nu"] # score in nu
    rl.nunu. <- adjustD2(DrLL[,"rl.nunu"], verbose=verbose) # -weight
    Newton.nu <- nu - rl.nu / rl.nunu. # Newton step

    ## concatenate Newton.nu and rl.nunu. to y, build formula, and estimate nu
    if("Newton.nu" %in% colnames(y))
        stop("y is not allowed to have a column named 'Newton.nu'")
    y. <- cbind(y, Newton.nu=Newton.nu, rl.nunu.=rl.nunu.)
    nu.formula <- update(nuFrhs, Newton.nu~.) # build formula Newton.nu ~ nuFrhs
    nu.obj <- gam(nu.formula, data=y., weights=-rl.nunu., ...) # updated nu object of type gamObject

    ## return list of two gamObject objects (for xi, nu)
    list(xi=xi.obj, nu=nu.obj, xi.weights=-rl.xixi., nu.weights=-rl.nunu.) # note: the naming (xi, nu) is not ideal but guarantees that colnames(param.old) = colnames(param.new) in gamGPDfit
}


### gamGPDfit ##################################################################

##' Semi-parametric estimation of GPD parameters via penalized maximum likelihood
##' estimation based on an iteratively reweighted least squares (backfitting)
##' algorithm
##'
##' @title Semi-parametric fitting of GPD parameters via penalized
##'        maximum likelihood estimation
##' @param x data.frame containing the losses (all other columns are treated
##'        as covariates)
##' @param threshold POT threshold above which losses are considered
##' @param nexceed number of exceedances
##' @param datvar name of the data column which contains the data to be modeled,
##'        for example, the losses
##' @param xiFrhs right-hand side of the formula for xi in the gam() call
##'        for fitting xi
##' @param nuFrhs right-hand side of the formula for nu in the gam() call
##'        for fitting nu
##' @param init bivariate vector containing initial values for (xi, beta)
##' @param niter maximal number of iterations in the backfitting algorithm
##' @param epsxi epsilon for stop criterion for xi
##' @param epsnu epsilon for stop criterion for nu
##' @param progress logical indicating whether progress information is displayed
##' @param verbose logical passed to gamGPDfitUp() (thus to DrLogL() and
##'        adjustD2())
##' @param ... additional arguments passed to gam() (called by gamGPDfitUp())
##' @return a list; see below
##' @author Marius Hofert
gamGPDfit <- function(x, threshold, nexceed=NULL, datvar,
                      xiFrhs, nuFrhs,
                      init=gpd.fit(x[, datvar], threshold=threshold, show=FALSE)$mle[2:1],
                      niter=128, epsxi=1e-5, epsnu=1e-5,
                      progress=TRUE, verbose=FALSE, ...)
{
    ## checks
    stopifnot(is.data.frame(x), length(init)==2, niter>=1, epsxi>0, epsnu>0,
              datvar %in% colnames(x))
    has.threshold <- !missing(threshold)
    has.nexceed <- !is.null(nexceed)
    if(has.threshold && has.nexceed) warning("Only one of 'threshold' and 'nexceed' is allowed -- will take 'threshold'") # both threshold and nexceed given
    if(!has.threshold && !has.nexceed) stop("Provide either 'threshold' or 'nexceed'") # none of threshold or nexceed given
    dim. <- dim(x) # dimension of x
    stopifnot((n <- dim.[1])>=1, dim.[2]>=2) # there should at least be one observation (actually, quite a bit more) and one column of covariates
    if(has.nexceed){ # nexceed given but no threshold
        stopifnot(0 < nexceed, nexceed <= n)
        threshold <- quantile(x[,datvar], probs=1-nexceed/n)
    } # => now we can work with threshold

    ## determine exceedances
    y. <- x[x[,datvar]>threshold,] # pick out exceedances; note: y. still contains covariates
    y.[,datvar] <- y.[,datvar] - threshold # replace exceedances by exceedance sizes
    n.ex <- nrow(y.) # number of exceedances

    ## initial values for xi, beta, and nu (nu = reparameterization of beta)
    xi.init <- init[1]
    beta.init <- init[2]
    nu.init <- log((1+xi.init) * beta.init)

    ## iteration ###############################################################

    iter <- 1 # iteration number
    updates <- list() # (empty) list of update (gamGPDfitUp) objects
    while(TRUE){

        ## update/fit parameters (xi, nu)
        param.old <- if(iter==1){
            matrix(rep(c(xi.init, nu.init), each=n.ex), ncol=2,
                       dimnames=list(rownames(y.), c("xi", "nu"))) # (n.ex,2)-matrix with cols "xi" and "nu"
        } else{
            param.new
        }
        updates[[iter]] <- gamGPDfitUp(y., xi.nu=param.old,
                                       xiFrhs=xiFrhs, nuFrhs=nuFrhs,
                                       yname=datvar, verbose=verbose, ...) # returns a list of two gam() objects containg the fitted (xi, nu)
        param.new <- sapply(updates[[iter]][c("xi", "nu")], fitted)
        ## note: param.old and param.new have the same rownames/colnames

        ## check
        if(any(dim(param.new)!=dim(param.old))) stop("gamGPDfitUp() returned an updated gamObject object of wrong dimension")
        if(progress){
            conv <- sapply(updates[[iter]][c("xi", "nu")], function(x) x$converged) # convergence status for (xi, nu)
            if(any(!conv)) warning("gam() in gamGPDfitUp() did not converge for ", paste(c("xi","nu")[!conv], collapse=", "))
        }

        ## tracing
        MRD.. <- colMeans(abs((param.old-param.new)/param.old)) # mean relative distance for (xi, nu)
        MRD. <- c(iter=iter, MRD..[1], MRD..[2])
        MRD <- if(iter==1) MRD. else rbind(MRD, MRD., deparse.level=0)
        if(progress){
            cat("Mean relative differences in iteration ", iter,
                " (xi, nu): (", sprintf("%g", MRD..[1]), ", ",
                sprintf("%g", MRD..[2]), ")\n", sep="")
        }

        ## check for "convergence"
        conv <- c(MRD..[1] <= epsxi, MRD..[2] <= epsnu)
        if(all(conv)) break
        if(iter >= niter){
            if(progress) warning("Reached 'niter' without the required precision for ",
                                 paste(c("xi","nu")[!conv], collapse=", "))
            break
        }
        iter <- iter + 1 # update iteration

    }

    ## finish and return #######################################################

    ## fitted parameters
    xi <- param.new[,"xi"] # fitted xi's
    nu <- param.new[,"nu"] # fitted nu's
    beta <- exp(nu)/(1+xi) # corresponding (fitted) beta

    ## log-likelihood
    logL <- rLogL(y=y.[,datvar], xi=xi, nu=nu) # reparameterized log-likelihood (reparameterization doesn't matter, it's the same *value* as the original log-likelihood)

    ## fitted gamObjects for xi and nu and standard errors
    ## (contained in updates but for convenience we treat this additionally)
    xiObj <- updates[[iter]]$xi
    nuObj <- updates[[iter]]$nu
    se.xi <- updates[[iter]]$xi.weights
    se.nu <- updates[[iter]]$nu.weights

    ## return list
    list(xi=xi, # estimated xi
         beta=beta, # estimated beta
         nu=nu, # estimated nu
         se.xi=se.xi, # standard error for xi
         se.nu=se.nu, # standard error for nu
         covar=y.[,-which(colnames(y.)==datvar), drop=FALSE], # covariates (corresponding to xi, beta, etc.)
         y=y.[,datvar], # exceedances
         res=log1p(y.[,datvar]*xi/beta)/xi, # residuals
         MRD=MRD, # mean relative distances between old/new (xi, nu) for all iterations
         logL=logL, # log-likelihood at the estimated parameters
         xiObj=xiObj, # gamObject for estimated xi (return object of mgcv::gam())
         nuObj=nuObj, # gamObject for estimated nu (return object of mgcv::gam())
         xiUpdates=lapply(updates, `[[`, "xi"), # updates for xi for each iteration (list of gamObject objects); contains xiObj as last element
         nuUpdates=lapply(updates, `[[`, "nu")) # updates for nu for each iteration (list of gamObject objects); contains nuObj as last element
}


### gamGPDboot #################################################################

##' Adjusted sampling for 1-element vectors
##'
##' @title Adjusted sampling for 1-element vectors
##' @param x see sample()
##' @param size see sample()
##' @param replace see sample()
##' @param prob see sample()
##' @return in case of length(x)==1, return x itself (rather than a sample from 1:x)
##'         otherwise return sample(x, ...)
##' @author Marius Hofert
sample.real <- function(x, size, replace=FALSE, prob=NULL)
    if(length(x)==1) x else sample(x, size=size, replace=replace, prob=prob)

##' Post-blackend bootstrap of Chavez-Demoulin and Davison (2005) for gamGPDfit()
##'
##' @title Post-blackend bootstrap of Chavez-Demoulin and Davison (2005) for
##'        gamGPDfit()
##' @param x see gamGPDfit()
##' @param B number of bootstrap replications
##' @param threshold see gamGPDfit()
##' @param nexceed see gamGPDfit()
##' @param datvar see gamGPDfit()
##' @param xiFrhs see gamGPDfit()
##' @param nuFrhs see gamGPDfit()
##' @param init see gamGPDfit()
##' @param niter see gamGPDfit()
##' @param epsxi see gamGPDfit()
##' @param epsnu see gamGPDfit()
##' @param boot.progress logical indicating whether progress information is displayed
##' @param progress see gamGPDfit() (only used if progress==TRUE)
##' @param verbose see gamGPDfit() (only used if progress==TRUE)
##' @param debug logical indicating whether initial fit is saved
##' @param ... see gamGPDfit()
##' @return a list; see below
##' @author Marius Hofert
gamGPDboot <- function(x, B, threshold, nexceed=NULL, datvar, xiFrhs, nuFrhs,
                       init=gpd.fit(x[, datvar], threshold=threshold, show=FALSE)$mle[2:1],
                       niter=128, epsxi=1e-5, epsnu=1e-5,
                       boot.progress=TRUE, progress=FALSE, verbose=FALSE,
                       debug=FALSE, ...)
{
    ## progress
    if(boot.progress) {
         if(progress) cat("\nStarting initial fit:\n") else {
             pb <- txtProgressBar(max=B+1, style=if(isatty(stdout())) 3 else 1) # setup progress bar
             on.exit(close(pb)) # close progress bar
        }
    }

    ## (major) fit using gamGPDfit()
    fit. <- gamGPDfit(x=x, threshold=threshold, nexceed=nexceed, datvar=datvar,
                      xiFrhs=xiFrhs, nuFrhs=nuFrhs, init=init, niter=niter,
                      epsxi=epsxi, epsnu=epsnu,
                      progress=if(!boot.progress) FALSE else progress,
                      verbose=if(!boot.progress) FALSE else verbose, ...)

    ## progress
    if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, 1)

    ## pick out fitted values
    xi <- fit.$xi # fitted xi
    nu <- fit.$nu # fitted nu
    beta <- exp(nu)/(1+xi) # fitted beta

    ## for debugging
    if(debug) save(fit., file="gamGPDboot_debug.rda")

    ## post-blackened bootstrap; see Chavez-Demoulin and Davison (2005)
    rnum <- 1 # iteration number
    bfit <- lapply(1:B, function(b){
        ## resample residuals within each group of (same) covariates,
        rr <- ave(fit.$res, fit.$covar,
                  FUN=function(r) sample.real(r, size=length(r), replace=TRUE))

        ## compute and put in corresponding exceedances (see Chavez-Demoulin and Davison (2005))
        ## Note: 1) they use a different reparameterization
        ##       2) solve rr=log(1+xi*y/beta)/xi w.r.t. y
        y. <- expm1(rr*xi)*beta/xi # reconstruct exceedances from residuals
        x. <- cbind(fit.$covar, y=y.) # add exceedances

        ## progress
        if(boot.progress && progress) cat("Starting fit in bootstrap run ",
                                          b, " of ", B, ":\n", sep="")

        ## call gamGPDfit()
        ## note: threshold=0 => we discard those exceedances which are equal to 0
        bfit. <- gamGPDfit(x=x., threshold=0, nexceed=NULL, datvar="y",
                           xiFrhs=xiFrhs, nuFrhs=nuFrhs,
                           init=gpd.fit(y., threshold=0, show=FALSE)$mle[2:1],
                           niter=niter, epsxi=epsxi, epsnu=epsnu,
                           progress=if(!boot.progress) FALSE else progress,
                           verbose=if(!boot.progress) FALSE else verbose, ...)

        ## progress
        if(boot.progress) if(progress) cat("\n") else setTxtProgressBar(pb, b+1)

        ## return fit
        bfit.
    })

    ## return list
    list(xi=cbind(fit.$xi, sapply(bfit, function(x) x$xi)), # estimated xi
         beta=cbind(fit.$beta, sapply(bfit, function(x) x$beta)), # estimated beta
         nu=cbind(fit.$nu, sapply(bfit, function(x) x$nu)), # estimated nu
         se.xi=cbind(fit.$se.xi, sapply(bfit, function(x) x$se.xi)), # standard error for xi
         se.nu=cbind(fit.$se.nu, sapply(bfit, function(x) x$se.nu)), # standard error for nu
         covar=fit.$covar, # covariates (corresponding to xi, beta, etc.; are the same for all bootstrap replications)
         y=cbind(fit.$y, sapply(bfit, function(x) x$y)), # exceedances
         res=cbind(fit.$res, sapply(bfit, function(x) x$res)), # residuals
         MRD=c(list(fit.$MRD), lapply(bfit, function(x) x$MRD)), # mean relative distances between old/new (xi, nu) for all iterations
         logL=c(fit.$logL, sapply(bfit, function(x) x$logL)), # log-likelihood at the estimated parameters
         xiObj=c(list(fit.$xiObj), lapply(bfit, function(x) x$xiObj)), # gamObject for estimated xi (return object of mgcv::gam())
         nuObj=c(list(fit.$nuObj), lapply(bfit, function(x) x$nuObj)), # gamObject for estimated nu (return object of mgcv::gam())
         xiUpdates=c(list(fit.$xiUpdates), lapply(bfit, function(x) x$xiUpdates)), # updates for xi for each iteration (list of gamObject objects); contains xiObj as last element
         nuUpdates=c(list(fit.$nuUpdates), lapply(bfit, function(x) x$nuUpdates))) # updates for nu for each iteration (list of gamObject objects); contains nuObj as last element
}


### functions to extract results ###############################################

##' Predict Lambda and compute (pointwise) 1-alpha confidence intervals
##'
##' @title Predict Lambda and compute (pointwise) 1-alpha confidence intervals
##' @param x object as returned by gam()
##' @param newdata 'newdata' object as required by predict(); named data.frame
##'        of type expand.grid(covar1=1:5, covar2=1:2)
##' @param alpha significance level
##' @return a (n, d+3) data.frame where n=nrow(newdata), d=ncol(newdata), and
##'         the last three columns contain the predicted Lambda and corresponding
##'         CIs
##' @author Marius Hofert
LambdaPredict <- function(x, newdata, alpha=0.05)
{
    LambdaPred <- predict(x, newdata=newdata, se.fit=TRUE) # predict object
    LambdaPr <- as.numeric(LambdaPred$fit) # predict Lambda
    LambdaSE <- as.numeric(LambdaPred$se.fit) # standard error
    qa2 <- qnorm(1-alpha/2) # 1-alpha/2 quantile of N(0,1)
    data.frame(newdata, # covariates as specified by newdata
               Lambda=exp(LambdaPr), # predicted Lambda
               CIlow=exp(LambdaPr-qa2*LambdaSE), # lower CI
               CIup=exp(LambdaPr+qa2*LambdaSE)) # upper CI
}

##' Predict xi and beta
##'
##' @title Predict xi and beta
##' @param x object as returned by gamGPDfit()
##' @param newdata 'newdata' object as required by predict(); named data.frame
##'        of type expand.grid(covar1=1:5, covar2=1:2)
##' @return a (n, d+2) data.frame where n=nrow(newdata), d=ncol(newdata), and
##'         the last two columns contain the corresponding predictions for
##'         xi and beta
##' @author Marius Hofert
##' Note: Standard errors would be available via predict(..., se.fit=TRUE), but
##'       only for xi (and nu), not for beta.
xibetaPredict <- function(x, newdata)
{
    xiPred <- as.numeric(predict(x$xiObj, newdata=newdata)) # predict xi
    nuPred <- as.numeric(predict(x$nuObj, newdata=newdata)) # predict nu
    data.frame(newdata, # covariates as specified by newdata
               xi=xiPred, # predicted xi
               beta=exp(nuPred)/(1+xiPred)) # predicted beta
}

##' Extract fits of xi and beta, and compute bootstrapped CIs from a gamGPDboot() object
##'
##' @title Extract fits of xi and beta, and compute bootstrapped CIs from a gamGPDboot() object
##' @param x object as returned by gamGPDboot()
##' @param alpha significance level
##' @return a (n, d+6) data.frame where n is the number of covariate combinations for
##'         which there are exceedances (x$covar), d the dimension of x$covar, and
##'         the columns contain the lower CI for xi, the fitted xi, the upper CI for xi,
##'         and the same for beta
##' @author Marius Hofert
xibetaFitCI <- function(x, alpha=0.05)
{
    ## define result object
    nr <- nrow(x$xi) # = nrow(x$beta) etc.
    res <- array(NA, dim=c(nr, 6),
                 dimnames=list(1:nr,
                               c("xi", "xiCIlow", "xiCIup", "beta", "betaCIlow", "betaCIup"))) # fitted xi's, beta's, and bootstrapped CIs

    ## fit
    xi <- x$xi # (nr, B+1) matrix; B = number of bootstrap replications
    beta <- x$beta # (nr, B+1) matrix; B = number of bootstrap replications
    res[,"xi"] <- xi[,1] # fitted xi
    res[,"beta"] <- beta[,1] # fitted beta

    ## CIs (derived from fitted values + bootstrapped ones)
    ## note: this can take some seconds
    qlow <- function(x, alpha) quantile(x, probs=alpha/2)[[1]] # lower CI
    qup <- function(x, alpha) quantile(x, probs=1-alpha/2)[[1]] # upper CI
    res[,"xiCIlow"] <- apply(xi, 1, qlow, alpha=alpha)
    res[,"xiCIup"] <- apply(xi, 1, qup, alpha=alpha)
    res[,"betaCIlow"] <- apply(beta, 1, qlow, alpha=alpha)
    res[,"betaCIup"] <- apply(beta, 1, qup, alpha=alpha)

    ## return
    data.frame(x$covar, res)
}

