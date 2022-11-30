#' Nonparametric test of flatness of a causal dose-response curve under exogeneity
#'
#' This function performs a hypothesis test that the causal dose-response curve theta(a) is flat on the support of the observed exposure A. The exposure may be discrete, continuous, or an arbitrary mixture of discrete and continuous components.  See the accompanying paper for details.
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of exposure values.
#' @param W \code{n x p} data.frame of potential confounders.
#' @param p vector of exponents to use in the norm of Omega. Defaults to 2.
#' @param control Optional list of control parameters. See \code{\link{causalNullTest.control}} for details.
#' @return \code{causalNullTest} returns a named list with the following elements:
#' \item{test}{A data.frame containing columns p (the exponent in the Lp norm used in the test), obs.stat (the observed Lp norm of the primitive Omega), p.val (the p-value corresponding to the test), ci.ll (the lower confidence limit corresponding to the Lp norm of the primitive), and ci.ul (the upper confidence limit corresponding to the Lp norm of the primitive).}
#' \item{control}{Controls used in fitting.}
#' If If \code{control$return.Omega == TRUE}, the following elements are also included in the output:
#' list(Omega.hat = data.frame(a=a.vals, Omega.hat), IF.vals = IF.vals, paths = paths))
#' \item{Omega.hat}{A data.frame with values of a and the corresponding estimated function Omega.hat(a).}
#' \item{IF.vals}{\code{n x m} matrix of influence function values of Omega.hat at every unique sorted value of the exposure. Each column is the influence function for a separate value of a, the values of which can be found in \code{Omega.hat}.}
#' \item{paths}{\code{control$n.sim x m} matrix of simulated paths from the limiting Gaussian process. Each row is an independent simulated path along the values of a.}
#'  If \code{save.nuis.fits = TRUE}, then the following elements are also included in the output:
#'  \item{mu.hat}{The estimated outcome regression function, as a list of fits if \code{control$cross.fit = TRUE}.}
#'  \item{g.hat}{The estimated propensity, possibly as a list of fits if \code{control$cross.fit = TRUE}.}
#'
#' @examples
#' # Sample data
#' n <- 1000
#' W <- data.frame(W1 = runif(n))
#' Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W$W1)))
#' A <- (1-Z) * rnorm(n, mean = W$W1, sd = abs(1 + W$W1))
#' Y <- rexp(n, rate = 1+abs(W$W1 * A))
#' causalNullTest(Y, A, W, p = c(1,2,Inf), control = list(cross.fit = FALSE, verbose=TRUE, g.n.bins = 2:5))

causalNullTest <- function(Y, A, W, p=2, control = list()) {

  call <- match.call(expand.dots = TRUE)
  control <- do.call("causalNullTest.control", control)

  .check.input(Y=Y, A=A, W=W, p=p, control=control)

  library(mvtnorm)
  n <- length(Y)
  a.vals <- sort(unique(A))

  if(control$cross.fit & is.null(control$folds)) {
    control$folds <- sample(rep(1:control$V, length.out = n), replace=FALSE)
  }

  if(is.null(control$mu.hat)) {
    library(SuperLearner)
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating outcome regressions...")
      control$mu.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v, "...")
        if(length(setdiff(Y, c(0,1))) == 0) {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
        } else {
          fit <- SuperLearner(Y = Y[control$folds != v], X = cbind(A, W)[control$folds != v,], SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
        }
        function(a, w) c(predict(fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating outcome regression...")
      if(length(setdiff(Y, c(0,1))) == 0) {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library, family = 'binomial', method = 'method.NNloglik')
      } else {
        mu.fit <- SuperLearner(Y = Y, X = cbind(A, W), SL.library = control$mu.SL.library, family = 'gaussian', method = 'method.NNLS')
      }
      control$mu.hat <- function(a, w) c(predict(mu.fit, newdata = cbind(A = a, w),onlySL = TRUE)$pred)
      if(control$verbose) cat("\n")
    }
  }

  if(is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(control$verbose) cat("Estimating propensities..")
      control$g.hat <- lapply(1:control$V, function(v) {
        if(control$verbose) cat("fold", v)
        fit <- cmdSuperLearner(A = A[control$folds != v], W = W[control$folds != v,,drop=FALSE], newA = A[control$folds == v], newW = W[control$folds == v,,drop=FALSE], control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose, saveFitLibrary = FALSE))
        c(fit$SL.densities)
        #function(a, w) c(predict.cmdSuperLearner(fit, newA = a, newW = w))
      })
      if(control$verbose) cat("\n")
    } else {
      if(control$verbose) cat("Estimating propensity...")
      g.fit <- cmdSuperLearner(A = A, W = W, control=list(SL.library = control$g.SL.library, n.bins = control$g.n.bins, verbose = control$verbose, saveFitLibrary = FALSE))
      control$g.hat <- g.fit$SL.densities
      rm(g.fit)
      #control$g.hat <- function(a, w) c(predict.cmdSuperLearner(g.fit, newA = a, newW = w))

      if(control$verbose) cat("\n")
    }
  }

  if(control$verbose) cat("Computing Omega...")
  if(!control$cross.fit) {
    ord <- order(A)
    A <- A[ord]
    Y <- Y[ord]
    W <- W[ord,,drop=FALSE]
    if(inherits(control$g.hat, "function")) {
      g.hats <- control$g.hat(A, W)
      control$g.hat <- NULL
    }
    else g.hats <- control$g.hat
    if(any(g.hats < control$g.trunc)) {
      warning("Truncating g.hats below. Possible positivity issues.")
      g.hats[g.hats < control$g.trunc] <- control$g.trunc
    }
    a.ecdf <- ecdf(A)
    a.weights <- sapply(a.vals, function(a0) mean(A == a0))
    A.a.val <- sapply(A, function(a0) which(a.vals == a0))
    u.vals <- a.ecdf(a.vals)
    mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat(a0, W)) #rows index W, columns index a.vals
    control$mu.hat <- NULL
    mu.hats <- mu.hats.a.vals[,A.a.val]
    theta.a.vals <- colMeans(mu.hats.a.vals)
    theta.A <- theta.a.vals[A.a.val]
    mu.hats.data <- diag(mu.hats)
    partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / n
    gamma.hat <- mean(mu.hats)
    Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A <= a0) * theta.A)) - gamma.hat * u.vals

    IF.vals <- sapply(a.vals, function(a0) {
      if(any(A <= a0))  mumean.vals <- partial.mu.means[,max(which(A <= a0))]
      else mumean.vals <- 0
      (as.numeric(A <= a0) - a.ecdf(a0)) * ((Y - mu.hats.data) / g.hats + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,n] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
    })

    Omega.hat <- colMeans(IF.vals) + Omega.a.vals

    if(control$verbose) cat("\nComputing covariance...\n")

    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(IF.vals[,s] * IF.vals[,t])
    }))
  }
  else {

    fold.Omega.hats <- matrix(NA, nrow = control$V, ncol = length(a.vals))
    IF.vals <- vector(length=control$V, mode='list')
    for(j in 1:control$V) {
      if(control$verbose) cat("fold", j, "...")
      Nv <- sum(control$folds == j)
      A.test <- A[control$folds == j]
      Y.test <- Y[control$folds == j]
      W.test <- W[control$folds == j,, drop=FALSE]
      ord <- order(A.test)
      A.test <- A.test[ord]
      Y.test <- Y.test[ord]
      W.test <- W.test[ord,, drop=FALSE]
      if(inherits(control$g.hat[[j]], "function")) {
        g.hats.test <- control$g.hat[[j]](a = A.test, w = W.test)
      }
      else g.hats.test <- control$g.hat[[j]]

      if(any(g.hats.test < control$g.trunc)) {
        warning("Truncating g.hats below. Possible positivity issues.")
        g.hats.test[g.hats.test < control$g.trunc] <- control$g.trunc
      }
      a.ecdf <- ecdf(A.test)
      a.weights <- sapply(a.vals, function(a0) mean(A.test == a0))
      A.a.val <- sapply(A.test, function(a0) which(a.vals == a0))
      u.vals <- a.ecdf(a.vals)
      mu.hats.a.vals <- sapply(a.vals, function(a0) control$mu.hat[[j]](a=a0, w=W.test)) #rows index W, columns index a.vals
      mu.hats <- mu.hats.a.vals[,A.a.val]
      theta.a.vals <- colMeans(mu.hats.a.vals)
      theta.A <- theta.a.vals[A.a.val]
      mu.hats.data <- diag(mu.hats)
      partial.mu.means <- t(apply(mu.hats, 1, cumsum)) / Nv
      gamma.hat <- mean(mu.hats)
      Omega.a.vals <- sapply(a.vals, function(a0) mean(as.numeric(A.test <= a0) * theta.A)) - gamma.hat * u.vals

      IF.vals[[j]] <- sapply(a.vals, function(a0) {
        if(any(A.test <= a0)) mumean.vals <- partial.mu.means[,max(which(A.test <= a0))]
        else mumean.vals <- 0
        (as.numeric(A.test <= a0) - a.ecdf(a0)) * ((Y.test - mu.hats.data) / g.hats.test + theta.A - gamma.hat) + mumean.vals - partial.mu.means[,ncol(partial.mu.means)] * a.ecdf(a0) - 2 * Omega.a.vals[which(a.vals == a0)]
      })

      fold.Omega.hats[j,] <- colMeans(IF.vals[[j]]) + Omega.a.vals
    }

    Omega.hat <- colMeans(fold.Omega.hats)
    if(control$verbose) cat("\nComputing covariance...\n")
    Sigma.hat <- sapply(1:length(a.vals), function(s) sapply(1:length(a.vals), function(t) {
      mean(unlist(lapply(IF.vals, function(IF) mean(IF[,s] * IF[,t]))))
    }))
  }

  if(control$verbose) cat("Simulating paths...\n")

  paths <- rmvnorm(control$n.sim, sigma=Sigma.hat)

  if(control$verbose) cat("Computing statistics...\n")

  a.weights <- sapply(a.vals, function(a) mean(A == a))
  ret <- t(sapply(p, function(pp) {
    stat <- ifelse(pp < Inf, (sum(abs(Omega.hat )^pp * a.weights))^{1/pp}, max(abs(Omega.hat)))

    if(pp < Inf) {
      stats <- (apply(abs(paths)^pp, 1, function(row) sum(row * a.weights)))^{1/pp}
    } else {
      stats <- apply(abs(paths), 1, max)
    }

    p.val <- mean(stats / sqrt(n) > stat)

    q <- quantile(stats, (1 - (1-control$conf.level) / 2))
    ci.ll <- max(stat - q / sqrt(n), 0)
    ci.ul <- stat + q / sqrt(n)

    res <- c(stat, p.val, ci.ll, ci.ul)

    res

  }))
  ret.df <- data.frame(p = p, obs.stat = ret[,1], p.val = ret[,2], ci.ll = ret[,3], ci.ul = ret[,4])
  ret.list <- list(test = ret.df)
  if(control$return.Omega) {
    ret.list <- c(ret.list, list(Omega.hat = data.frame(a=a.vals, Omega.hat), IF.vals = IF.vals, paths = paths))
  }
  if(control$save.nuis.fits) {
    ret.list <- c(ret.list, mu.hat = control$mu.hat, g.hat = control$g.hat)
    if(control$cross.fit) ret.list <- c(ret.list, folds = control$folds)
  }

  return(ret.list)
}


#' Initialize control parameters for causalNullTest
#'
#' This function initializes the control parameters for use in \code{causalNullTest}. The outcome regression function mu is by default estimated using \code{\link[SuperLearner]{SuperLearner}}, and the propensity is estimated using the conditional mixed density method implemented in \code{\link{cmdSuperLearner}}. Alternatively, the estimation process can be overriden by providing predictions from pre-fit nuisance estimators.
#'
#' @param mu.SL.library Library of candidate learners for the outcome regression to be passed on to \code{\link[SuperLearner]{SuperLearner}}. Ignored if \code{mu.hats} is provided.
#' @param g.SL.library Library of candidate learners for the outcome regression to be passed on to \code{\link{cmdSuperLearner}}.
#' @param g.n.bins Numeric vector of number of bins to use for estimation of the propensity. Passed on to \code{\link{cmdSuperLearner}}.
#' @param cross.fit Logical indicating whether to cross-fit nuisance parameters. Defaults to \code{TRUE}.
#' @param V Positive integer number of folds for cross-fitting. Defaults to 10.
#' @param folds Optional \code{n x 1} numeric vector indicating which fold each observation is in.
#' @param save.nuis.fits Logical indicating whether to save the fitted nuisance objects.
#' @param mu.hat Optional pre-fit outcome regression. If \code{cross.fit} is \code{FALSE}, then a function that takes arguments \code{a} (a vector) and \code{w} (a data.frame) and returns predictions of the outcome regression function. If \code{cross.fit} is \code{TRUE}, then a list of functions of length \code{V} with the fitted outcome regression functions on each of the \code{V} training sets. If provided as a list of functions, then \code{folds} must be provided. If \code{mu.SL.library = NULL}, then \code{mu.hats} must be specified.
#' @param g.hat Optional pre-fit treatment propensities. If \code{cross.fit} is \code{FALSE}, then a function that takes arguments \code{a} (a vector) and \code{w} (a data.frame) and returns predictions of the standardized propensity function OR an \code{n x 1} vector of fitted propensities. If \code{cross.fit} is \code{TRUE}, then a list of functions of length \code{V} with the fitted standardized propensity functions on each of the \code{V} training sets OR a list of length \code{V} of fitted out-of-sample propensities on each of the folds. If provided as a list, then \code{folds} must be provided. If \code{g.SL.library = NULL}, then \code{g.hats} must be specified.
#' @param g.trunc Value at which to truncate predicted propensities from below. Any propensity values less than \code{g.trunc} will be set to \code{g.trunc}.
#' @param n.sim Number of simulations to use for the limiting Gaussian process in computing approximate quantiles.
#' @param return.Omega Logical indicating whether to return the estimated primitive function Omega.
#' @param conf.level Optional confidence level to use for computing confidence bands for Omega.
#' @param verbose Logical indicating whether to print progress to the command line.
#' @return Named list containing the control options.

causalNullTest.control <- function(mu.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                   g.SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
                                   g.n.bins = 2:(length(unique(A))/50),
                                   cross.fit = TRUE,
                                   V = 10,
                                   folds = NULL,
                                   save.nuis.fits = FALSE,
                                   mu.hat = NULL,
                                   g.hat = NULL,
                                   g.trunc = .001,
                                   n.sim = 1e4,
                                   return.Omega = FALSE,
                                   conf.level = .95,
                                   verbose = FALSE) {
  list(mu.SL.library = mu.SL.library, g.SL.library = g.SL.library, g.n.bins = g.n.bins, cross.fit = cross.fit, V = V,
       folds = folds, mu.hat = mu.hat,  g.hat = g.hat, g.trunc = g.trunc, n.sim = n.sim, return.Omega = return.Omega,
       save.nuis.fits = save.nuis.fits, conf.level = conf.level, verbose = verbose)
}

.check.input <- function(Y, A, W, p, control) {
  if(length(Y) != length(A)) stop("Y and A must have the same length")
  if(!is.null(control$mu.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$mu.hat) < length(unique(control$folds))) {
        stop("mu.hat provided and cross.fit=TRUE, but mu.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("mu.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(!is.null(control$g.hat)) {
    if(control$cross.fit) {
      if(is.null(control$folds)) {
        stop("g.hat provided and cross.fit=TRUE, but folds not provided.")
      }
      if(length(control$g.hat) < length(unique(control$folds))) {
        stop("g.hat provided and cross.fit=TRUE, but g.hats is not a list with the same length as the number of folds.")
      }
    } else {
      if(!is.null(control$folds)) {
        stop("g.hat provided and cross.fit=FALSE, but folds were provided.")
      }
    }
  }
  if(is.null(control$mu.SL.library)) {
    if(is.null(control$mu.hats)) {
      stop("mu.hats must be provided if mu.SL.library is not specified.")
    }
  } else {
    if(is.null(control$mu.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(is.null(control$g.SL.library)) {
    if(is.null(control$g.hats)) {
      stop("g.hats must be provided if g.SL.library is not specified.")
    }
  } else {
    if(is.null(control$g.hats)) {
      if(control$cross.fit & is.null(control$folds) & is.null(control$V)) {
        stop("cross.fit = TRUE, but number of folds not specified.")
      }
    }
  }
  if(any(is.na(Y) | is.na(A) | is.na(W))) {
    stop("Missing outcome, treatment, or confounders detected; missing data not allowed.")
  }
}
