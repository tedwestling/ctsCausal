% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixed.dens.R
\name{cmdSuperLearner.onebin}
\alias{cmdSuperLearner.onebin}
\title{cmdSuperLearner for a specific number of bins}
\usage{
cmdSuperLearner.onebin(A, W, b, SL.library, verbose, validRows)
}
\arguments{
\item{A}{\code{n x 1} numeric vector of exposure values.}

\item{W}{\code{n x p} data.frame of covariate values to condition upon.}

\item{b}{Integer number of bins >= 2.}

\item{SL.library}{Library to use for bin-specific probabilities.}

\item{verbose}{Logical indicating whether to print progress reports to the command line.}

\item{validRows}{List of rows in each CV fold.}
}
\value{
Returns a named list with the following elements:
\item{bins}{List of length \code{b} containing the sets used for each bin.}
\item{bin.fits}{List of length \code{b} containing the estimated SuperLearner objects for each bin.}
\item{a.ecdf}{Empirical CDF of the exposure.}
\item{SL.bin.probs}{SuperLearner conditional probabilities of being in each bin.}
\item{SL.densities}{SuperLearner conditional standardized mixed density correspondint to each bin.}
\item{cv.library.densities}{Cross-validated library conditional standardized mixed density corresponding to each bin.}
\item{library.densities}{Library conditional standardized mixed density corresponding to each bin fit on the full data.}
\item{alg.names}{Algorithm names.}
}
\description{
This function estimates the conditional mixed density using a given number of bins \code{b}.
}
\examples{
# Define parameters
n <- 300
W <- data.frame(matrix(rnorm(3 * n), ncol = 3))
Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-W[,1] + W[,2])))
A <- (1-Z) * rnorm(n, mean = W[,2] - W[,3], sd = abs(1 + W[,1]))
validRows <- cmdCVFolds(n = n, cvControl = list(V = 10, shuffle=TRUE, validRows = NULL))
bin.fit <- cmdSuperLearner.onebin(A, W, b = 2, SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"), verbose=TRUE, validRows = validRows)
}
