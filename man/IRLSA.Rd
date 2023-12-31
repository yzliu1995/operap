% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{IRLSA}
\alias{IRLSA}
\title{Iterative Re-weighted Least Squares: All Iterations}
\usage{
IRLSA(
  coeff,
  theta,
  cen,
  y,
  Z,
  cov,
  Cnstrn,
  coeffLambda,
  eps = 1e-04,
  maxiter = 10L,
  seed = 0L,
  withCov = FALSE,
  type = "surv"
)
}
\arguments{
\item{coeff}{A vector of initial coefficients for nodes.}

\item{theta}{A vector of initial coefficients for covariates that need adjusting.}

\item{cen}{A vector of censoring statuses or binary outcomes.}

\item{y}{A vector of outcomes.}

\item{Z}{A matrix of nodes.}

\item{cov}{A matrix of covariates.}

\item{Cnstrn}{A matrix of constraints.}

\item{coeffLambda}{A vector of tuning parameter lambda.}

\item{eps}{A numeric value denoting the accuracy level.}

\item{maxiter}{An integer denoting the maximum number of iterations to obtain the MLEs.}

\item{seed}{The seed used for generating the simulation dataset.}

\item{withCov}{Whether any covariates need adjusting.}

\item{type}{A string denoting the outcome type, either "survival" or "binary" outcome.}
}
\value{
The estimates
}
\description{
This function performs all iterations of the Iterative Re-weighted Least Squares (IRLS) algorithm to obtain the Maximum Likelihood Estimates (MLEs) based on the given parameters.
}
