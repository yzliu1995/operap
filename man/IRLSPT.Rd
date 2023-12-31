% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{IRLSPT}
\alias{IRLSPT}
\title{Perform Parameter Tuning for Each Iteration of Iterative Re-weighted Least Squares in the Pruning Step}
\usage{
IRLSPT(
  coeff,
  theta,
  cen,
  y,
  Z,
  cov,
  LCnstrn,
  qCnstrn,
  mCnstrn,
  minCoeff,
  maxCoeff,
  eps = 1e-04,
  maxiter = 10L,
  withCov = FALSE,
  seed = 0L,
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

\item{LCnstrn}{A matrix of linear partial ordering constraints.}

\item{qCnstrn}{A matrix of quadratic constraints.}

\item{mCnstrn}{An integer denoting how many more constraints are needed.}

\item{minCoeff}{The lower bound for all the coefficients.}

\item{maxCoeff}{The upper bound for all the coefficients.}

\item{eps}{A numeric value denoting the accuracy level.}

\item{maxiter}{An integer denoting the maximum number of iterations.}

\item{withCov}{Whether any covariates need adjusting.}

\item{seed}{The seed used for generating the simulation dataset.}

\item{type}{A string indicating the type of outcome, either "survival" or "binary" outcome.}
}
\value{
The coefficients after tuning the parameter.
}
\description{
This function performs parameter tuning for each iteration of the Iterative Re-weighted Least Squares (IRLS) algorithm in the pruning step based on the given parameters.
}
