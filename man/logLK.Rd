% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{logLK}
\alias{logLK}
\title{Calculate Negative Log (Partial) Likelihood}
\usage{
logLK(Z, cen, beta, y, cov, theta, type = "surv", withCov = FALSE)
}
\arguments{
\item{Z}{A matrix of nodes.}

\item{cen}{A vector of censoring statuses or binary outcomes.}

\item{beta}{A vector of node coefficients.}

\item{y}{A vector of outcomes.}

\item{cov}{A matrix of covariates.}

\item{theta}{A vector of covariate coefficients.}

\item{type}{A string indicating the type of outcome, currently either "survival" or "binary".}

\item{withCov}{Whether any covariates need adjusting.}
}
\value{
The negative log likelihood.
}
\description{
This function calculates the negative log (partial) likelihood based on the given parameters.
}
