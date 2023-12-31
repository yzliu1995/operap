% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OPERA.R
\name{as.param}
\alias{as.param}
\title{Constructing an object having the coefficients}
\usage{
as.param(beta.vec, lambda = NULL)
}
\arguments{
\item{beta.vec}{The values of the coefficients. If there is a single tuning parameter, it should be a vector. If there are multiple tuning parameters, it should be a matrix, with each row associated with one tuning parameter.}

\item{lambda}{The tuning parameter associated with this beta. It can be a scalar or a vector.}
}
\value{
The object that stores the coefficients.
}
\description{
This function constructs an object called "param" used for better storage, printing, and visualization of the coefficients.
}
