% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OPERA.R
\name{summaryC}
\alias{summaryC}
\title{Calculating the Kaplan-Meier survival probabilities}
\usage{
summaryC(obj, times)
}
\arguments{
\item{obj}{A \code{survfit} object representing the survival data.}

\item{times}{A vector of time points at which the survival probabilities are to be computed.}
}
\value{
The K-M survival probabilities for the specified \code{times}.
}
\description{
The function computes the K-M survival probabilities for the specified \code{times}. The K-M survival probabilities estimate the probability of survival at each time point based on the observed survival data.
}
