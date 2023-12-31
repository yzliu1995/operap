% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OPERA.R
\name{edgeMisclassfication}
\alias{edgeMisclassfication}
\title{Calculating the edge misclassification rate}
\usage{
edgeMisclassfication(edgeMatrix, trueS, estimateS)
}
\arguments{
\item{edgeMatrix}{A vector of partial ordering relationship node variables.}

\item{trueS}{The true stages.}

\item{estimateS}{The estimated stages.}
}
\value{
The edge misclassification rate.
}
\description{
The function computes the edge misclassification rate, which measures the accuracy of the estimated edges compared to the true edges. It quantifies the proportion of misclassified edges in relation to the total number of edges.
}
