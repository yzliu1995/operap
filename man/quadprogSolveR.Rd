% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{quadprogSolveR}
\alias{quadprogSolveR}
\title{C++ function to call the R function \code{\link[quadprog:solve.QP]{quadprog::solve.QP()}}}
\usage{
quadprogSolveR(Dmat, dvec, Amat)
}
\arguments{
\item{Dmat}{A matrix appearing in the quadratic function to be minimized.}

\item{dvec}{A vector appearing in the quadratic function to be minimized.}

\item{Amat}{A matrix defining the constraints under which we want to minimize the quadratic function.}
}
\value{
A list containing the solution obtained from the R function \code{\link[quadprog:solve.QP]{quadprog::solve.QP()}}.
}
\description{
This function serves as a bridge between C++ and R, allowing you to call the R function \code{\link[quadprog:solve.QP]{quadprog::solve.QP()}} from C++ code.
}
\seealso{
\code{\link[quadprog:solve.QP]{quadprog::solve.QP()}}
}
