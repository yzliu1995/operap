% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OPERA.R
\name{Mode}
\alias{Mode}
\title{Calculating the mode}
\usage{
Mode(x)
}
\arguments{
\item{x}{A vector.}
}
\value{
The mode of the input vector.
}
\description{
This function takes a vector as input and calculates the mode, which is the most frequently occurring value in the vector. The function achieves this by utilizing two auxiliary functions: \code{\link[base:tabulate]{base::tabulate()}}, which calculates the frequency of each unique value in the vector, and \code{\link[base:match]{base::match()}}, which returns the position of each unique value when matched.
}
\examples{
Mode(c(1, 2, 3, 2, 4, 2))
# Output: The mode is 2

}
\seealso{
\code{\link[base:tabulate]{base::tabulate()}}: Calculates the frequency for each unique value in a vector.

\code{\link[base:match]{base::match()}}: Returns the position of each unique value when matched in a vector.

Other functions for finding statistical measures: \code{\link[base:mean]{base::mean()}}, \code{\link[base:Extremes]{base::min()}}, \code{\link[base:Extremes]{base::max()}}, \code{\link[stats:median]{stats::median()}}
}
