% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main_function.R
\name{makePositiveDefinite}
\alias{makePositiveDefinite}
\title{Positive definite matrix}
\usage{
makePositiveDefinite(M1, M2 = NULL, threshold = 0.1)
}
\arguments{
\item{M1, M2}{A squared numeric matrix, typically a correlation or a covariance matrix. It must be symmetric.}

\item{threshold}{A correction factor.}
}
\value{
A list with the corrected input matrices and the correction \code{threshold}-lambda.
}
\description{
Adjust the diagonal of a symmetric square matrix, by the smallest eigenvalue method, in order to make it positive definite.
}
\details{
Finds the smallest eigenvalue lambda of \code{M1} (or \code{M1} and \code{M2} if supplied) and adds (threshold-lambda) to the diagonal to make it positive definite.
}
