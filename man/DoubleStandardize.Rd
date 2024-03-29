% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Preliminary_Functions.R
\name{DoubleStandardize}
\alias{DoubleStandardize}
\title{Function that double standardizes a matrix}
\usage{
DoubleStandardize(X, tol = .Machine$double.eps^0.5, maxIter = 500)
}
\arguments{
\item{X}{Original matrix to be double-standardized}

\item{tol}{Tolerance. Default is square root of machine precision}

\item{maxIter}{Maximum iteration. Default is 500}
}
\value{
A list that contains:
\item{Result}{The matrix after double standardization.}
\item{Iter}{Number of iterations the function runs.}
}
\description{
Function that double standardizes a matrix
}
\details{
After double standardization, the matrix will have mean zero and (population) variance one for each row and column.
}
