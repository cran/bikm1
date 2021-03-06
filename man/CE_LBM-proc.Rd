% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bikm1_MLBM_Binary_functions.R
\name{CE_LBM}
\alias{CE_LBM}
\title{CE_LBM function for agreement between co-clustering partitions}
\usage{
CE_LBM(v,w,vprime,wprime)
}
\arguments{
\item{v}{numeric vector  specifying the class of  rows.}

\item{w}{numeric vector  specifying the class of  columns.}

\item{vprime}{numeric vector  specifying another partition of  rows.}

\item{wprime}{numeric vector  specifying another partition of  columns.}
}
\value{
ce_vw: the value of the index (between 0 and 1). A value of 0 corresponds to a perfect match.
}
\description{
Produce a measure of agreement between two pairs of partitions for co-clustering using CE_simple on columns and rows of a matrix. A value of 1 means a perfect match.
}
\examples{
\donttest{
require(bikm1)
set.seed(42)
v=floor(runif(4)*2)
vprime=floor(runif(4)*2)
w=floor(runif(4)*3)
wprime=floor(runif(4)*3)
error=CE_LBM(v,w,vprime,wprime)
}
}
