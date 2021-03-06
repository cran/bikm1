% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bikm1_LBM_Binary_functions.R
\name{BinBlocVisuResum_LBM}
\alias{BinBlocVisuResum_LBM}
\title{BinBlocVisuResum_LBM function  for visualization of binary matrix data-sets}
\usage{
BinBlocVisuResum_LBM(x,z,v)
}
\arguments{
\item{x}{binary matrix of observations.}

\item{z}{a numeric vector specifying the class of rows.}

\item{v}{a numeric vector specifying the class of columns.}
}
\value{
a \pkg{plot} object.
}
\description{
Produce a plot object representing the resumed co-clustered data-sets.
}
\examples{
require(bikm1)
set.seed(42)
n=200
J=120
g=3
h=2
theta=list()
theta$pi_g=1/g *matrix(1,g,1)
theta$rho_h=1/h *matrix(1,h,1)
theta$alpha_gh=matrix(runif(6),ncol=h)
data=BinBlocRnd_LBM(n,J,theta)
BinBlocVisuResum_LBM(data$x,data$xrow,data$xcol)
}
