% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bikm1_LBM_Binary_functions.R
\name{BinBlocICL_LBM}
\alias{BinBlocICL_LBM}
\title{BinBlocICL_LBM function  for computation of the ICL criterion in the Binary LBM}
\usage{
BinBlocICL_LBM(a,b,x,z1,v1)
}
\arguments{
\item{a}{an hyperparameter for priors on the mixing proportions. By default, a=4.}

\item{b}{an hyperparameter for prior on the Bernoulli parameter. By default, b=1.}

\item{x}{contingency matrix of observations.}

\item{z1}{a numeric vector specifying the class of rows.}

\item{v1}{a numeric vector specifying the class of columns.}
}
\value{
a value of the ICL criterion.
}
\description{
Produce a value of the ICL criterion in the Binary LBM.
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
BinBlocICL_LBM(a=4,b=1,data$x, data$xrow,data$xcol)
}
