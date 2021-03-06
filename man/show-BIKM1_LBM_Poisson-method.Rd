% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIKM1_LBM_Poisson-class.R
\name{show,BIKM1_LBM_Poisson-method}
\alias{show,BIKM1_LBM_Poisson-method}
\title{Show method for a BIKM1_LBM_Poisson object}
\usage{
\S4method{show}{BIKM1_LBM_Poisson}(object)
}
\arguments{
\item{object}{a BIKM1_LBM_Poisson object}
}
\description{
show method for a \code{\linkS4class{BIKM1_LBM_Poisson}} object
}
\examples{
\donttest{require(bikm1)
J=200
K=120
h=3
l=2
theta=list()
theta$rho_h=1/h *matrix(1,h,1)
theta$tau_l=1/l *matrix(1,l,1)
theta$gamma_hl=matrix(c(1, 6,4, 1, 7, 1),ncol=2)
data=PoissonBlocRnd(J,K,theta)
res=BIKM1_LBM_Poisson(data$x,4,4,4,init_choice='random')
show(res)}

}
