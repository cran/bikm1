##' bikm1 package
##'
##' This package is designed to co-cluster a contingency (resp. binary) matrix, or  double binary matrices in blocks respectively under the (normalized or not) Poisson (resp binary) Latent Block Model and the Multiple Latent Block Model. It enables to automatically select the number of row and column clusters and to compare partition estimations with reference partitions.
##'
##' @section Features: Package for the segmentation of the rows and columns inducing a co-clustering and automatically select the number of row and column clusters.
##'
##' @section Model 1 : \code{\link{BIKM1_LBM_Poisson}} . This fitting procedure produces a \code{\linkS4class{BIKM1_LBM_Poisson}} object.
##' @section Model 2 : \code{\link{BIKM1_LBM_Binary}} . This fitting procedure produces a \code{\linkS4class{BIKM1_LBM_Binary}} object.
##' @section Model 3: \code{\link{BIKM1_MLBM_Binary}} . This fitting procedure produces a \code{\linkS4class{BIKM1_MLBM_Binary}} object.
##'@section Technical remarks: Display of the result with \code{\link{plot,BIKM1_LBM_Poisson-method}} and
##'
##' with \code{\link{show,BIKM1_LBM_Poisson-method}}, with \code{\link{summary,BIKM1_LBM_Poisson-method}} and with \code{\link{print,BIKM1_LBM_Poisson-method}}.
##'
##'Display of the result with \code{\link{plot,BIKM1_LBM_Binary-method}} and
##'
##' with \code{\link{show,BIKM1_LBM_Binary-method}}, with \code{\link{summary,BIKM1_LBM_Binary-method}} and with \code{\link{print,BIKM1_LBM_Binary-method}}.
##'
##'Display of the result with \code{\link{plot,BIKM1_MLBM_Binary-method}} and
##'
##' with \code{\link{show,BIKM1_MLBM_Binary-method}}, with \code{\link{summary,BIKM1_MLBM_Binary-method}} and with \code{\link{print,BIKM1_MLBM_Binary-method}}.
##' @name bikm1-package
##' @docType package
##' @author Valerie Robert \email{valerie.robert.math@@gmail.com}
##' @references  Keribin, Celeux and Robert, The Latent Block Model: a useful model for high dimensional data. https://hal.inria.fr/hal-01658589/document
##'
##'   Govaert and Nadif. Co-clustering, Wyley (2013).
##'
##' Keribin, Brault and Celeux. Estimation and Selection for the Latent Block Model on Categorical Data, Statistics and Computing (2014).
##'
##' Robert. Classification croisee pour l'analyse de bases de donnees de grandes dimensions de pharmacovigilance. Thesis, Paris Saclay (2017).
##'
##' Robert, Vasseur and Brault. Comparing high dimensional partitions with the Co-clustering Adjusted Rand Index, Journal of Classification, 38(1), 158-186 (2021).
##'
##' @import ade4
##' @import pracma
##' @import methods
##' @import ggplot2
##' @import reshape2
##' @import grid
##' @importFrom lpSolve lp.assign
##' @importFrom parallel mclapply detectCores
##' @importFrom grDevices dev.new
##' @importFrom graphics boxplot lines par title mtext
##' @importFrom stats rgamma rpois runif rbinom
##' @importFrom gtools permutations
##'
NULL


