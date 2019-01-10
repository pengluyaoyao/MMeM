#' MMeM: Estimating the variance covariance components of the multivariate mixed effects model
#'
#'
#' @author Luyao Peng \email{luyaopeng.cn@gmail.com}
#' @author Rui Yang \email{ray.cn.us@gmail.com}
#' @description
#' This package analyzes data under multivariate mixed effects model using multivariate REML and multivariate Henderson3 methods.
#' Currently, it only supports multivariate mixed effects model with one fixed effects and one random effects and two response variates.
#' @references
#' \itemize{
#' \item Meyer, K. (1985). Maximum Likelihood Estimation of Variance Components for a Multivariate Mixed Model with Equal Design Matrices.
#'  Biometrics, 41(1), 153-165. <doi:10.2307/2530651>
#' \item Wesolowska Janczarek, M. T. (1984). "Estimation of covariance matrices in unbalanced random and mixed multivariate models."
#' Biometrical journal 26(6), 665,674. <doi:10.1002/bimj.4710260613>
#' }
#' @docType package
#' @name MMeM
"_PACKAGE"

#' simulated datasets
#' @title
#' simulated bivariate data
#'
#' @description
#' This is a simulated data with 2 dependent variables and one fixed effects
#' and one random effects
#'
#' @docType data
#' @name simdata
#' @usage data(simdata)
#'
#'
#' @keywords datasets
#'
NULL
