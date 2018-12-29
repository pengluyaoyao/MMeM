###HendersonIII

#' Multivariate Henderson3 method
#' @param fml two-sided linear formula object describing both the fixed-effects and random-effects  parts  of  the  model,
#'  with  the  response  on  the  left  of  a ~ operator. For univariate response, put variable name directly; for multivariate responses
#'  combine variables using concatenate operator, for example, for bivariate responses, c(var1, var2).  The predictor terms are separated  by + operators,  on  the  right.   Random-effects  terms  are
#' distinguished by vertical bars '|' separating expressions for design matrices from grouping factors.
#' @param data data frame containing the variables named in formula.
#' @param factor_X (logical) indicating whether predictor is a factor or continuous. By default is TRUE
#' @examples
#' \dontrun{
#' data(simdata)
#' results_henderson <- MMeM_henderson3(fml = c(V1,V2) ~ X_vec + (1|Z_vec),
#' data = simdata, factor_X = TRUE)
#' }
#' @return \code{T.estimates} is the estimated matrix of the variance covariance matrix of the block random effects with corresponding standard errors;
#' \code{E.estimates} is the estimated matrix of the variance covariance matrix of the residuals with corresponding standard errors;
#'
#' @references Wesolowska‐Janczarek, M. T. "Estimation of covariance matrices in unbalanced random and mixed multivariate models." Biometrical journal 26.6 (1984): 665-674.
#' @importFrom psych tr
#' @importFrom stringr str_extract_all
#' @export
MMeM_henderson3 <- function(fml, data, factor_X){
  data_matrix = MMeM_terms(fml ,data, factor_X = factor_X)
  X = data_matrix$X
  Y = data_matrix$Y
  Z = data_matrix$Z
  N = data_matrix$N
  I = data_matrix$I
  q = data_matrix$q
  DV = data_matrix$DV

  HX<-X%*%MASS::ginv(t(X)%*%X)%*%t(X)
  XZ<- cbind(X,Z)
  HXZ<-XZ%*%MASS::ginv(t(XZ)%*%XZ)%*%t(XZ)
  C_E = (I - HXZ)
  C_D = HXZ - HX
  C_T = C_D-C_E*psych::tr(C_D)/psych::tr(C_E)



  E_aov<-(t(Y)%*%C_E%*%Y)/psych::tr(C_E)
  cov_E = (2/psych::tr(C_E))*(E_aov%x%E_aov)
  E.variance<- unique(diag(cov_E))
  E.estimates = E_aov[upper.tri(E_aov,diag=TRUE)]
  E_results = rbind(E.estimates, E.variance)

  T_aov = (t(Y)%*%C_T%*%Y)/psych::tr(Z%*%t(Z)%*%(I-HX))
  D_aov = (t(Y)%*%C_D%*%Y)
  cov_D = (2/psych::tr(C_D))*(D_aov%x%D_aov)

  cov_T<- (cov_D + cov_E*(psych::tr(C_D))^2)/(psych::tr(Z%*%t(Z)%*%(I-HX)))^2
  T.variance = unique(diag(cov_T))
  T.estimates = T_aov[upper.tri(T_aov,diag=TRUE)]
  T_results = rbind(T.estimates, T.variance)

  Tnames = c()
  Enames = c()
  for(i in 1:q){
    for(j in i:q){
      Y_names = paste(DV[i:j], collapse = " ")
      Tnames = c(Tnames, paste('T:', Y_names))
      Enames = c(Enames, paste('E:', Y_names))
    }
  }

  colnames(T_results) = Tnames
  colnames(E_results) = Enames

  return(list(T.estimates = T_results, E.estimates = E_results))
}
##utility functions
MMeM_terms <- function(fml, data, factor_X){

  Fml = stats::formula(fml)
  any_RE <- length( lme4::findbars(Fml))

  if(any_RE == 0){
    stop('No random effects in the model')
  } else{
    df = stats::get_all_vars(Fml, data)

    terms = attr(stats::terms.formula(Fml), 'variables')
    DVs = all.names(terms[2])
    if(length(DVs) == 3){
      message(paste('Bivariate response:', DVs[2], 'and', DVs[3]))
      Y = as.matrix(df[,match(DVs[2:length(DVs)], colnames(df))])
      DV = c(DVs[2], DVs[3])
    }else if(length(DVs) == 1){
      message(paste('Univariate response:', DVs[1:length(DVs)]))
      Y = as.matrix(df[,match(DVs[1:length(DVs)], colnames(df))])
      DV = DVs[1]
    }else{
      stop('Dependent variables should be univariate or bivariate.')
    }

    re_term = stringr::str_extract_all(format(terms[length(terms)]), "(?<=\\|).+?(?=\\))")[[1]]
    re_strip = gsub(" ", "", re_term, fixed = TRUE)
    re_data = df[,match(re_strip, colnames(df))]
    Z = stats::model.matrix(~ -1 + factor(re_data))

    IV = all.names(terms[-c(1,2,length(terms))])
    IV_data = df[,match(IV, colnames(df))]
    if(factor_X == TRUE){
      X = stats::model.matrix(~ factor(IV_data))
    }else{
      X = stats::model.matrix(~ IV_data)
    }


    N = nrow(X)
    I = diag(1, N)
    q = ncol(as.matrix(Y))
  }

  return(list(X = X, Z = Z, Y = Y, N =N, I =I, q=q, DV = DV))
}





