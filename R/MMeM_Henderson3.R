library(MASS)
library(Matrix)
library(psych)


###HendersonIII

MMeM_henderson3 <- function(fml, data, factor_X){


  data_matrix = MMeM_terms(fml ,data, factor_X = factor_X)
  X = data_matrix$X
  Y = data_matrix$Y
  Z = data_matrix$Z
  N = data_matrix$N
  I = data_matrix$I
  q = data_matrix$q
  DV = data_matrix$DV

  HX<-X%*%ginv(t(X)%*%X)%*%t(X)
  XZ<- cbind(X,Z)
  HXZ<-XZ%*%ginv(t(XZ)%*%XZ)%*%t(XZ)
  C_E = (I - HXZ)
  C_D = HXZ - HX
  C_T = C_D-C_E*psych::tr(C_D)/psych::tr(C_E)



  E_aov<-(t(Y)%*%C_E%*%Y)/psych::tr(C_E)
  cov_E = (2/psych::tr(C_E))*(E_aov%x%E_aov)
  v_E<- unique(diag(cov_E))
  E.estimate = rbind(E_aov[upper.tri(E_aov,diag=TRUE)], v_E)

  T_aov = (t(Y)%*%C_T%*%Y)/psych::tr(Z%*%t(Z)%*%(I-HX))
  D_aov = (t(Y)%*%C_D%*%Y)
  cov_D = (2/psych::tr(C_D))*(D_aov%x%D_aov)

  cov_T<- (cov_D - cov_E*(psych::tr(C_D))^2)/(psych::tr(Z%*%t(Z)%*%(I-HX)))^2
  v_T = unique(diag(cov_T))
  T.estimate = rbind(T_aov[upper.tri(T_aov,diag=TRUE)], v_T)

  Tnames = c()
  Enames = c()
  for(i in 1:q){
    for(j in i:q){
      Y_names = paste(DV[i:j], collapse = " ")
      Tnames = c(Tnames, paste('T:', Y_names))
      Enames = c(Enames, paste('E:', Y_names))
    }
  }
  w = q*(q+1)
  names = double(w)
  names[odd(1:w)] = Tnames
  names[even(1:w)] = Enames
  rownames(Vcov) = names
  colnames(Vcov) = names

  return(list(T.estimate = , E.estimate = E_aov[upper.tri(E_aov,diag=TRUE)], VCOV = Vcov))
}
