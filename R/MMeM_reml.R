##multivariate REML

#' finding Q and checking Q
#'
#' @param T
#' @param E
#'
#' @example
#'
#' @import matrixcalc
#' @import jointDiag
#'
#' @export
find_Q <- function(T,E, q){
  C = array(NA, dim=c(q,q,2))

  #  if(det(T) <= 0){
  #   T = as.matrix(nearPD(T)$mat)
  # }
  #   if(det(E) <= 0){
  #   E = as.matrix(nearPD(E)$mat)
  # }

  C[,,1] = T
  C[,,2] = E
  Q_est = jointDiag::jadiag(C)$B
  Q_adj = Q_est%*%C[,,2]%*%t(Q_est)
  Q_adj[upper.tri(Q_adj)] = 0
  Q_adj[lower.tri(Q_adj)]= 0
  Q_est_adj = solve(sqrt(Q_adj))%*%Q_est
  return(list(Q = Q_est_adj))#, T=T, E=E))
}

check_Q <- function(Q,T,E){
  QTQ = Q%*%T%*%t(Q)
  QEQ = Q%*%E%*%t(Q)
  if(matrixcalc::is.diagonal.matrix(QTQ) == TRUE & matrixcalc::is.diagonal.matrix(QEQ) == TRUE){
    print('----- Q is good -----')
  }else{
    print('----- Q is not good -----')
    return('fail')
  }
}

odd <- function(x) x%%2 != 0
even <- function(x) x%%2 == 0
##### Main function of multivariate mixed effects model_REML#####

MMeM_reml <- function(fml, data, factor_X, T.start, E.start, maxit=50, tol = 0.000000001){

  T = T.start
  E = E.start

  data_matrix = MMeM_terms(fml ,data, factor_X = factor_X)
  X = data_matrix$X
  Y = data_matrix$Y
  Z = data_matrix$Z
  N = data_matrix$N
  I = data_matrix$I
  q = data_matrix$q

  if(q >1){
    if(det(T.start)<=0 | det(E.start)<=0){
      stop('Initial variance-covariance matrix should be positive definite')
    }
    if(isSymmetric(T.start) == FALSE | isSymmetric(E.start) == FALSE){
      stop('Initial variance-covariance matrix should be symmetric')
    }
  }

  ###ZHZ
  H = I-X%*%ginv(t(X)%*%X)%*%t(X)
  h = ncol(X)
  s = ncol(Z)
  ZH = t(Z)%*%H
  ZHZ = t(Z)%*%H%*%Z

  for(iter in 1:maxit){
    print(paste('iteration:',iter))

    Q_results = find_Q(T, E, q)
    Q = Q_results$Q

    Check_Q = check_Q(Q, T, E)
    if(Check_Q =='fail'){
      break
    }

    Q_inv = solve(Q)

    ####ZHy, Lambda, yHy, traces of C and CC, u_c
    Lambda= diag(Q%*%T%*%t(Q))

    Cc = list()
    Cc_trace = double(q)
    for(i in 1:q){
      Cc[[i]] = solve(ZHZ+(1/Lambda[i])*diag(s))
      Cc_trace[i] = tr(Cc[[i]])
    }

    CcCc_trace = matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        CcCc_trace[i,j] = tr(Cc[[i]]%*%Cc[[j]])
      }
    }

    ZHy_c = list()
    y_c = list()
    u_c = list()
    for(i in 1:q){
      if(q == 1){
        y_c[[i]] = Q[i,1]*Y[,1]
      }else{
        y_c[[i]] = Q[i,1]*Y[,1]+ Q[i,2]*Y[,2]
      }
      ZHy_c[[i]] = ZH%*%y_c[[i]]
      u_c[[i]] = Cc[[i]]%*%ZHy_c[[i]]
    }

    ycHyc = matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        ycHyc[i,j] = t(y_c[[i]])%*%H%*%y_c[[j]]
      }
    }

    ###REML equations
    BTTc=matrix(0,q,q)
    BTEc=matrix(0,q,q)
    BEEc=matrix(0,q,q)
    dTc=matrix(0,q,q)
    dEc=matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        BTTc[i,j] = 1/(Lambda[i]*Lambda[j]) * (s - ((1/Lambda[i])*Cc_trace[i]+(1/Lambda[j])*Cc_trace[j])+ (1/(Lambda[i]*Lambda[j]))*CcCc_trace[i,j])
        BTEc[i,j] = 0.5*(1/(Lambda[i]*Lambda[j])) *(Cc_trace[i]+Cc_trace[j]-(1/Lambda[i]+1/Lambda[j])*CcCc_trace[i,j])
        BEEc[i,j] = (N-h-s)+(1/(Lambda[i]*Lambda[j])) * CcCc_trace[i,j]
        dTc[i,j] = 1/(Lambda[i]*Lambda[j]) * (t(u_c[[i]])%*%u_c[[j]])
        dEc[i,j] = ycHyc[i,j]-t(u_c[[i]])%*%ZHy_c[[j]]-(1/Lambda[i])*t(u_c[[i]])%*%u_c[[j]]
      }
    }

    oldT = T
    oldE = E

    Bc = bdiag()
    dc = c()
    for(i in 1:q){
      for(j in i:q){
        if(i == j){
          Bc = bdiag(Bc, matrix(c(BTTc[i,j],BTEc[i,j],BTEc[i,j], BEEc[i,j]),2,2))
          dc = c(dc, c(dTc[i,j], dEc[i,j]))
        }else{
          Bc = bdiag(Bc, 2*matrix(c(BTTc[i,j],BTEc[i,j],BTEc[i,j], BEEc[i,j]),2,2))
          dc = c(dc, 2*c(dTc[i,j], dEc[i,j]))
        }

      }
    }
    thetas <- solve(Bc)%*%as.matrix(dc)

    T.new = Matrix(0,q,q)
    T.new[upper.tri(T.new,diag=TRUE)] =thetas[odd(1:length(thetas))]
    T.new = forceSymmetric(T.new, uplo = 'U')
    E.new = Matrix(0,q,q)
    E.new[upper.tri(E.new,diag=TRUE)] =thetas[even(1:length(thetas))]
    E.new = forceSymmetric(E.new, uplo = 'U')

    T = as.matrix(Q_inv%*%T.new%*%t(Q_inv))
    E = as.matrix(Q_inv%*%E.new%*%t(Q_inv))

    if(max(abs(T-oldT)) < sqrt(tol) & max(abs(E-oldE))< sqrt(tol)){
      break
    }
  }

  if(iter == maxit){
    message('Note: convergence is not reached, results may not be valid.')
  }

  m = 0
  deriv = c()
  for(i in 1:q){
    for(j in i:q){
      m = m+1
      D = matrix(0,q,q)
      D[i,j] = 1
      D = forceSymmetric(D, uplo = 'U')
      deriv[m] = (Q%*%D%*%t(Q))[i,j]
    }
  }

  Vcov = solve(Bc%*%diag(rep(deriv^2,each = 2)))*2

  Tnames = c()
  Enames = c()
  for(i in 1:q){
    for(j in i:q){
      Y_names = paste(colnames(data)[i:j], collapse = " ")
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

  return(list(T.estimates = T[upper.tri(T,diag=TRUE)], E.estimates = E[upper.tri(E, diag = TRUE)], VCOV = Vcov))
}



