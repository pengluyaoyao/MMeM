##multivariate REML

#' Multivariate REML Method
#' @description Estimating the variance components under the multivariate mixed effects model using REML methods
#' @param fml a two-sided linear formula object describing both the fixed-effects and random-effects  parts  of  the  model,
#'  with  the  response  on  the  left  of  a ~ operator. For univariate response, put variable name directly; for multivariate responses
#'  combine variables using concatenate operator, for example, for bivariate responses, c(var1, var2).  The predictor terms are separated  by + operators,  on  the  right.   Random-effects  terms  are
#' distinguished by vertical bars '|' separating expressions for design matrices from grouping factors.
#' @param data data frame containing the variables named in formula.
#' @param factor_X (logical) indicating whether predictor is a factor or continuous. By default is TRUE
#' @param T.start the starting matrix for the variance covariance matrix of the block random effects, it has to be positive definite q by q symmetric matrix.
#' @param E.start the starting matrix for the variance covariance matrix of the block random effects, it has to be positive definite q by q symmetric matrix.
#' @param maxit the maximum number of iterations
#' @param tol the convergence tolerance
#' @details
#' Suppose n observational units, q variates, p fixed effects coefficients and s random effects units.
#' The model supports multivariate mixed effects model for one-way randomized block design with equal design matrices:
#' \deqn{Y = XB + ZU + E}
#' where Y is n by q reponse matrix;
#' X is n by p design matrix for the fixed effects;
#' B is p by q coefficients matrix for the fixed effects;
#' Z is n by s design matrix for the random effects;
#' U is s by q matrix for the random effects;
#' E is n by q random errors matrix.
#'
#' The model also supports simple OLS multivariate regression:
#' \deqn{y = Xb + Zu + e}
#' where y is n by 1 response vector;
#' b is p by 1 coefficients vector for the fixed effects;
#' u is s by 1 matrix for the random effects.
#'
#' @return the output \code{T.estimates} is the estimated matrix of the variance covariance matrix of the block random effects;
#' \code{E.estimates} is the estimated matrix of the variance covariance matrix of the residuals; \code{VCOV} is the asymptotic
#' dispersion matrix of the estimated variance covariance components.
#' @examples
#' \dontrun{
#' data(simdata)
#' T.start <- matrix(c(10,5,5,15),2,2)
#' E.start <- matrix(c(10,1,1,3),2,2)
#' results_reml <- MMeM_reml(fml = c(V1,V2) ~ X_vec + (1|Z_vec), data = simdata,
#' factor_X = TRUE, T.start = T.start, E.start = E.start, maxit = 10)}
#' }
#'
#' @references Meyer, K. A. R. I. N. "Maximum likelihood estimation of variance components for a multivariate mixed model with equal design matrices." Biometrics 1985: 153-165.
#'
#' @importFrom matrixcalc is.diagonal.matrix
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom jointDiag jadiag
#' @importFrom psych tr
#' @importFrom Matrix bdiag
#' @importFrom Matrix Matrix
#' @importFrom stringr str_extract_all
#' @export
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
  DV = data_matrix$DV

  if(q >1){
    if(det(T.start)<=0 | det(E.start)<=0){
      stop('Initial variance-covariance matrix should be positive definite')
    }
    if(isSymmetric(T.start) == FALSE | isSymmetric(E.start) == FALSE){
      stop('Initial variance-covariance matrix should be symmetric')
    }
  }

  ###ZHZ
  H = I-X%*%MASS::ginv(t(X)%*%X)%*%t(X)
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
      Cc_trace[i] = psych::tr(Cc[[i]])
    }

    CcCc_trace = matrix(0,q,q)
    for(i in 1:q){
      for(j in 1:q){
        CcCc_trace[i,j] = psych::tr(Cc[[i]]%*%Cc[[j]])
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

    Bc = Matrix::bdiag()
    dc = c()
    for(i in 1:q){
      for(j in i:q){
        if(i == j){
          Bc = Matrix::bdiag(Bc, matrix(c(BTTc[i,j],BTEc[i,j],BTEc[i,j], BEEc[i,j]),2,2))
          dc = c(dc, c(dTc[i,j], dEc[i,j]))
        }else{
          Bc = Matrix::bdiag(Bc, 2*matrix(c(BTTc[i,j],BTEc[i,j],BTEc[i,j], BEEc[i,j]),2,2))
          dc = c(dc, 2*c(dTc[i,j], dEc[i,j]))
        }

      }
    }

    if(matrixcalc::is.singular.matrix(as.matrix(Bc))){
      stop('Information matrix is not invertible, please increase the levels of the rando effects')
    }else{
      thetas <- solve(Bc)%*%as.matrix(dc)
    }


    T.new = matrix(0,q,q)
    T.new[upper.tri(T.new,diag=TRUE)] =thetas[odd(1:length(thetas))]
    T.new = Matrix::forceSymmetric(T.new, uplo = 'U')
    E.new = matrix(0,q,q)
    E.new[upper.tri(E.new,diag=TRUE)] =thetas[even(1:length(thetas))]
    E.new = Matrix::forceSymmetric(E.new, uplo = 'U')

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
      D = Matrix::forceSymmetric(D, uplo = 'U')
      deriv[m] = (Q%*%D%*%t(Q))[i,j]
    }
  }

  Vcov = Matrix::Matrix(solve(Bc%*%diag(rep(deriv^2,each = 2)))*2, sparse = TRUE)

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

  return(list(T.estimates = T[upper.tri(T,diag=TRUE)], E.estimates = E[upper.tri(E, diag = TRUE)], VCOV = Vcov))
}
## utility functions
find_Q <- function(T,E,q){
  C = array(NA, dim=c(q,q,2))
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

MMeM_terms <- function(fml, data, factor_X){

  Fml = formula(fml)
  any_RE <- length( lme4::findbars(Fml))

  if(any_RE == 0){
    stop('No random effects in the model')
  } else{
    df = get_all_vars(Fml, data)

    terms = attr(terms.formula(Fml), 'variables')
    DVs = all.names(terms[2])
    if(length(DVs) == 3){
      message(paste('Bivariate response:', DVs[2], 'and', DVs[3]))
      Y = as.matrix(df[,match(DVs[2:length(DVs)], colnames(df))])
    }else if(length(DVs) == 1){
      message(paste('Univariate response:', DVs[1:length(DVs)]))
      Y = as.matrix(df[,match(DVs[1:length(DVs)], colnames(df))])
    }else{
      stop('Dependent variables should be univariate or bivariate.')
    }

    re_term = stringr::str_extract_all(format(terms[length(terms)]), "(?<=\\|).+?(?=\\))")[[1]]
    re_strip = gsub(" ", "", re_term, fixed = TRUE)
    re_data = df[,match(re_strip, colnames(df))]
    Z = model.matrix(~ -1 + factor(re_data))

    IV = all.names(terms[-c(1,2,length(terms))])
    IV_data = df[,match(IV, colnames(df))]
    if(factor_X == TRUE){
      X = model.matrix(~ factor(IV_data))
    }else{
      X = model.matrix(~ IV_data)
    }


    N = nrow(X)
    I = diag(1, N)
    q = ncol(as.matrix(Y))
  }

  return(list(X = X, Z = Z, Y = Y, N =N, I =I, q=q))
}







