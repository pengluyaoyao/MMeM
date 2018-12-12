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
find_Q <- function(T,E){
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
#####main function of multivariate mixed effects model_REML#####

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

  ###ZHZ
  H = I-X%*%ginv(t(X)%*%X)%*%t(X)
  h = ncol(X)
  s = ncol(Z)
  ZH = t(Z)%*%H
  ZHZ = t(Z)%*%H%*%Z

  for(iter in 1:maxit){
    print(iter)

    Q_results = find_Q(T, E)
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
  # rownames(Vcov) = c('T:y1&y1', 'E: y1&y1', 'T:y1&y2','E:y1&y2', 'T:y2&y2', 'E:y2&y2')
  # colnames(Vcov) = c('T:y1&y1', 'E: y1&y1', 'T:y1&y2','E:y1&y2', 'T:y2&y2', 'E:y2&y2')
  return(list(T.estimates = T[upper.tri(T,diag=TRUE)], E.estimates = E[upper.tri(E, diag = TRUE)], VCOV = Vcov))
}

############## TEST ON BIVARIATE CASE SIMULATION

z1 = c(rep(1,45))
z2 = c(rep(1,37))
z3 = c(rep(1,40))
z4 = c(rep(1,62))
z5 = c(rep(1,46))
z6 = c(rep(1,53))
z7 = c(rep(1,50))
z8 = c(rep(1,45))
z9 = c(rep(1,49))
z10= c(rep(1,47))
z11 = c(rep(1,47))
z12 = c(rep(1,53))
x1 = rbind(matrix(rep(c(1,1,0),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,1),18), 18,3, byrow = TRUE),matrix(rep(c(1,0,0),13), 13,3, byrow = TRUE))
x2 = rbind(matrix(rep(c(1,1,0),11), 11,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),13), 13,3, byrow = TRUE))
x3 = rbind(matrix(rep(c(1,1,0),11), 11,3, byrow = TRUE),matrix(rep(c(1,0,1),15), 15,3, byrow = TRUE),matrix(rep(c(1,0,0),14), 14,3, byrow = TRUE))
x4 = rbind(matrix(rep(c(1,1,0),21), 21,3, byrow = TRUE),matrix(rep(c(1,0,1),21), 21,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
x5 = rbind(matrix(rep(c(1,1,0),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,1),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,0),17), 17,3, byrow = TRUE))
x6 = rbind(matrix(rep(c(1,1,0),25), 25,3, byrow = TRUE),matrix(rep(c(1,0,1),19), 19,3, byrow = TRUE),matrix(rep(c(1,0,0),9), 9,3, byrow = TRUE))
x7 = rbind(matrix(rep(c(1,1,0),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,1),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,0),21), 21,3, byrow = TRUE))
x8 = rbind(matrix(rep(c(1,1,0),12), 12,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
x9 = rbind(matrix(rep(c(1,1,0),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,1),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,0),15), 15,3, byrow = TRUE))
x10 = rbind(matrix(rep(c(1,1,0),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,1),19), 19,3, byrow = TRUE),matrix(rep(c(1,0,0),11), 11,3, byrow = TRUE))
x11 = rbind(matrix(rep(c(1,1,0),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
x12 = rbind(matrix(rep(c(1,1,0),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,1),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,0),23), 23,3, byrow = TRUE))
#ind = cumsum(c(1,14,18,13,11,13,13,11,15,14,21,21,20,13,16,17,25,19,9,13,16,21,12,13,20,17,17,15,17,19,11,14,13,20,16,14,23))
N = 574
Z = as.matrix(bdiag(z1, z2, z3,z4,z5,z6,z7,z8,z9,z10,z11,z12))
I = diag(1, N)
X = rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
each = c(45,37,40,62,46,53,50,45,49,47,47,53)
Z_vec = c()
for(i in 1:12){
  Z_vec = c(Z_vec, rep(i, each[i]))
}

X_vec = X[,2]+X[,3]*2

B<-matrix(c(2, 3, 1, 2.5, 3.5, 1.5),3,2)
T.true<-matrix(c(100,5,5,10),2,2)
E.true<-matrix(c(60,13,13,40),2,2)

n = 1000
s = ncol(Z)
h = ncol(X)

#initial E and T
T.start = matrix(c(10,5,5,15),2,2)
E.start = matrix(c(6,1,1,3),2,2)

T_ = matrix(0, n,3)
E_ = matrix(0, n,3)
VC = list()
for (i in 1:n){
  U<-mvrnorm(n = s, mu=c(0,0), Sigma=T.true)
  e<-mvrnorm(n = N, mu=c(0,0), Sigma=E.true)
  Y<-X%*%B+Z%*%U+e

  data = as.data.frame(cbind(Y, X_vec, Z_vec))

  T.start = matrix(c(10,5,5,15),2,2)
  E.start = matrix(c(6,1,1,3),2,2)
  results = MMeM_reml(c(V1,V2) ~ X_vec + (1|Z_vec), data, factor_X = TRUE, T.start, T.start)

  T_[i,] = results$T.estimates
  E_[i,] = results$E.estimates
  VC[[i]] = results$VCOV
}

colMeans(T_)
colMeans(E_)
apply(T_, 2, sd)
apply(E_, 2, sd)

l = matrix(0,6,6)
for(i in 1:n){l = l+VC[[i]]}
sqrt(diag(l/n))


# Q <- matrix(c(0.001608,0.000975,-0.039382,0.023869),2,2)
# Q%*%E%*%t(Q)
# Q%*%T%*%t(Q)

##### TEST ON UNIVARIATE CASE REAL DATA == UNIVARIATE REML

# y1 = as.matrix(getME(mod1, 'y'))
# X = as.matrix(getME(mod1, 'X'))
# Z = as.matrix(getME(mod1, 'Z'))
# N = 246
# I = diag(1, N)
# T.start = 3
# E.start = 4
# q = 1


library(nlme)
library(lme4)
library(msm)
alcohol1 <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")
attach(alcohol1)
model.c <- lme(alcuse ~ coa, data=alcohol1, random= ~ 1 | id)
mod1<-lmer(alcuse ~ coa  +(1|id) ,alcohol1,REML=1)
var <-model.c$apVar

results = MMeM_reml(alcuse ~ coa + (1|id), alcohol1, T.start, T.start)

model.matrix( c(alcuse, age_14) ~ coa +(1|id), data=alcohol1)

