library(MASS)
library(Matrix)
library(psych)


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

Z = as.matrix(bdiag(z1, z2, z3,z4,z5,z6,z7,z8,z9,z10,z11,z12))

X = rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)

B<-matrix(c(2, 3, 1, 2.5, 3.5, 1.5),3,2)
T.true<-matrix(c(100,5,5,10),2,2)
E.true<-matrix(c(60,13,13,40),2,2)

N = nrow(X)
s = ncol(Z)
h = ncol(X)
I = diag(1, N)

n<-1000
# E_result<-matrix(0,n,4)
# D_result<-matrix(0,n,4)
# B_result<-matrix(0,n,4)
E_anova<-matrix(0,n,4)
T_anova<-matrix(0,n,4)
sd_E<-matrix(0,n,4)
sd_T<-matrix(0,n,4)

HX<-X%*%ginv(t(X)%*%X)%*%t(X)
XZ<- cbind(X,Z)
HXZ<-XZ%*%ginv(t(XZ)%*%XZ)%*%t(XZ)
C_E = (I - HXZ)
C_ET = HXZ - HX
C_T = C_ET-C_E*tr(C_ET)/tr(C_E)

for (i in 1:n){
  U<-mvrnorm(n = s, mu=c(0,0), Sigma=T.true)
  e<-mvrnorm(n = N, mu=c(0,0), Sigma=E.true)
  Y<-X%*%B+Z%*%U+e

  # Pz<-Z%*%ginv(t(Z)%*%Z)%*%t(Z)
  # Qz<-I-Pz
  # QzY<-Qz%*%(X%*%B+Z%*%U+E)
  # PzY<-Pz%*%(X%*%B+Z%*%U+E)
  # QzX<-Qz%*%X
  # PzX<-Pz%*%X
  # QzH<-QzX%*%ginv(t(QzX)%*%QzX)%*%t(QzX)
  # PzH<-PzX%*%ginv(t(PzX)%*%PzX)%*%t(PzX)
  # SSE_Q<-t(QzY)%*%QzY-t(QzY)%*%QzH%*%QzY
  # SSE_P<-t(PzY)%*%PzY-t(PzY)%*%PzH%*%PzY
  #
  # Lambda_E<-SSE_Q/(tr(Qz)-rankMatrix(QzX)[1])
  # Lambda_U<-SSE_P/(tr(Pz)-rankMatrix(PzX)[1])
  # D_tilda<-(Lambda_U-Lambda_E)/3
  # B_tilda<-ginv(t(PzX)%*%PzX)%*%t(PzX)%*%PzY
  # E_result[i,]<-c(Lambda_E)
  # D_result[i,]<-c(D_tilda)
  # B_result[i,]<-c(B_tilda)

  ###HendersonIII
  #SSE<-t(Y)%*%Y-t(Y)%*%HXZ%*%Y
  E_aov<-(t(Y)%*%C_E%*%Y)/tr(C_E)
  E_anova[i,]<-E_aov
  cov_E = (2/tr(C_E))*(E_aov%x%E_aov)
  sd_E[i,]<- sqrt(diag(cov_E))


  #SSD<-t(Y)%*%HXZ%*%Y-t(Y)%*%HX%*%Y
  #T_aov<-(SSD-E_aov*tr(HXZ-HX))/tr(Z%*%t(Z)%*%(I-HX))
  T_aov = (t(Y)%*%C_T%*%Y)/tr(Z%*%t(Z)%*%(I-HX))
  T_anova[i,]<-T_aov
  ET_aov = (t(Y)%*%C_ET%*%Y)
  cov_ET = (2/tr(C_ET))*(ET_aov%x%ET_aov)

  cov_T<- (cov_ET - cov_E*(tr(C_ET)^2))/(tr(Z%*%t(Z)%*%(I-HX))^2)
  sd_T[i,] = sqrt(diag(cov_T))
}


rbind(colMeans(E_anova),apply(E_anova,2, sd))
rbind(colMeans(T_anova),apply(T_anova,2, sd))
colMeans(sd_E)
colMeans(sd_T)
