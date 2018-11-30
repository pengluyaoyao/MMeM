library(MASS)
library(Matrix)
library(psych)


Z<-matrix(c(1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1),9,3)
ZtZ<-Z%*%t(Z)
X<-cbind(rep(1,9),rep(c(1,2,3),3))
B<-matrix(c(2,3,2.5,3.5),2,2)
T<-matrix(c(1,0.5,0.5,1),2,2)
E<-matrix(c(1,0.3,0.3,1),2,2)
I<-diag(1,9)
N<-nrow(I)

n<-10000
E_result<-matrix(0,n,4)
D_result<-matrix(0,n,4)
B_result<-matrix(0,n,4)
E_anova<-matrix(0,n,4)
D_anova<-matrix(0,n,4)
for (i in 1:n){
  U<-mvrnorm(n=3,mu=c(0,0), Sigma=T)
  E<-mvrnorm(n=9,mu=c(0,0), Sigma=E)

  Y<-X%*%B+Z%*%U+E
  HX<-X%*%ginv(t(X)%*%X)%*%t(X)
  XZ<- cbind(X,Z[,-3])
  HXZ<-XZ%*%ginv(t(XZ)%*%XZ)%*%t(XZ)
  Pz<-Z%*%ginv(t(Z)%*%Z)%*%t(Z)
  Qz<-I-Pz
  QzY<-Qz%*%(X%*%B+Z%*%U+E)
  PzY<-Pz%*%(X%*%B+Z%*%U+E)
  QzX<-Qz%*%X
  PzX<-Pz%*%X
  QzH<-QzX%*%ginv(t(QzX)%*%QzX)%*%t(QzX)
  PzH<-PzX%*%ginv(t(PzX)%*%PzX)%*%t(PzX)
  SSE_Q<-t(QzY)%*%QzY-t(QzY)%*%QzH%*%QzY
  SSE_P<-t(PzY)%*%PzY-t(PzY)%*%PzH%*%PzY

  Lambda_E<-SSE_Q/(tr(Qz)-rankMatrix(QzX)[1])
  Lambda_U<-SSE_P/(tr(Pz)-rankMatrix(PzX)[1])
  D_tilda<-(Lambda_U-Lambda_E)/3
  B_tilda<-ginv(t(PzX)%*%PzX)%*%t(PzX)%*%PzY
  E_result[i,]<-c(Lambda_E)
  D_result[i,]<-c(D_tilda)
  B_result[i,]<-c(B_tilda)

  ###HendersonIII
  SSE<-t(Y)%*%Y-t(Y)%*%HXZ%*%Y
  E_aov<-SSE/sum(diag(I-HXZ))
  E_anova[i,]<-E_aov
  SSD<-t(Y)%*%HXZ%*%Y-t(Y)%*%HX%*%Y
  D_aov<-(SSD-E_aov*sum(diag(HXZ-HX)))/tr(Z%*%t(Z)%*%(I-HX))
  D_anova[i,]<-D_aov
}

rbind(colMeans(E_result),apply(E_result,2, sd))
rbind(colMeans(D_result),apply(D_result,2, sd))
#colMeans(B_result)
rbind(colMeans(E_anova),apply(E_anova,2, sd))
rbind(colMeans(D_anova),apply(D_anova,2, sd))

b<-sum(diag(Qz-QzH))
a<-sum(diag(Pz-PzH))
t<-3
c<-tr(I-HXZ)
d<-tr(HXZ-HX)
e<-tr(Z%*%t(Z)%*%(I-HX))
A<-(b*Pz%*%(I-PzH)%*%Pz-a*Qz%*%(I-QzH)%*%Qz)/(a*b*t)
B<-(HXZ-HX-d*(I-HXZ)/c)/e
