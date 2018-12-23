############## TEST ON BIVARIATE CASE SIMULATION

# z1 = c(rep(1,45))
# z2 = c(rep(1,37))
# z3 = c(rep(1,40))
# z4 = c(rep(1,62))
# z5 = c(rep(1,46))
# z6 = c(rep(1,53))
# z7 = c(rep(1,50))
# z8 = c(rep(1,45))
# z9 = c(rep(1,49))
# z10= c(rep(1,47))
# z11 = c(rep(1,47))
# z12 = c(rep(1,53))
# x1 = rbind(matrix(rep(c(1,1,0),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,1),18), 18,3, byrow = TRUE),matrix(rep(c(1,0,0),13), 13,3, byrow = TRUE))
# x2 = rbind(matrix(rep(c(1,1,0),11), 11,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),13), 13,3, byrow = TRUE))
# x3 = rbind(matrix(rep(c(1,1,0),11), 11,3, byrow = TRUE),matrix(rep(c(1,0,1),15), 15,3, byrow = TRUE),matrix(rep(c(1,0,0),14), 14,3, byrow = TRUE))
# x4 = rbind(matrix(rep(c(1,1,0),21), 21,3, byrow = TRUE),matrix(rep(c(1,0,1),21), 21,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
# x5 = rbind(matrix(rep(c(1,1,0),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,1),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,0),17), 17,3, byrow = TRUE))
# x6 = rbind(matrix(rep(c(1,1,0),25), 25,3, byrow = TRUE),matrix(rep(c(1,0,1),19), 19,3, byrow = TRUE),matrix(rep(c(1,0,0),9), 9,3, byrow = TRUE))
# x7 = rbind(matrix(rep(c(1,1,0),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,1),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,0),21), 21,3, byrow = TRUE))
# x8 = rbind(matrix(rep(c(1,1,0),12), 12,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
# x9 = rbind(matrix(rep(c(1,1,0),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,1),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,0),15), 15,3, byrow = TRUE))
# x10 = rbind(matrix(rep(c(1,1,0),17), 17,3, byrow = TRUE),matrix(rep(c(1,0,1),19), 19,3, byrow = TRUE),matrix(rep(c(1,0,0),11), 11,3, byrow = TRUE))
# x11 = rbind(matrix(rep(c(1,1,0),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,1),13), 13,3, byrow = TRUE),matrix(rep(c(1,0,0),20), 20,3, byrow = TRUE))
# x12 = rbind(matrix(rep(c(1,1,0),16), 16,3, byrow = TRUE),matrix(rep(c(1,0,1),14), 14,3, byrow = TRUE),matrix(rep(c(1,0,0),23), 23,3, byrow = TRUE))
#
# Z = as.matrix(Matrix::bdiag(z1, z2, z3,z4,z5,z6,z7,z8,z9,z10,z11,z12))
# N = nrow(Z)
# X = rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
# each = c(45,37,40,62,46,53,50,45,49,47,47,53)
# Z_vec = c()
# for(i in 1:12){
#   Z_vec = c(Z_vec, rep(i, each[i]))
# }
#
# X_vec = X[,2]+X[,3]*2
#
# B<-matrix(c(2, 3, 1, 2.5, 3.5, 1.5),3,2)
# T.true<-matrix(c(100,5,5,10),2,2)
# E.true<-matrix(c(60,13,13,40),2,2)
#
# s = ncol(Z)
# h = ncol(X)
#
# #initial E and T
# T.start = matrix(c(10,5,5,15),2,2)
# E.start = matrix(c(6,1,1,3),2,2)
#
#
# U<-MASS::mvrnorm(n = s, mu=c(0,0), Sigma=T.true)
# e<-MASS::mvrnorm(n = N, mu=c(0,0), Sigma=E.true)
# Y<-X%*%B+Z%*%U+e
#
# data = as.data.frame(cbind(Y, X_vec, Z_vec))
# write.csv(data, file = "/Users/pengluyao/Documents/R_Package_Dvlp/MMeM/data/simulated_data.csv")

# data = read.csv('/Users/pengluyao/Documents/R_Package_Dvlp/MMeM/data/simulated_data.csv')
# T.start = matrix(c(10,5,5,15),2,2)
# E.start = matrix(c(10,1,1,3),2,2)
# results_reml = MMeM_reml(c(V1,V2) ~ X_vec + (1|Z_vec), data, factor_X = TRUE, T.start, E.start, maxit = 10)
# results_henderson = MMeM_henderson3(c(V1,V2) ~ X_vec + (1|Z_vec), data, factor_X = TRUE)

#save(data, file="/Users/pengluyao/Documents/R_Package_Dvlp/MMeM/data/mydata.rda")
# # ##### TEST ON UNIVARIATE CASE REAL DATA == UNIVARIATE REML
#
# alcohol1 <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")
# attach(alcohol1)
# # model.c <- nlme::lme(alcuse ~ coa, data=alcohol1, random= ~ 1 | id)
# # mod1<-lme4::lmer(alcuse ~ age  +(1|id) ,alcohol1,REML=1)
# # summary(mod1)
# # stats::vcov(mod1, full =TRUE)
# # var <-model.c$apVar
#
# T.start = 3
# E.start = 4
# results = MMeM_reml(alcuse ~ age + (1|id), alcohol1, factor_X = FALSE, T.start, E.start)


