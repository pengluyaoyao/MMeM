# MMeM (Multivariate Mixed-effects Model)

# Description

Estimating the variance covariance components matrix under the multivariate mixed effects model. Currently this package supports multivariate mixed effects model with two response variables, one fixed effects and one random effects.

# Estimation Methods

- Multivariate REML method  
- Multivariate Henderson3 method

# Examples

- bivariate mixed effects model:

```
 library(MMeM)
 data(simdata)
 T.start = matrix(c(10,5,5,15),2,2)
 E.start = matrix(c(10,1,1,3),2,2)
 results_henderson = MMeM_henderson3(fml = c(V1,V2) ~ X_vec + (1|Z_vec), data = simdata, factor_X = TRUE)
 results_reml = MMeM_reml(fml = c(V1,V2) ~ X_vec + (1|Z_vec), data = simdata, factor_X = TRUE, T.start = T.start, E.start =      E.start, maxit = 10)
```
- univariate mixed effects model:
```
alcohol1 <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")
attach(alcohol1)
mod1<-lme4::lmer(alcuse ~ age  +(1|id) ,alcohol1,REML=1)
summary(mod1)

library(merDeriv)
vcov(mod1, full =TRUE)
# # var <-model.c$apVar
#
T.start = 3
E.start = 4
results = MMeM_reml(alcuse ~ age + (1|id), alcohol1, factor_X = FALSE, T.start, E.start)
```
# Values
MMeM_reml: 
- T.estimates: the estimated matrix of the variance covariance matrix of the block random effects 
- E.estimates is the estimated matrix of the variance covariance matrix of the residuals 
- VCOV is the asymptotic dispersion matrix of the estimated variance covariance components

MMeM_henderson3: 
- T.estimates: the estimated matrix of the variance covariance matrix of the block random effects with corresponding standard errors
- E.estimates is the estimated matrix of the variance covariance matrix of the residuals with corresponding standard errors


# References

Meyer, K. A. R. I. N. "Maximum likelihood estimation of variance components for a multivariate mixed model with equal design matrices." Biometrics 1985: 153-165

Wesolowska‐Janczarek, M. T. "Estimation of covariance matrices in unbalanced random and mixed multivariate models." Biometrical journal 26.6 (1984): 665-674.
