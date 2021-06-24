# unweighted quantile regression without considering random effects
UW1=function(Dat, tau){
  R=rq (y ~ x, tau=tau, data=Dat)
  list(coefficients=coef(R))
}
# unweighted quantile regression considering random effects
UW2=function(Dat, tau){
  R=lqmm(y ~ x, random = ~ 1, group = clusterID, tau = tau, data = Dat)
  list(coefficients=coef(R), sigmav2=R$scale)
}

# weighted quantile regression without considering random effects
WT=function(Dat, tau){
  d=Dat$d1*Dat$d2
  R=rq (y ~ x, tau=tau, data=Dat, weights=d)
  list(coefficients=coef(R))
}