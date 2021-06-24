New=function(Dat, tau,  maxIter=50, perc=0.5){
  clusterU=unique(Dat$clusterID) #unique cluster IDs
  n=length(clusterU) #numberof clusters
  m=as.vector(table(Dat$clusterID))#sample cluster size
  Mhat=aggregate(Dat$d2, list(Dat$clusterID), sum)[,2]
  #initial values of beta_tau and sigma_v^2 using QR for LME models without weights
  R=lqmm(y ~ x, random = ~ 1, group = clusterID, tau = tau, data = Dat)
  beta_tau_old=coef(R)
  sigma_v2_old=R$scale
	cat("beta_tau= ", beta_tau_old, "sigma_v2= ", sigma_v2_old, "\n")
	esp_sigma=1
	esp_beta=1
	iterN=0
while( (esp_sigma>1e-3 | esp_beta>1e-3) & iterN < maxIter){
	iterN=iterN+1
      vhat=rep(NA,n) #random effects estimates
      V=rep(NA,n) #variance of vhat|v
      for (i in 1:n){
      	Dati=subset(Dat, clusterID==clusterU[i])#data in cluster i
      	Xi=cbind(1, Dati$x)
      	epsiloni=Dati$y-Xi%*%beta_tau_old
      	#1. update random effects v
      	vhat[i]=coef(rq(epsiloni~1, tau=tau, weights=Dati$d2))
     	 	#2. compute variance of vhat|v
      	#normal kernel density
      	epsilonbar=sum(Dati$d2*epsiloni)/Mhat[i]
      	sigmai2=sum(Dati$d2*(epsiloni-epsilonbar)^2)/(Mhat[i]-1)
      	hi=1.06*sqrt(sigmai2)*Mhat[i]^(-1/5)
      	fi=sum(Dati$d2*dnorm((vhat[i]-epsiloni)/hi))/(Mhat[i]*hi) #fi is a scalor
      	A1i=sum(fi*Dati$d2)
      	#the formula for A2i depends on sampling design
      	#the following formula assumes a stratified SRSWR
      	Dati$q=tau-((epsiloni-vhat[i])<0)
      	si2=aggregate(Dati$q, list(Dati$strata), var)[,2]#variance of qi by strata
      	Mihat=aggregate(Dati$d2, list(Dati$strata), sum)[,2]#Mihat by strata
      	mi=aggregate(rep(1,m[i]), list(Dati$strata), sum)[,2]#sample size by strata
      	A2i=sum(Mihat^2/mi*(1-mi/Mihat)*si2)
      	V[i]=max(A2i/(A1i^2),1e-3)
    }

    precision1=1/V
    precision2=1/sigma_v2_old
    precision=precision1+precision2
    mu=precision1/precision*vhat #mean of v|vhat by cluster
    sigma2=1/precision #variance of v|vhat by cluster

	K=200
	for(i in 1:n){
		Dati=subset(Dat, clusterID==clusterU[i])
		hvi=rnorm(K, mean=mu[i], sd=sqrt(sigma2[i]))
    		hvi_v=rep(hvi,rep(nrow(Dati),K))
		Dati_n=do.call(rbind, replicate(K, Dati, simplify=FALSE))
   		Dati_n$ny=Dati_n$y - hvi_v
		Dati_n$ny=jitter(Dati_n$ny)
		if(i==1){
			Dat_n=Dati_n
		}else{
			Dat_n=rbind(Dat_n,Dati_n)
		}
	}


  	Dat_n$d=Dat_n$d1*Dat_n$d2

	s_inx=sample(1:nrow(Dat_n), round(nrow(Dat_n)*perc))
	Dat_n_in=Dat_n[s_inx,]
	R=rq(ny ~ x, tau=tau, data=Dat_n_in, weights=d, method="fn")
  	beta_tau_i=R$coefficients

	esp_beta=sum(beta_tau_i-beta_tau_old)
	beta_tau_old=beta_tau_i

    	u2=mu^2+sigma2 #E(v^2|vhat)
    	d1_u=sqldf("select d1 from Dat group by clusterID")
    	sigma_v2_new=sum(d1_u*u2)/sum(d1_u) - (sum(d1_u*mu)/sum(d1_u))^2

	esp_sigma=sum(sigma_v2_new-sigma_v2_old)

	sigma_v2_old=sigma_v2_new

	cat("beta_tau= ", beta_tau_old, "sigma_v2= ", sigma_v2_old, "\n")
}
return(list(beta_tau=beta_tau_old, sigma=sigma_v2_old))

}
