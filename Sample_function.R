spFunction <- function (date_in) {
  return(tryCatch(stratified(date_in, group="strata", size=c("0"=size.in.c[1], "1"=size.in.c[2])), error=function(e) NULL))
}

Sample=function(N, n, beta, sigma_e2, sigma_v2, cluster.max, size.in.c){
  #this function generates one clustered stratified sample
  #N: Number of Clusters in the population
  #n: Number of clusters in the Sample
  #beta: beta0 and beta1
  #sigma_e2: variance of error
  #sigma_v2: variance of v_i
  #cluster.max: maximum population cluster size
  #size.in.c: common stratified sample size in each cluster (length 2 vector)
  
  #population
  v=rnorm(N, mean=0, sd=sqrt(sigma_v2)) #population cluster random effects
  top=exp(2.5+v)
  u=top/(1+top)
  M=round(cluster.max*u) #population cluster size
  #M=pmax(M,60)
  inx=M>0			##Remove any possible 0 sample cluster
  M=M[inx]
  clusterID=rep(1:length(M), M)
  person=unlist(sapply(M, function (x) 1:x))
  personID=paste(clusterID, person, sep='-')
  popsize=sum(M) #population size
  x=rnorm(popsize)
  e=rnorm(popsize, mean=0,sd=sqrt(sigma_e2))
  V=rep(v[inx],M)
  y=beta[1]+beta[2]*x+V+e
  strata=(e>0)+0
  data.pop=data.frame(clusterID, personID, strata,x, y)
  
  #sample
  #stage 1 (PPS)
  Pi=n*M/sum(M)#probability proportional to size (PPS)
  d1=1/Pi #stage 1 weights
  In=UPrandomsystematic(Pi,eps=1e-6) #cluster selection: PSU
  select=(In==1)
  clusterID.s=(1:N)[select]#selected cluster ID
  data.sp=subset(data.pop, clusterID %in% clusterID.s )#select clustered samples
  data.sp$d1=d1[data.sp$clusterID]# stage 1 weights
  #stage 2 (Stratified SRSWOR within cluster)
  data.sp1=split(data.sp, data.sp$clusterID) #split into a list by clusterID

  st=lapply(data.sp1, spFunction)
  #stage 2 data
  data.sp2=do.call(rbind, st)
  #stage 2 weights
  data.sp.in=data.sp[data.sp$clusterID %in% unique(data.sp2$clusterID),]
  size1=table(data.sp.in$clusterID, data.sp.in$strata)
  wts=as.vector(t(size1))/size.in.c
  size2=rep(size.in.c, length(unique(data.sp2$clusterID)))
  d2=rep(wts,size2)
  data.sp2$d2=d2
  data.sp2
}

