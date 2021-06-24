library(quantreg)
library(MASS)
library(splines)
library(gam)
library(sampling)
library(lqmm)
library(splitstackshape)
library(sqldf)
library(doParallel)

nc=2			##number of cpus/threads can be used for parallel computing
if(nc>1){
	noParallel=FALSE
}else{
	noParallel=TRUE
}

rpath <- "Q:/simulation/";
setwd(rpath);
source("New_function.R")
source("Sample_function.R")

perc=1
maxIter=50


tau=0.4
beta=c(0.5,1)
sigma_v2=1
true_incp=beta[1]+qnorm(tau)

ns=10

flines <- paste("Q:/simulation/result/Perc", perc, "_maxIter", maxIter, "_tau", tau, ".txt", sep="");

if(noParallel==T){
bias_incpi=c()
bias_betai=c()
bias_sigmai=c()
for(i in 1:ns){
        f_in=paste(rpath, "data/", i, ".txt", sep="")
        Dat=read.table(f_in, header=T, sep=",")
        res=New(Dat, tau, perc=perc, maxIter=maxIter)

        bias_incpi[i]=res$beta_tau[1]-true_incp
        bias_betai[i]=res$beta_tau[2]-beta[2]
        bias_sigmai[i]=res$sigma - sigma_v2
}
}

cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)
bias_est_b=foreach(i = 1:ns, .combine=c, .export=c(), .packages=c('quantreg', 'sqldf', 'lqmm', 'sampling', 'splitstackshape')) %dopar% {
        f_in=paste(rpath, "data/", i, ".txt", sep="")
        Dat=read.table(f_in, header=T, sep=",")
        res=New(Dat, tau, perc=perc, maxIter=maxIter)

        bias_incp=res$beta_tau[1]-true_incp
        bias_beta=res$beta_tau[2]-beta[2]
        bias_sigma=res$sigma - sigma_v2
        return(data.frame(bias_incp, bias_beta, bias_sigma))

}
parallel::stopCluster(cl)

bias_est=matrix(unlist(bias_est_b), ncol = 3, byrow = TRUE)
colnames(bias_est)=c("Intercep", "Beta", "Sigma_v2")
write.table(bias_est, flines, sep="\t", quote=F, col.names=T, row.names=F)

