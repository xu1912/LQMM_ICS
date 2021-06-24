library(quantreg)
library(MASS)
library(splines)
library(gam)
library(sampling)
library(lqmm)
library(splitstackshape)
library(sqldf)
library(doParallel)

rpath <- "Q:/simulation/";
setwd(rpath);
source("New_function.R")
source("Sample_function.R")
source("Comparison_functions.R")

perc=1
maxIter=50

tau=0.4
beta=c(0.5,1)
sigma_v2=1
true_incp=beta[1]+qnorm(tau)

ns=200

flines <- paste("Q:/simulation/result/Perc", perc, "_maxIter", maxIter, "_tau", tau, ".txt", sep="");

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

	  res1=UW1()
	  res2=UW2()
	  res3=WT1()
}


bias_est=data.frame(bias_incpi, bias_betai, bias_sigmai)
colnames(bias_est)=c("Intercep", "Beta", "Sigma_v2")
write.table(bias_est, flines, sep="\t", quote=F, col.names=T, row.names=F)

