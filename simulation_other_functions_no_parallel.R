library(quantreg)
library(MASS)
library(splines)
library(gam)
library(sampling)
library(lqmm)
library(splitstackshape)
library(sqldf)

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

flines <- paste("Q:/simulation/result/Other_methods_tau", tau, ".txt", sep="");

bias_incpi=matrix(0, nrow=ns, ncol=3)
bias_betai=matrix(0, nrow=ns, ncol=3)
bias_sigmai=c()
for(i in 1:ns){
        f_in=paste(rpath, "data/", i, ".txt", sep="")
        Dat=read.table(f_in, header=T, sep=",")
	  res1=UW1(Dat, tau)
	  res2=UW2(Dat, tau)
	  res3=WT(Dat, tau)
        bias_incpi[i,1]=res1$coefficients[1]-true_incp
        bias_betai[i,1]=res1$coefficients[2]-beta[2]
        bias_incpi[i,2]=res2$coefficients[1]-true_incp
        bias_betai[i,2]=res2$coefficients[2]-beta[2]
        bias_incpi[i,3]=res3$coefficients[1]-true_incp
        bias_betai[i,3]=res3$coefficients[2]-beta[2]
        bias_sigmai[i]=res2$sigmav2 - sigma_v2

}


bias_est=data.frame(bias_incpi, bias_betai, bias_sigmai)
colnames(bias_est)=c("Intercep-UW1", "Intercep-UW2", "Intercep-WT", "Beta-UW1", "Beta-UW2", "Beta-WT", "Sigma_v2")
write.table(bias_est, flines, sep="\t", quote=F, col.names=T, row.names=F)

