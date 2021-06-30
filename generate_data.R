library(MASS)
library(splines)
library(gam)
library(sampling)
library(splitstackshape)

set.seed(998)

sigma_v2=14/3
sigma_e2=1
beta=c(0.5, 1)
N=10000
n=100
cluster.max=500
size.in.c=30/3*(2:1)

rpath <- paste("/work/cxu2/ChenXu/QR/data_sigmav2_", round(sigma_v2), sep="");
if(!dir.exists(rpath)){
        dir.create(rpath, showWarnings = FALSE)
}
setwd (rpath);
source("../Sample_function.R")

ns=1000

for(i in 1:ns){

        f_lines=paste(i, ".txt", sep="")
        if(file.exists(f_lines)){
                next
        }
        Dat=Sample(N, n, beta, sigma_e2, sigma_v2, cluster.max, size.in.c)
        write.table(Dat, f_lines,quote=F, row.names=F, col.names=T, sep=",")

}
