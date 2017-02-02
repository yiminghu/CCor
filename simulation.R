######## simulation to compare CCor, PCor and Gaussian graphical model ########
library(ROCR)
library(GeneNet)
library(space)
#library(WGCNA)
### simulating five-module case with different sample size ###
### calculate mean and sd of FP FN and AUC of all replicates ###
##########################################################
### functions for extracting within and between pairs ###
##########################################################
within <- function(M,NN){
  NN <- c(0,NN)
  res <- -1
  for(i in 1:(length(NN)-1)){
    res <- c(res,M[(NN[i]+1):NN[i+1],(NN[i]+1):NN[i+1]][lower.tri(M[(NN[i]+1):NN[i+1],(NN[i]+1):NN[i+1]])])
  }
  res[-1]
}
between <- function(M,NN){
  NN <- c(0,NN)
  res <- -1
  for(i in 1:(length(NN)-2)){
    res <- c(res,as.numeric(M[(NN[i]+1):NN[i+1],(NN[i+1]+1):NN[length(NN)]]))
  }
  res[-1]
}
##########################################################

repl <- 100 ## number of replicates
# N <- c(60,20,70,10,10,20,50,60)
N <- c(50,30,70,90,60)
NN <- cumsum(N)
W <- symMatrix(c(0.5,0.3,0.3,0.2,0.3,0.6,0.3,0.2,0.3,0.55,0.3,0.2,0.55,0.25,0.5),5,byrow=T,upper=T) 
A <- covmat(W,N)

## check what's the precision matrix like in this case ##
pA <- solve(A)
npA <- -diag(1/sqrt(diag(pA)))%*%pA%*%diag(1/sqrt(diag(pA)))
pAw <- within(npA,NN)
pAb <- between(npA,NN)
hist(pAw)
hist(pAb)
#########################################################

nsample <- sum(N)/c(30,20,10,5,2,1) ## possible sample size
auc1 <- matrix(0,repl,length(nsample)) ## auc of CCor
auc2 <- matrix(0,repl,length(nsample)) ## auc of PCor
auc3 <- matrix(0,repl,length(nsample)) ## auc of PPCor1, penalized GeneNet
level <- 0
k <- 1
pb <- txtProgressBar(0, repl*length(nsample), style=3)
for(i in 1:length(nsample)){
  for(j in 1:repl){
    Sp <- gen(A,nsample[i]) + matrix(rnorm(nsample[i]*sum(N),mean=0,sd=level),nsample[i],sum(N))
    Ccr <- CCor(t(Sp))
    Pcr <- cor(Sp)
    PPcr1 <- ggm.estimate.pcor(Sp,method='static')
    Cw <- within(Ccr,NN); Cb <- between(Ccr,NN)
    Pw <- within(Pcr,NN); Pb <- between(Pcr,NN)
    PPw1 <- within(PPcr1,NN); PPb1 <- between(PPcr1,NN)
    pred1 <- prediction(c(Cw,Cb),c(rep(1,length(Cw)),rep(0,length(Cb))))
    auc1[j,i] <- unlist(performance(pred1,'auc')@y.values)
    pred2 <- prediction(c(Pw,Pb),c(rep(1,length(Pw)),rep(0,length(Pb))))
    auc2[j,i] <- unlist(performance(pred2,'auc')@y.values)
    pred3 <- prediction(c(PPw1,PPb1),c(rep(1,length(PPw1)),rep(0,length(PPb1))))
    auc3[j,i] <- unlist(performance(pred3,'auc')@y.values)
    k <- k + 1
    setTxtProgressBar(pb, k)
  }
}

cbind(apply(auc1,2,mean),apply(auc1,2,sd))
cbind(apply(auc2,2,mean),apply(auc2,2,sd))
cbind(apply(auc3,2,mean),apply(auc3,2,sd))



repl <- 100 ## number of replicates
N <- c(60,25,70,15,20,35,75,90,50,60)
NN <- cumsum(N)
thd1 <- -0.15
thd2 <- 0.3
cor_str <- c(runif(1,0.5,1),runif(9,thd1,thd2),runif(1,0.5,1),runif(8,thd1,thd2),runif(1,0.5,1),runif(7,thd1,thd2),runif(1,0.5,1),runif(6,thd1,thd2),runif(1,0.5,1),runif(5,thd1,thd2),runif(1,0.5,1),runif(4,thd1,thd2),runif(1,0.5,1),runif(3,thd1,thd2),runif(1,0.5,1),runif(2,thd1,thd2),runif(1,0.5,1),runif(1,thd1,thd2),runif(1,0.5,1))
W <- symMatrix(cor_str,10,byrow=T,upper=T)
W
A <- covmat(W,N)
tst <- gen(A,10)
nsample <- c(10,25,50,100,200,300)
auc1 <- matrix(0,repl,length(nsample)) ## auc of CCor
auc2 <- matrix(0,repl,length(nsample)) ## auc of PCor
auc3 <- matrix(0,repl,length(nsample)) ## auc of PPCor1, penalized GeneNet
level <- 2
k <- 1
pb <- txtProgressBar(0, repl*length(nsample), style=3)
for(i in 1:length(nsample)){
  for(j in 1:repl){
    Sp <- gen(A,nsample[i]) + matrix(rnorm(nsample[i]*sum(N),mean=0,sd=level),nsample[i],sum(N))
    Ccr <- CCor(t(Sp))
    Pcr <- cor(Sp)
    PPcr1 <- ggm.estimate.pcor(Sp,method='static')
    Cw <- within(Ccr,NN); Cb <- between(Ccr,NN)
    Pw <- within(Pcr,NN); Pb <- between(Pcr,NN)
    PPw1 <- within(PPcr1,NN); PPb1 <- between(PPcr1,NN)
    pred1 <- prediction(c(Cw,Cb),c(rep(1,length(Cw)),rep(0,length(Cb))))
    auc1[j,i] <- unlist(performance(pred1,'auc')@y.values)
    pred2 <- prediction(abs(c(Pw,Pb)),c(rep(1,length(Pw)),rep(0,length(Pb))))
    auc2[j,i] <- unlist(performance(pred2,'auc')@y.values)
    pred3 <- prediction(c(PPw1,PPb1),c(rep(1,length(PPw1)),rep(0,length(PPb1))))
    auc3[j,i] <- unlist(performance(pred3,'auc')@y.values)
    k <- k + 1
    setTxtProgressBar(pb, k)
  }
}

cbind(apply(auc1,2,mean),apply(auc1,2,sd))
cbind(apply(auc2,2,mean),apply(auc2,2,sd))
cbind(apply(auc3,2,mean),apply(auc3,2,sd))

