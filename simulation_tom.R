require('miscTools')
covmat = function(W,N){
  ##N: a vector indicating the number of elements in each cluster
  ##W: a matrix indicating the correlation within and between cluster (assume correlations between elements from same 
  ## cluster are the same)
  
  n <- length(N)
  pos <- c(0,cumsum(N))
  Sig <- matrix(0,sum(N),sum(N))
  ##diagnal blocks## 
  for(i in 1:n){
    Sig[(pos[i]+1):pos[i+1],(pos[i]+1):pos[i+1]] <- W[i,i]
  }
  ##off diagnal blocks
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      Sig[(pos[i]+1):pos[i+1],(pos[j]+1):pos[j+1]] <- W[i,j]
      Sig[(pos[j]+1):pos[j+1],(pos[i]+1):pos[i+1]] <- W[j,i]
    }
  }
  diag(Sig) <- 1
  Sig
}



gen = function(Sig,p){
  A = chol(Sig)
  A <- t(A)
  n = ncol(A)
  Sample = matrix(NA,p,n)
  for(i in 1:p){
    Sample[i,] = A%*%rnorm(n)
  }
  Sample
}

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
library(WGCNA)
library(ROCR)
repl <- 100 ## number of replicates
N <- c(50,30,70,90,60)
NN <- cumsum(N)
W <- symMatrix(c(0.5,0.3,0.3,0.2,0.3,0.6,0.3,0.2,0.3,0.55,0.3,0.2,0.55,0.25,0.5),5,byrow=T,upper=T) 
A <- covmat(W,N)
auc4 <- matrix(0,repl,length(nsample)) ## auc of Tom
level <- 0
k <- 1
pb <- txtProgressBar(0, repl*length(nsample), style=3)
for(i in 1:length(nsample)){
  for(j in 1:repl){
    Sp <- gen(A,nsample[i]) + matrix(rnorm(nsample[i]*sum(N),mean=0,sd=level),nsample[i],sum(N))
    tom1 <- TOMsimilarityFromExpr(Sp)
    Tomw <- within(tom1,NN); Tomb <- between(tom1,NN)    
    pred4 <- prediction(c(Tomw,Tomb),c(rep(1,length(Tomw)),rep(0,length(Tomb))))
    auc4[j,i] <- unlist(performance(pred4,'auc')@y.values)
    k <- k + 1
    setTxtProgressBar(pb, k)
  }
}


cbind(apply(auc4,2,mean),apply(auc4,2,sd))
