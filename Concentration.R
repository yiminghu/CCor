####functions for simulation####
##data generation
##covariance matrix
require('miscTools')
covmat = function(W,N){
  ##generating covariance matrix of block structure
  ##Arguments
  ##N: a vector indicating the number of elements in each module
  ##W: a matrix indicating the correlation within and between modules (assume correlations between elements from same 
  ## module are the same)
  ##Values
  ##A matrix with specified block structure with diagnal blocks corresponding to within-module and off-diagnal blocks corresponding to between-module covariance
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
  ##generating p samples with mean 0 and covariance matrix Sig
  ##Arguments
  ##Sig: covariance matrix
  ##p: sample size
  ##Values
  ##a p by n matrix, each row is a sample with n genes
  A = chol(Sig)
  A <- t(A)
  n = ncol(A)
  Sample = matrix(NA,p,n)
  for(i in 1:p){
    Sample[i,] = A%*%rnorm(n)
  }
  Sample
}


concentration <- function(sig,p,i,j){
  ##calculating the concentration value of CCor of a pair of genes
  ##Arguments
  ##sig: the covariance matrix of all the genes
  ##p: sample size
  ##i,j: the labels of a gene pair
  ##Values
  ##concentrating value of CCor between gene i and j
  n <- nrow(sig)
  all <- (1:n)[-c(i,j)]
  mu1 <- sig[i,all]
  mu2 <- sig[j,all]
  m <- n - 2
  #print('calculating covariance matrix of U and V...')
  S11 <- matrix(0,m,m)
  for(s in 1:m){
    for(t in s:m){
      S11[s,t] <- sig[all[s],all[t]] + sig[i,all[s]]*sig[i,all[t]]
      S11[t,s] <- S11[s,t]
    }
  }
  S22 <- matrix(0,m,m)
  for(s in 1:m){
    for(t in s:m){
      S22[s,t] <- sig[all[s],all[t]] + sig[j,all[s]]*sig[j,all[t]]
      S22[t,s] <- S22[s,t]
    }
  }
  S12 <- matrix(0,m,m)
  for(s in 1:m){
    for(t in 1:m){
      S12[s,t] <- sig[i,j]*sig[all[s],all[t]] + sig[i,all[t]]*sig[j,all[s]]
    }
  }
  S <- rbind(cbind(S11,S12),cbind(t(S12),S22))
  A <- diag(rep(1,m))-rep(1,m)%*%t(rep(1,m))/m
  #print('calculating D...')
  decS <- eigen(S/p)
  S.root <- decS$vectors%*%sqrt(diag(decS$values))%*%t(decS$vectors)
  S.root.inv <- decS$vectors%*%sqrt(diag(1/decS$values))%*%t(decS$vectors)
  decB <- eigen(S.root%*%rbind(cbind(matrix(0,m,m),A),cbind(A,matrix(0,m,m)))%*%S.root)
  M <- t(decB$vectors)
  C1 <- M%*%S.root.inv%*%c(mu1,mu2)
  Lam1 <- decB$values
  ED <- sum(Lam1*(1+as.vector(C1)^2))/2
  #print('calculating H...')
  decS11 <- eigen(S11/p)
  S11.root <- decS11$vectors%*%sqrt(diag(decS11$values))%*%t(decS11$vectors)
  S11.root.inv <- decS11$vectors%*%sqrt(diag(1/decS11$values))%*%t(decS11$vectors)
  decB <- eigen(S11.root%*%A%*%S11.root)
  P <- t(decB$vectors)
  C2 <- P%*%S11.root.inv%*%mu1
  Lam2 <- decB$values
  EH <- sum(Lam2*(1+as.vector(C2)^2))
  #print('calculating G...')
  decS22 <- eigen(S22/p)
  S22.root <- decS22$vectors%*%sqrt(diag(decS22$values))%*%t(decS22$vectors)
  S22.root.inv <- decS22$vectors%*%sqrt(diag(1/decS22$values))%*%t(decS22$vectors)
  decB <- eigen(S22.root%*%A%*%S22.root)
  Q <- t(decB$vectors)
  C3 <- Q%*%S22.root.inv%*%mu2
  Lam3 <- decB$values
  EG <- sum(Lam3*(1+as.vector(C3)^2))
  ED/sqrt(EH*EG)
}
