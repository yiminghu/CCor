
source('./concentration.R')
source('./CCor.mCCor.R')
require("miscTools")

##examples shown in Section 3 of the paper

#######################two-module example#########################
two_module <- function(w1,w2,n1=100,n2=100,p=100,repl=200){
  ##an illustrating example of two modules
  ##Arguments
  ##w1: within-module correlation
  ##w2: between-module correlation
  ##n1: number of genes in first module
  ##n2: number of genes in second module
  ##p: sample size
  ##repl: number of replicates for simulation (for estimation of the concentrating value of CCor)
  ##Values
  ##W: a matrix indicating the correlation within and between modules
  ##W.ccen: a matrix of the concentrating values of CCor
  ##Wccen.est: a matrix of the average of CCor from simulation
  ##Wpcor.est: a matrix of the average of PCor from simulation
  N <- c(n1,n2)
  NN <- c(0,cumsum(N))
  W <- symMatrix(c(w1,w2,w1),2,byrow=T,upper=T)
  #W <- symMatrix(c(0.5,0.4,0.4,0.2,0.3,0.45,0.3,0.2,0.3,0.55,0.2,0.2,0.6,0.25,0.5),5,byrow=T,upper=T)
  A <- covmat(W,N)
  #A <- A + matrix(rnorm(nrow(A)*ncol(A),mean=0,sd=0.001),nrow(A),ncol(A)
  sig <- 1
  W.ccen <- W
  for(i in 1:nrow(W)){
    for(j in i:ncol(W)){
      if(i == j){
        W.ccen[i,j] <- concentration(sig=A,p=p,i=NN[i]+1,j=NN[i]+2)
      }else{
        W.ccen[i,j] <- concentration(sig=A,p=p,i=NN[i]+1,j=NN[j]+1)
        W.ccen[j,i] <- W.ccen[i,j] 
      }
    }
  }
  Wij <- matrix(0,length(N)*(1+length(N))/2,2)
  k <- 1
  for(i in 1:length(N)){
    for(j in i:length(N)){
      if(i == j){
        Wij[k,] <- c(NN[i]+1,NN[i]+2)
      }else{
        Wij[k,] <- c(NN[i]+1,NN[j]+1)
      }
      k <- k + 1
    }
  }
  Wc.est <- matrix(0,repl,length(N)*(1+length(N))/2)
  Wp.est <- matrix(0,repl,length(N)*(1+length(N))/2)
  for(i in 1:repl){
    Sp <- gen(sig^2*A,p)
    Ccr <- CCor(Sp)
    Pcr <- PCor(t(Sp))
    for(j in 1:ncol(Wc.est)){
      Wc.est[i,j] <- Ccr[Wij[j,1],Wij[j,2]]
      Wp.est[i,j] <- Pcr[Wij[j,1],Wij[j,2]]
    }
  }
  Wccen.est <- symMatrix(apply(Wc.est,2,mean),length(N),byrow=T,upper=T)
  Wpcor.est <- symMatrix(apply(Wp.est,2,mean),length(N),byrow=T,upper=T)
  return(list(W=W,W.ccen=W.ccen,Wccen.est=Wccen.est,Wpcor.est=Wpcor.est))
}









#######################five-module example#########################
##number of genes in each module
N <- c(50,10,70,70,60) 
NN <- c(0,cumsum(N))
##the within and between module correlation of these 5 modules
W <- symMatrix(c(0.5,0.35,0.35,0.2,0.3,0.45,0.3,0.2,0.3,0.55,0.3,0.2,0.6,0.25,0.5),5,byrow=T,upper=T) 
A <- covmat(W,N)
p = 100

####concentration and simulation of CCor and mCCor####
###calculating concentrating values of CCor###
W.ccen <- W
for(i in 1:nrow(W)){
  for(j in i:ncol(W)){
    if(i == j){
      W.ccen[i,j] <- concentration(sig=A,p=p,i=NN[i]+1,j=NN[i]+2)
    }else{
      W.ccen[i,j] <- concentration(sig=A,p=p,i=NN[i]+1,j=NN[j]+1)
      W.ccen[j,i] <- W.ccen[i,j] 
    }
  }
}
W.ccen

###calculating mean of CCor and mCCor using simulation
Wij <- matrix(0,length(N)*(1+length(N))/2,2)
k <- 1
for(i in 1:length(N)){
  for(j in i:length(N)){
    if(i == j){
      Wij[k,] <- c(NN[i]+1,NN[i]+2)
    }else{
      Wij[k,] <- c(NN[i]+1,NN[j]+1)
    }
    k <- k + 1
  }
}

repl <- 200
Wc.est <- matrix(0,repl,length(N)*(1+length(N))/2)
Wm.est <- matrix(0,repl,length(N)*(1+length(N))/2)


for(i in 1:repl){
  Sp <- gen(sig^2*A,p)
  Ccr <- CCor(Sp)
  mCcr <- mCCor(Sp,t=0.4)
  for(j in 1:ncol(Wc.est)){
    Wc.est[i,j] <- Ccr[Wij[j,1],Wij[j,2]]
    Wp.est[i,j] <- mCcr[Wij[j,1],Wij[j,2]]
  }
}

Wccor.est <- symMatrix(apply(Wc.est,2,mean),length(N),byrow=T,upper=T)
Wmccor.est <- symMatrix(apply(Wp.est,2,mean),length(N),byrow=T,upper=T)

W
round(Wccor.est,2)
round(Wmccor.est,2)


