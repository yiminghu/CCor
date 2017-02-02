CCor <- function(X){
  ##calculating CCor
  ##Arguments
  ##X: p by n matrix
  ##Values
  ##n by n matrix of CCor
  begin = Sys.time()
  dyn.load("./ccor.so")
  X <- t(X)
  Pcor <- cor(X)  
  n <- nrow(Pcor)
  n = as.integer(n)
  Ans <- matrix(0,n,n)
  storage.mode(Pcor) <- 'double'
  storage.mode(Ans) <- 'double'
  res = .C("ccor",Pcor,Ans,n,DUP=F)
  print(Sys.time()-begin)
  ##print time for calculating CCor
  return(Ans)
}


mCCor <- function(X, threshold=0.3, w1){
  ##calculating mCCor
  ##Arguments
  ##X: p by n matrix
  ##threshold: threshold for mCCor
  ##w1: weight for mCCor
  ##Values
  ##n by n matrix of mCCor
  begin = Sys.time()
  dyn.load("./mccor.so")
  X <- t(X)
  Pcor <- cor(X)
  n <- nrow(Pcor)
  n = as.integer(n)
  threshold = as.double(threshold)
  Ans <- matrix(0,n,n)
  w <- rep(0,n)
  storage.mode(Pcor) <- 'double'
  storage.mode(Ans) <- 'double'
  storage.mode(w) <- 'double'
  storage.mode(w1) <- 'double'
  res = .C("mccor",Pcor,Ans,n,threshold,w,w1,DUP=F)
  print(Sys.time()-begin)
  return(Ans)
}

