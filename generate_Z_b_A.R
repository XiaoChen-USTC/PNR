get_Z <- function(N,p){    #covariates
  Sigma=matrix(0, nrow=p, ncol=p)
  for(i in 1:p)
    for(j in 1:p)
      Sigma[i,j] = 0.5^(abs(i-j))
  return(mvrnorm(N,rep(0,p),Sigma))
}

block_N <- function(N,K){  #grouping
  b=NULL
  for(i in 1:K)
    if(i!=K) {b[((i-1)*round(N/K)+1):(i*round(N/K))]=i} else
      b[((i-1)*round(N/K)+1):N]=i
    return(b)
}

get_A <- function(N,b,P){        #Stochastic Block Model
  A=matrix(0, nrow=N, ncol=N)
  m=NULL
  for(i in 1:(N-1))
    for(j in (i+1):N){
      m=runif(1,0,1)
      if(m<(P[b[i],b[j]]))
        A[i,j]=A[j,i]=1 
    }
  for(i in 1:N)
    if(sum(A[i,])==0) A[i,(i+1)%%N]=A[(i+1)%%N,i]=1
  return(A)
}