get_A_ER_d <- function(N){        #dense ER Model
  A=matrix(0, nrow=N, ncol=N)
  m=matrix(0, nrow=N, ncol=N)
  for(i in 2:N)
    for(j in 1:(i-1)){
      m[i,j]=runif(1,0,1)
      if(m[i,j]<(1-(10/N+1/(3*N^(0.5)))))
        A[i,j]=A[j,i]=1 
    }
  for(i in 1:N)
    if(sum(A[i,]==0)) A[i,(i+1)%%N]=A[(i+1)%%N,i]=1
  return(A)
}