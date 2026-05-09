gen_y_lin <- function(N,T,par,W,S){
  epsilon=matrix(rnorm(N*T,sd=0.4),ncol=T)
  al=matrix(par[1],N,T)
  X=matrix(rnorm(N*T),ncol=T)
  Y=matrix(0,N,T)
  for(i in 1:T){
    Y[,i] = al[,i] + par[2] * W %*% X[,i] + epsilon[,i]
  }
  
  Xt=cbind(rep(1, T),as.vector(S%*%W%*%X[,]))
  invXX = solve(crossprod(Xt))
  Yvec=as.vector(S%*%Y[,])
  rhoEst = as.numeric(invXX%*%colSums(Xt*Yvec)) 
  sigmaHat2 = mean((Yvec - Xt%*%rhoEst)^2)
  return(list(Y=Y,X=X,rhoEst=rhoEst,sigmaHat2=sigmaHat2))
}

gen_y_lin_group <- function(N,T,par,W,m,K,b){
  epsilon=matrix(rnorm(N*T,sd=0.4),ncol=T)
  # al=matrix(par[1],N,T)
  X=matrix(rnorm(N*T),ncol=T)
  Y=matrix(0,N,T)
  for(i in 1:T){
    for(k in 1:K){
      Y[b==k,i] = par[k,1] + par[k,2] * W[b==k,] %*% X[,i] + epsilon[b==k,i]
    }
  }
  
  al = matrix(0,nrow=m,ncol=2)
  for(i in 1:m){
    X_it=cbind(rep(1, T),as.vector(W[i,]%*%X[,]))
    invXX = solve(crossprod(X_it))
    al[i,] = as.numeric(invXX%*%t(X_it)%*%Y[i,]) 
  }
  result <- kmeans(al, centers=K)
  g = result$cluster

  al1 = matrix(0,nrow=K,ncol=2)
  for(k in 1:length(unique(result$cluster))){
    X_t=matrix(0,nrow = 2, ncol = 2)
    Y_w = c(0,0)
    for(i in 1:m){
      if(g[i]==k){
        X_it = cbind(rep(1, T),as.vector(W[i,]%*%X[,]))
        X_t=X_t + t(X_it)%*%X_it
        Y_w = Y_w + t(X_it)%*%Y[i,]
      }
    }
    invXX = solve(X_t)
    al1[k,] = as.numeric(invXX%*%Y_w) 
  }
  return(list(Y=Y,X=X,rhoEst=al1,khat=length(unique(result$cluster)),g=g))
}



