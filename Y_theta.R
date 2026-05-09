gen_y <- function(N,T,par,W,S){
  epsilon=matrix(rnorm(N*T,sd=0.4),ncol=T)
  al=matrix(par[1],N,T)
  I=diag(rep(1,N))
  Q=solve(I-par[2]*W)
  Y=Q%*%(epsilon+al)
  
  I = S%*%rep(1,N)
  Xt=cbind(rep(I, T),as.vector(S%*%W%*%Y))
  invXX = solve(crossprod(Xt))
  Yvec=as.vector(S%*%Y)
  rhoEst = as.numeric(invXX%*%colSums(Xt*Yvec)) 
  sigmaHat2 = mean((Yvec - Xt%*%rhoEst)^2)
  return(list(Y=Y,rhoEst=rhoEst,sigmaHat2=sigmaHat2))
}