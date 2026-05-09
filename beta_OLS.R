beta_OLS <- function(Y,W,S){
  N=dim(Y)[1]
  T=dim(Y)[2]
  I = S%*%rep(1,N)
  Xt=cbind(rep(I, T),as.vector(S%*%W%*%Y))
  invXX = solve(crossprod(Xt))
  Yvec=as.vector(S%*%Y)
  rhoEst = as.numeric(invXX%*%colSums(Xt*Yvec)) 
  sigmaHat2 = mean((Yvec - Xt%*%rhoEst)^2)
  return(list(Y=Y,rhoEst=rhoEst,sigmaHat2=sigmaHat2))
}