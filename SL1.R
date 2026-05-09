source('main.R')
Rcpp::sourceCpp('main.cpp')

SL <- function(X,Y,i,N,m,A){
  cvfit0=cv.glmnet(x=X,y=Y[i,],nfolds=5,intercept=TRUE,type.measure = "mse")
  fit.l=glmnet(x=X,y=Y[i,],lambda =cvfit0$lambda.1se,alpha = 1,intercept = TRUE)
  beta0=coef.glmnet(fit.l)[-1]
  
  fit=CV.SIGNAL(Y=Y[i,],X=X,beta0 = beta0,nfolds=5,nlambda1=cvfit0$lambda[1:10],alpha=c(seq(0.1,0.9,0.1),1/seq(0.1,0.9,0.1)), 
                weights = 1,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.s=SIGNAL_l(Y=Y[i,],X=X,beta0 = beta0,lambda1=fit$lambda.1se[1],
                 lambda2 =fit$lambda.1se[2],weights=1,iter_max = 2000,delta=1e-7)
  beta1=fit.s$Beta
  beta1[beta1<0.3]=0
  beta1[beta1>0.3]=1
  L = mean((Y[i,]-X%*%beta1)^2) 
  if(i==1) {beta1=c(0,beta1)} else 
    if(i==(N-m)) {beta1=c(beta1,0)} else
      beta1=c(beta1[1:(i-1)],0,beta1[i:(N-m-1)])
  M=beta1+A[i+m,-c(1:m)]
  acc=sum(beta1==A[i+m,-c(1:m)])/(N-m)
  TPR=sum(M==2)/sum(A[i+m,-c(1:m)]==1)
  TNR=sum(M==0)/sum(A[i+m,-c(1:m)]==0)
  return(list(beta=beta1,acc=acc,TPR=TPR,TNR=TNR,L=L))
}

SL0<- function(X,Y,i,N,m,A){
  cvfit0=cv.glmnet(x=X,y=Y[i,],nfolds=5,intercept=TRUE,type.measure = "mse")
  fit.l=glmnet(x=X,y=Y[i,],lambda =cvfit0$lambda.1se,alpha = 1,intercept = TRUE)
  beta0=coef.glmnet(fit.l)[-1]
  beta0[beta0>0.3]=1
  beta0[beta0<0.3]=0
  L_l = mean((Y[i,]-X%*%beta0)^2) 
  if(i==1) {beta0=c(0,beta0)} else 
    if(i==(N-m)) {beta0=c(beta0,0)} else
      beta0=c(beta0[1:(i-1)],0,beta0[i:(N-m-1)])
  M=beta0+A[i+m,-c(1:m)]
  acc=sum(beta0==A[i+m,-c(1:m)])/(N-m)
  TPR=sum(M==2)/sum(A[i+m,-c(1:m)]==1)
  TNR=sum(M==0)/sum(A[i+m,-c(1:m)]==0)
  return(list(beta=beta0,acc=acc,TPR=TPR,TNR=TNR,L_l=L_l))
}



