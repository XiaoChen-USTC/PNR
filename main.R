
####the correlation function######
Tpm=function(p,rho){
  A=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){A[i,j]=1}
      else{A[i,j]=rho^(abs(i-j))}
    }
    
  }
  return(A)
}
########The signal lasso Function lambda1*|X|+lambda2*|X-1|########
SIGNAL_l=function(Y,X,beta0,m=NULL,lambda1,lambda2,weights=1,constant=TRUE,iter_max,delta,
                  XX=NULL, XY=NULL, X_mean=NULL, Y_mean=NULL){
  p=ncol(X)
  N=nrow(X)
  lambda3=weights*lambda2
  ptm=proc.time()
  
  if(is.null(m)){
    # Compute pre-calculated matrices if not provided
    if (is.null(XX) || is.null(XY)) {
      if (constant) {
        X_mean = colMeans(X)
        Y_mean = mean(Y)
        X_c = scale(X, center=X_mean, scale=FALSE)
        Y_c = Y - Y_mean
      } else {
        X_mean = rep(0, p)
        Y_mean = 0
        X_c = X
        Y_c = Y
      }
      XX = crossprod(X_c)
      XY = crossprod(X_c, Y_c)
    }
    
    Re2 = diag(XX)
    ep1 = (lambda1 + lambda3) / Re2
    ep2 = (lambda1 - lambda3) / Re2
    
    # Call C++ covariance update function
    fit.c = Signal_c_cov(XY=XY, XX=XX, beta0=beta0, Re2=Re2, 
                         ep1=ep1, ep2=ep2, iter_max=iter_max, p=p, delta=delta)
    beta1 = fit.c$Beta
    mu1 = Y_mean - sum(X_mean * beta1)
    iters = fit.c$iters
    
  }else{
    n=N/m
    iter=1
    res=1
    while (iter<iter_max && res>delta) {
      Y1=t(matrix(Y-X%*%beta0,m,n))
      sig=cov(Y1)
      eig_sig=eigen(sig)
      sig_in=eig_sig$vectors%*%diag((eig_sig$values)^(-1/2))%*%t(eig_sig$vectors)
      l12=replicate(n,sig_in, simplify=FALSE)
      sig_block=bdiag(l12) 
      Y_hat=as.vector(sig_block%*%Y)
      X_hat=as.matrix(sig_block%*%X)
      
      if (constant) {
        X_mean_hat = colMeans(X_hat)
        Y_mean_hat = mean(Y_hat)
        X_c_hat = scale(X_hat, center=X_mean_hat, scale=FALSE)
        Y_c_hat = Y_hat - Y_mean_hat
      } else {
        X_mean_hat = rep(0, p)
        Y_mean_hat = 0
        X_c_hat = X_hat
        Y_c_hat = Y_hat
      }
      
      XX_hat = crossprod(X_c_hat)
      XY_hat = crossprod(X_c_hat, Y_c_hat)
      Re2_hat = diag(XX_hat)
      
      ep1 = (lambda1 + lambda3) / Re2_hat
      ep2 = (lambda1 - lambda3) / Re2_hat
      
      fit.c = Signal_c_cov(XY=XY_hat, XX=XX_hat, beta0=beta0, Re2=Re2_hat, 
                           ep1=ep1, ep2=ep2, iter_max=iter_max, p=p, delta=delta)
      beta1 = fit.c$Beta
      mu1 = Y_mean_hat - sum(X_mean_hat * beta1)
      
      res=max(abs(beta1-beta0))
      beta0=beta1
      iter=iter+1
    }
    iters = iter
  }
  
  mse=mean((Y-mu1-X%*%beta1)^2)
  return(list(Mu=mu1, Beta=beta1, Mse=mse, Iter=iters, Times=proc.time()-ptm))
}

#######The CV function for signal lasso####
CV.SIGNAL=function(Y,X,beta0,nfolds,nlambda1,alpha=seq(0.3,0.7,length.out=5),weights=1,constant=TRUE,iter_max,delta){
  n=length(Y)
  folds=cv.folds(n,nfolds)
  
  XX_list = list(); XY_list = list()
  Y_mean_list = list(); X_mean_list = list()
  X1_list = list(); Y1_list = list()
  X2_list = list(); Y2_list = list()
  beta_list = list()
  
  # 1. Pre-calculate covariance matrices and center data before the grid search
  for (k in 1:nfolds) {
    train_idx = as.vector(unlist(folds[-k]))
    test_idx = as.vector(unlist(folds[k]))
    
    X1_train = X[train_idx, , drop=FALSE]
    Y1_train = Y[train_idx]
    
    if (constant) {
      X_mean_list[[k]] = colMeans(X1_train)
      Y_mean_list[[k]] = mean(Y1_train)
      X1_c = scale(X1_train, center=X_mean_list[[k]], scale=FALSE)
      Y1_c = Y1_train - Y_mean_list[[k]]
    } else {
      X_mean_list[[k]] = rep(0, ncol(X))
      Y_mean_list[[k]] = 0
      X1_c = X1_train
      Y1_c = Y1_train
    }
    
    XX_list[[k]] = crossprod(X1_c)          
    XY_list[[k]] = crossprod(X1_c, Y1_c)    
    
    X1_list[[k]] = X1_train
    Y1_list[[k]] = Y1_train
    X2_list[[k]] = X[test_idx, , drop=FALSE]
    Y2_list[[k]] = Y[test_idx]
    beta_list[[k]] = beta0
  }
  
  # 2. Smooth the parameter grid by sorting lambda1 in descending order
  param_grid = expand.grid(lambda2 = nlambda1, alpha = alpha)
  param_grid$lambda1 = param_grid$alpha * param_grid$lambda2
  param_grid = param_grid[order(param_grid$lambda1, decreasing = TRUE), ]
  
  re = rep(NA, nrow(param_grid)) # Initialize with NA
  lambda3_res = matrix(0, nrow(param_grid), 2)
  
  # Early Stopping Parameters
  min_mse = Inf
  stop_threshold = 1.5 # If current MSE > 1.5 * min_mse, we stop
  
  # 3. Main grid search with Early Stopping
  for (row in 1:nrow(param_grid)) {
    lambda1 = param_grid$lambda1[row]
    lambda2 = param_grid$lambda2[row]
    se = 0
    
    for (k in 1:nfolds) {
      fit = SIGNAL_l(Y=Y1_list[[k]], X=X1_list[[k]], beta0=beta_list[[k]], m=NULL,
                     lambda1=lambda1, lambda2=lambda2, weights=weights,
                     constant=constant, iter_max=iter_max, delta=delta,
                     XX=XX_list[[k]], XY=XY_list[[k]], 
                     X_mean=X_mean_list[[k]], Y_mean=Y_mean_list[[k]])
      
      beta_list[[k]] = fit$Beta
      se = se + mean((Y2_list[[k]] - fit$Mu - X2_list[[k]]%*%fit$Beta)^2)
    }
    
    current_mse = se / nfolds
    re[row] = current_mse
    lambda3_res[row, ] = c(lambda1, lambda2)
    
    # Update minimum MSE found so far
    if (current_mse < min_mse) {
      min_mse = current_mse
    }
    
    # --- Early Stopping Logic ---
    # If the error starts to explode (overfitting or numerical instability), 
    # break the loop to save time.
    if (current_mse > min_mse * stop_threshold) {
      # cat("Early stopping triggered at row", row, "due to MSE increase.\n")
      break
    }
  }
  
  # Filter out rows that weren't calculated due to early stopping
  valid_idx = !is.na(re)
  re = re[valid_idx]
  lambda3_res = lambda3_res[valid_idx, , drop=FALSE]
  
  index = which.min(re)
  return(list(lambda.1se=lambda3_res[index,], Re=re))   
}


###############The product penalty  function   lambda*|X|*|X-1|  ########
PNR_l=function(Y,X,beta0, m=NULL,lambda,constant=TRUE,iter_max,delta){
  #Y: the response vector
  #X: the covariate matrix
  #m:if NULL : error type is independent, else: m is the number of observations of each subject
  #beta0: initial value of the estimated parameter
  #lambda1: the first penalty parameter
  #lambda2: the second penalty parameter
  #iter_max: the maximal iterations
  #delta: control the accuracy
  p=ncol(X)
  N=nrow(X)
  n=N/m
  ptm=proc.time()
  if(is.null(m)){
    Re2=diag(t(X)%*%X) 
    con=1*constant
    fit.c=PNR_c(Y=Y, X=X, beta0=beta0, Re2=Re2,
                lambda=lambda, constant=con,iter_max=iter_max, p=p, delta=delta)
  }
  else{
    iter=1
    res=1
    while (iter<iter_max && res>delta) {
      
      Y1=t(matrix(Y-X%*%beta0,m,n))
      sig=cov(Y1)
      eig_sig=eigen(sig)
      sig_in=eig_sig$vectors%*%diag((eig_sig$values)^(-1/2))%*%t(eig_sig$vectors)
      l12=replicate(n,sig_in, simplify=FALSE)
      sig_block=bdiag(l12) 
      Y_hat=as.vector(sig_block%*%Y)
      X_hat=as.matrix(sig_block%*%X)
      Re2=diag(t(X_hat)%*%X_hat) 
      con=1*constant
      fit.c=PNR_c(Y=Y_hat, X=X_hat, beta0=beta0, Re2=Re2,
                    lambda=lambda, constant=con,iter_max=iter_max, p=p, delta=delta)
      beta1=fit.c$Beta
      res=max(abs(beta1-beta0))
      beta0=beta1
      iter=iter+1
    }
  }
  mse=mean((Y-fit.c$Mu-X%*%fit.c$Beta)^2)
  return(list(Mu=fit.c$Mu,Beta=fit.c$Beta, Mse=mse, Iter=fit.c$iters,Times=proc.time()-ptm))
}
