# simulation_group.R
# Simulation for Partial Network Recovery with Group Structure

# Get script directory from command args or use current directory
args <- commandArgs(trailingOnly = FALSE)
script_path <- grep("^--file=", args, value = TRUE)
if (length(script_path) > 0) {
  script_path <- sub("^--file=", "", script_path)
  sim_dir <- dirname(normalizePath(script_path))
} else {
  sim_dir <- getwd()
}
setwd(sim_dir)
code_dir <- file.path(sim_dir, "..")

# Source required R functions
source('../generate_Z_b_A.R')
source('../Y_theta.R')
source('../SMC1.R')
source('../SL1.R')
source('../main.R')
source('../get_y_lin.R')
source('../get_A_ER.R')

# Compile C++ functions
Rcpp::sourceCpp('../main.cpp')

# Load required packages
required_packages <- c("parallel", "devtools", "lars", "MASS", "glmnet", "ncvreg", "pbmcapply", "lpSolve")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(parallel)
library(devtools)
library(lars)
library(MASS)
library(glmnet)
library(ncvreg)
library(pbmcapply)
library(lpSolve)

# Simulation parameters
N <- 100
T <- 20
ran_num <- 100
q <- 1 / (3 * N^(0.45))
K <- 2
par <- matrix(c(0.2, -0.2, 2, 2), nrow = K)
P <- matrix(c(3 * q, q, q, 3 * q), 2, 2)


#' Simulation for Partial Network Recovery with Group Structure
#'
#' @param N Number of nodes
#' @param K Number of blocks
#' @param P Block probability matrix
#' @param T Number of time periods
#' @param par Parameter matrix: rows for groups, columns for alpha and beta
#' @param q Sparsity parameter
#' @return List containing accuracy, TPR, TNR, prediction error for different methods
simulation_group <- function(N,K,P,T,par,q){
  # Obtain group and adjacency matrix A
  b=block_N(N,K)
  Tr=TRUE
  while (Tr) {
    A=get_A(N,b,P)

    # Reorganize A
    SMC=SMC1(A,1)
    A=SMC$A
    m=SMC$m
    m2=SMC$m2
    k1=length(unique(b[m2]))
    b=c(b[m2],b[-m2])
    A12=A[1:m,-c(1:m)]
    S=SMC$S
    W=A/rowSums(A)
    n_t=rowSums(A)
    n_t=n_t[-c(1:m)]
    nr1=nr2=nr3=SMC$nr
    acc_SMC=SMC$acc
    TPR_SMC=SMC$TPR
    TNR_SMC=SMC$TNR

    M1=gen_y_lin_group(N,T+1,par,W,m,K,b)
    Y=M1$Y
    X=M1$X
    rhoEst=M1$rhoEst
    khat = M1$khat
    g=M1$g

    if(sum(rowSums(A[-c(1:m),-c(1:m)])==0)==0 && khat==K && k1==K){
      Tr=FALSE
    }
  }

  g_hat = g_hat_l = c(rep(1,N-m))
  L = L_l = c(rep(Inf,N-m))
  acc_e1=TPR_e1=TNR_e1=matrix(0,nrow=K,ncol=(N-m))
  acc_l1=TPR_l1=TNR_l1=matrix(0,nrow=K,ncol=(N-m))
  # Save A_e and A_l for each kk
  A_e_all = array(0, c(K, N-m, N))
  A_l_all = array(0, c(K, N-m, N))
  for(kk in 1:K){
    L_kk = L_kk_l = NULL
    for(k in 1:3){
      acc_e=TPR_e=TNR_e=NULL
      acc_l=TPR_l=TNR_l=NULL
      A_e=matrix(0,N-m,N)
      A_l=matrix(0,N-m,N)
      A_e[,c(1:m)]=t(A12)
      A_l[,c(1:m)]=t(A12)

      Y_e=Y[(m+1):N,-(T+1)]-rhoEst[kk,2]*(1/nr2)*t(A12)%*%(X[1:m,-(T+1)])-rhoEst[kk,1]
      Y_l=Y[(m+1):N,-(T+1)]-rhoEst[kk,2]*(1/nr3)*t(A12)%*%(X[1:m,-(T+1)])-rhoEst[kk,1]
      for(i in 1:(N-m)){

        # Estimation
        X_e=t(rhoEst[kk,2]*(1/nr2[i])*(X[,-(T+1)]))
        X_e=X_e[,-c(1:m,i+m)]
        SL_e=SL(X_e,Y_e,i,N,m,A)
        L_kk[i] = SL_e$L

        # Lasso
        X_l=t(rhoEst[kk,2]*(1/nr3[i])*(X[,-(T+1)]))
        X_l=X_l[,-c(1:m,i+m)]
        SL_l=SL0(X_l,Y_l,i,N,m,A)
        L_kk_l[i] = SL_l$L_l

        A_e[i,c((m+1):N)]=SL_e$beta
        acc_e[i]=SL_e$acc
        TPR_e[i]=SL_e$TPR
        TNR_e[i]=SL_e$TNR

        A_l[i,c((m+1):N)]=SL_l$beta
        acc_l[i]=SL_l$acc
        TPR_l[i]=SL_l$TPR
        TNR_l[i]=SL_l$TNR
      }
      acc_e1[kk,]=acc_e
      TPR_e1[kk,]=TPR_e
      TNR_e1[kk,]=TNR_e
      acc_l1[kk,]=acc_l
      TPR_l1[kk,]=TPR_l
      TNR_l1[kk,]=TNR_l
      nr2=rowSums(A_e)
      nr3=rowSums(A_l)
      nr2[nr2==0]=q*(N-m)
      nr3[nr3==0]=q*(N-m)

      # Save current kk's A_e and A_l
      A_e_all[kk,,] = A_e
      A_l_all[kk,,] = A_l
    }
    g_hat[L_kk<L]=kk
    g_hat_l[L_kk_l<L_l]=kk
    L[L_kk<L]=L_kk[L_kk<L]
    L_l[L_kk_l<L_l]=L_kk_l[L_kk_l<L_l]
  }

  # Update A_e and A_l based on final g_hat and g_hat_l
  for(i in 1:(N-m)){
    A_e[i,] = A_e_all[g_hat[i], i,]
    A_l[i,] = A_l_all[g_hat_l[i], i,]
  }

  acc_e = TPR_e = TNR_e = 0
  acc_l = TPR_l = TNR_l = 0
  for(i in 1:(N-m)){
    acc_e = acc_e + acc_e1[g_hat[i],i]/(N-m)
    TPR_e = TPR_e + TPR_e1[g_hat[i],i]/(N-m)
    TNR_e = TNR_e + TNR_e1[g_hat[i],i]/(N-m)
    acc_l = acc_l + acc_l1[g_hat[i],i]/(N-m)
    TPR_l = TPR_l + TPR_l1[g_hat[i],i]/(N-m)
    TNR_l = TNR_l + TNR_l1[g_hat[i],i]/(N-m)
  }

  g_hat1 = c(g, g_hat)
  g_hat_l1 = c(g,g_hat_l)
  gpn_e=matrix(0,nrow=K,ncol=K)
  for(j in 1:K){
    for(k in 1:K){
      gpn_e[j,k]=sum(g_hat1[b==k]==j)
    }
  }
  gpr_e = lp.assign(gpn_e,direction = "max")$objval/N

  gpn_l=matrix(0,nrow=K,ncol=K)
  for(j in 1:K){
    for(k in 1:K){
      gpn_l[j,k]=sum(g_hat_l1[b==k]==j)
    }
  }
  gpr_l = lp.assign(gpn_l,direction = "max")$objval/N

  rhoEst_e = diag(rhoEst[g_hat,2])
  alpha_e = rhoEst[g_hat,1]
  yhat_e=alpha_e + (1/nr2)*rhoEst_e%*%A_e%*%(X[,T+1])
  prederr_e=mean((Y[c((m+1):N),T+1]-yhat_e)^2)
  rhoEst_l = diag(rhoEst[g_hat_l,2])
  alpha_l = rhoEst[g_hat_l,1]
  yhat_l=alpha_l + (1/nr3)*rhoEst_l%*%A_l%*%(X[,T+1])
  prederr_l=mean((Y[c((m+1):N),T+1]-yhat_l)^2)

  return(list(acc_e=acc_e,TPR_e=TPR_e,TNR_e=TNR_e,
              acc_l=acc_l,TPR_l=TPR_l,TNR_l=TNR_l,
              prederr_e=prederr_e,prederr_l=prederr_l,
              gpr_e,gpr_l))
}

# Run simulation with parallel processing
cat("Simulation start time:", as.character(Sys.time()), "\n")

total <- ran_num

# Create progress log file
progress_log <- "progress_group.txt"
if (file.exists(progress_log)) file.remove(progress_log)

result <- pbmclapply(1:ran_num, function(l) {
  res <- simulation_group(N, K, P, T, par, q)
  cat("1", file = progress_log, append = TRUE)
  done_count <- file.info(progress_log)$size
  message(paste0("Done: index ", l, " | Completed: ", done_count, "/", ran_num, " | ", Sys.time()))
  return(res)
}, mc.cores = 20, mc.preschedule = FALSE)

# Remove progress log file
if (file.exists(progress_log)) file.remove(progress_log)

# Aggregate results
result1 <- 0
result2 <- matrix(0, ran_num, 10)
for (i in 1:ran_num) {
  result2[i, ] <- as.numeric(result[[i]])
}
result1 <- colMeans(result2, na.rm = TRUE)

write.table(result1,
  file = paste0("result_group_N=", N, "_T=", T, ".txt"),
  sep = "\t", row.names = FALSE, col.names = FALSE
)