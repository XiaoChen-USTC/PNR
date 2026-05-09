# Install and load required packages if not already available

# Check and install MASS if needed
if (!requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
library(MASS)

# Check and install matrixcalc if needed
if (!requireNamespace("matrixcalc", quietly = TRUE)) {
  install.packages("matrixcalc")
}
library(matrixcalc)

# Install and load StructureMC from GitHub
if (!require("StructureMC", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  suppressWarnings(remotes::install_github("celehs/StructureMC", force = TRUE))
}
library(StructureMC)

SMC1 <- function(A,i){
  #partial matrix
  m=c(i,which(A[i,]==1))
  N=ncol(A)
  m1=m2=length(m)
  A11=A[m,m]
  A12=A[m,-m]
  A22_t=A[-m,-m]
  A=rbind(cbind(A11,A12),cbind(t(A12),A22_t))
  S=diag(A[i,])
  Arecovery = rbind(cbind(A11,A12),cbind(t(A12),matrix(NA,nrow=N-m1,ncol=N-m2)))
  A22 = smc.FUN(Arecovery, 2, "True", m1, m2) 
  A22_e=A22
  A22_e[A22_e>0.5]=1
  A22_e[A22_e<0.5]=0
  acc_SMC=sum(colSums(A22_t==A22_e))/(N-m1)^2
  M=A22_t+A22_e
  TPR_SMC=sum(colSums(M==2))/sum(colSums(A22_t==1))
  TNR_SMC=sum(colSums(M==0))/sum(colSums(A22_t==0))
  Ahat=rbind(cbind(A11,A12),cbind(t(A12),A22))
  nr=rowSums(Ahat)
  nr=nr[-c(1:m1)]
  nr[nr==0]=mean(nr)
  return(list(acc=acc_SMC,TPR=TPR_SMC,TNR=TNR_SMC,S=S,A=A,m=m1,nr=nr,m2=m))
}