# urban_crime.R
# Real data analysis: Urban crime data

Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

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

# Load required packages
required_packages <- c("readxl", "devtools", "glmnet", "ncvreg", "lars", "ggplot2", "cowplot", "dplyr", "reshape2", "MASS", "tidyr", "lpSolve", "StructureMC")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(readxl)
library(devtools)
library(glmnet)
library(ncvreg)
library(lars)
library(ggplot2)
library(cowplot)
library(dplyr)
library(reshape2)
library(MASS)
library(tidyr)
library(lpSolve)
library(StructureMC)

# Source required R functions
source('../SL1.R')

# Load urban crime data
data_dir <- file.path(sim_dir, "real data", "urban crime")
par1 <- read.table(file.path(data_dir, "part1.txt"))[, -1]
par1 <- as.matrix(par1)
par2 <- read.table(file.path(data_dir, "part2.txt"))[, -1]
par2 <- as.matrix(par2)

# Wn
Ws <- read.table(file.path(data_dir, "W.txt"))[-1, -1]
Ws <- as.matrix(Ws)
As <- Ws
As[As != 0] <- 1

# x0
socio <- read.table(file.path(data_dir, "socio.txt"))
socio <- as.matrix(socio)

tract <- socio[, 1]

# log(1+c) transformed counts of Part I and Part II
ln_par1 <- log(par1 + 1)
ln_par2 <- log(par2 + 1)
df1 <- data.frame(time = c(1:72), par1 = as.vector(colMeans(ln_par1)), par2 = as.vector(colMeans(ln_par2)))

# Plot average crime trends over time
df1 <- data.frame(time = c(1:72), par1 = as.vector(colMeans(ln_par1)), par2 = as.vector(colMeans(ln_par2)))
p1 <- ggplot(df1, aes(x = time)) +
  geom_line(aes(y = par1, color = "Part I"), size = 1) +
  geom_line(aes(y = par2, color = "Part II"), size = 1) +
  labs(x = "Month T",
       y = "Averaged logarithmic transformed crimes",
       color = "Type") +
  theme_minimal()

# Plot by census
df2 <- data.frame(Census = c(1:138), par1 = as.vector(rowMeans(ln_par1)), par2 = as.vector(rowMeans(ln_par2)))
data_long <- tidyr::pivot_longer(df2, cols = c("par1", "par2"), names_to = "Type", values_to = "Value")

p2 <- ggplot(data_long, aes(x = Census, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Census", y = "Averaged logarithmic transformed crimes") +
  theme_minimal()

# Save plots to file (avoid graphics device issues)
cowplot::save_plot("crime_trends.png", p1)
cowplot::save_plot("crime_by_census.png", p2)


N = as.numeric(dim(par1)[1])
T = as.numeric(dim(par1)[2])

m1=c(1,which(As[1,]==1))
m=length(m1)
A11=As[m1,m1]
A12=As[m1,-m1]
A21=As[-m1,m1]
A22_t=As[-m1,-m1]
A=rbind(cbind(A11,A12),cbind(A21,A22_t))
W=A/rowSums(A)
n_t = nr2 = nr3 = as.vector(rowSums(A))[-(1:m)]

part1 = rbind(ln_par1[m1,], ln_par1[-m1,])
part2 = rbind(ln_par2[m1,], ln_par2[-m1,])

b_i = matrix(0,nrow=m,ncol=3)
for(i in 1:m){
  X_it=cbind(rep(1, T-1),as.vector(W[i,]%*%part1[,-T]),as.vector(W[i,]%*%part2[,-1]))
  invXX = solve(crossprod(X_it))
  b_i[i,] = as.numeric(invXX%*%t(X_it)%*%part1[i,-1]) 
}

# Elbow plot for k selection
wss <- function(k) {
  kmeans(b_i, k, nstart = 10)$tot.withinss
}

k.values <- 1:5
wss_values <- sapply(k.values, wss)

df <- data.frame(
  k = k.values,
  wss = wss_values
)

p3 <- ggplot(data = df, aes(x = k, y = wss)) +
  geom_line(color = "#1E90FF", size = 0.5) +
  geom_point(color = "gold", size = 2) +
  labs(x = "Community number K", y = "Total Within-Cluster Sum of Squares") +
  theme_minimal()

cowplot::save_plot("elbow_plot.png", p3)

result = kmeans(b_i,centers=3)
g = result$cluster

rhoEst = matrix(0,nrow=3,ncol=3)
for(k in 1:3){
  X_t=matrix(0,nrow = 3, ncol = 3)
  Y_w = c(0,0,0)
  for(i in 1:m){
    if(g[i]==k){
      X_it = cbind(rep(1, T-1),as.vector(W[i,]%*%part1[,-T]), as.vector(W[i,]%*%part2[,-1]))
      X_t=X_t + t(X_it)%*%X_it
      Y_w = Y_w + as.vector(t(X_it)%*%part1[i,-1])
    }
  }
  invXX = solve(X_t)
  rhoEst[k,] = as.numeric(invXX%*%Y_w) 
}

g_hat = g_hat_l = c(rep(1,N-m))
L = L_l = c(rep(Inf,N-m))
acc_e1=TPR_e1=TNR_e1=matrix(0,nrow=3,ncol=(N-m))
acc_l1=TPR_l1=TNR_l1=matrix(0,nrow=3,ncol=(N-m))
for(kk in 1:3){
  print(kk)
  L_kk = L_kk_l = NULL
  for(k in 1:1){
    acc_e=TPR_e=TNR_e=NULL
    acc_l=TPR_l=TNR_l=NULL
    A_e=matrix(0,N-m,N)
    A_l=matrix(0,N-m,N)
    A_e[,c(1:m)]=t(A12)
    A_l[,c(1:m)]=t(A12)
    
    Y_e=part1[(m+1):N,-1]-rhoEst[kk,2]*(1/nr2)*t(A12)%*%(part1[1:m,-T]) -rhoEst[kk,3]*(1/nr2)*t(A12)%*%(part2[1:m,-1])-rhoEst[kk,1]
    Y_l=part1[(m+1):N,-1]-rhoEst[kk,2]*(1/nr3)*t(A12)%*%(part1[1:m,-T]) -rhoEst[kk,3]*(1/nr3)*t(A12)%*%(part2[1:m,-1])-rhoEst[kk,1]
    for(i in 1:(N-m)){
      # print(i)
      
      #estimation
      X_e=t(rhoEst[kk,2]*(1/nr2[i])*(part1[,-T]) + rhoEst[kk,3]*(1/nr2[i])*(part2[,-1]))
      X_e=X_e[,-c(1:m,i+m)]
      SL_e=SL(X_e,Y_e,i,N,m,A)
      L_kk[i] = SL_e$L
      
      #lasso
      X_l=t(rhoEst[kk,2]*(1/nr3[i])*(part1[,-T]) + rhoEst[kk,3]*(1/nr3[i])*(part2[,-1]))
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
  }
  g_hat[L_kk<L]=kk
  g_hat_l[L_kk_l<L_l]=kk
  L[L_kk<L]=L_kk[L_kk<L]
  L_l[L_kk_l<L_l]=L_kk_l[L_kk_l<L_l]
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

# Save clustering results
g_hat1 <- c(g, g_hat)

g_hat2 <- cbind(tract, g_hat1)
write.csv(g_hat2, file = "group.csv", row.names = FALSE)

b_i1 = matrix(0,nrow=N,ncol=3)
for(i in 1:N){
  X_it=cbind(rep(1, T-1),as.vector(W[i,]%*%part1[,-T]),as.vector(W[i,]%*%part2[,-1]))
  invXX = solve(crossprod(X_it))
  b_i1[i,] = as.numeric(invXX%*%t(X_it)%*%part1[i,-1]) 
}

result1 <- kmeans(b_i1, centers=3)
g1 = result1$cluster
centers = result1$centers

print(sum(g1 == 1))
print(sum(g1 == 2))
print(sum(g1 == 3))
print(centers)






