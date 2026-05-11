# SSE180.R
# Real data analysis: Shanghai Stock Exchange 180 stocks

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
required_packages <- c("readxl", "devtools", "glmnet", "ncvreg", "lars", "ggplot2", "cowplot", "dplyr", "reshape2", "MASS")
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

# Install and load StructureMC from GitHub
if (!require("StructureMC", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  suppressWarnings(remotes::install_github("celehs/StructureMC", force = TRUE))
}
library(StructureMC)

# Source required R functions
source('../SMC1.R')
source('../SL1.R')
source('../beta_OLS.R')
source('../RRCN.R')
source('../main.R')

# Compile C++ functions
Rcpp::sourceCpp('../main.cpp')

# Load SSE180 data
data_dir <- file.path(sim_dir, "real data", "SSE180")
shareholds <- read_excel(file.path(data_dir, "上证180股东.xlsx"))
shareholds <- shareholds[-1, -c(1, 2)]
Y <- read_excel(file.path(data_dir, "上证180.xlsx"), col_names = FALSE)
stockname <- read_excel(file.path(data_dir, "名称.xlsx"))

stockname <- as.data.frame(stockname)
stockname <- rbind("浦发银行", stockname)
names(stockname) <- NULL

# Build adjacency matrix
n <- dim(shareholds)[1]
A <- matrix(0, n, n)

shareholds[is.na(shareholds)] <- "NA"

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    for (k in 1:10) {
      for (l in 1:10) {
        if (shareholds[i, k] != "NA" && shareholds[j, l] != "NA") {
          if (shareholds[i, k] == shareholds[j, l]) {
            A[i, j] <- A[i, j] + 1
            A[j, i] <- A[i, j]
          }
        }
      }
    }
  }
}
A1 <- A
A1[A >= 1] <- 1
A1[A <= 0] <- 0

# Remove stocks with no connections
l <- NULL
for (i in 1:n) {
  if (sum(A1[i, ]) == 0) {
    l <- c(l, i)
  }
}

Y1 <- t(Y[-(1:4), -1])
l <- c(l, which(apply(Y1, 1, function(row) any(is.na(row)))))
names(l) <- NULL
l <- unique(l)
l <- l[order(l)]

A2 <- A1[-l, -l]
n1 <- dim(A2)[1]
stockname <- as.vector(stockname[-l, ])

# Compute log return
Y1 <- matrix(as.numeric(Y1[-l, ]), nrow = n1)
Y1 <- log(Y1 / 100 + 1)

# Run SMC method
SMC <- SMC1(A2, 4)
A3 <- SMC$A
m <- SMC$m
m2 <- SMC$m2
A12 <- A3[1:m, -c(1:m)]
S <- SMC$S
n_t <- rowSums(A3)
n_t <- n_t[-c(1:m)]
Y11 <- Y1[m2, ]
Y12 <- Y1[-m2, ]
Y1 <- rbind(Y11, Y12)


#' Run PNR analysis on real data
#'
#' @param Y1 Price data matrix
#' @param A3 Adjacency matrix
#' @param A12 Partial adjacency matrix
#' @param m Number of known nodes
#' @param n1 Total number of nodes
#' @return Matrix of results
run_SSE180 <- function(Y1, A3, A12, m, n1) {
  result <- matrix(0, 10, 6)
  T <- 100
  for (t in 1:10) {
    Y2 <- Y1[, (t - 1) * 20 + 1:(T + (t - 1) * 20)]

    M1 <- beta_OLS(Y2, A3, S)
    rhoEst <- M1$rhoEst

    # Signal lasso
    acc_e <- TPR_e <- TNR_e <- 0
    A_e <- matrix(0, n1 - m, n1 - m)
    for (i in 1:(n1 - m)) {
      Y_e <- Y2[(m + 1):n1, ] - rhoEst[2] * t(A12) %*% (Y2[1:m, ]) - rhoEst[1]
      X_e <- t(rhoEst[2] * Y2[, ])
      X_e <- X_e[, -c(1:m, i + m)]
      SL_e <- SL(X_e, Y_e, i, n1, m, A3)
      A_e[i, ] <- SL_e$beta
      acc_e <- acc_e + SL_e$acc
      TPR_e <- TPR_e + SL_e$TPR
      TNR_e <- TNR_e + SL_e$TNR
    }
    acc_e <- acc_e / (n1 - m)
    TPR_e <- TPR_e / (n1 - m)
    TNR_e <- TNR_e / (n1 - m)

    # Lasso
    acc_l <- TPR_l <- TNR_l <- 0
    A_l <- matrix(0, n1 - m, n1 - m)
    for (i in 1:(n1 - m)) {
      Y_l <- Y2[(m + 1):n1, ] - rhoEst[2] * t(A12) %*% (Y2[1:m, ]) - rhoEst[1]
      X_l <- t(rhoEst[2] * Y2[, ])
      X_l <- X_l[, -c(1:m, i + m)]
      SL_l <- SL0(X_l, Y_l, i, n1, m, A3)
      A_l[i, ] <- SL_l$beta
      acc_l <- acc_l + SL_l$acc
      TPR_l <- TPR_l + SL_l$TPR
      TNR_l <- TNR_l + SL_l$TNR
    }
    acc_l <- acc_l / (n1 - m)
    TPR_l <- TPR_l / (n1 - m)
    TNR_l <- TNR_l / (n1 - m)

    result[t, ] <- c(acc_e, acc_l, NA, TPR_e, TPR_l, NA)
  }

  # RRCN method
  for (t in 1:10) {
    Y2 <- Y1[, (t - 1) * 20 + 1:(T + (t - 1) * 20)]

    acc_rrcn <- TPR_rrcn <- TNR_rrcn <- 0
    A_rrcn <- matrix(0, n1 - m, n1)
    for (i in 1:(n1 - m)) {
      result_rrcn <- RRCN(Y2, i, n1 - m, n1, m, A3)
      A_rrcn[i, ] <- result_rrcn$beta
      acc_rrcn <- acc_rrcn + result_rrcn$acc
      TPR_rrcn <- TPR_rrcn + result_rrcn$TPR
      TNR_rrcn <- TNR_rrcn + result_rrcn$TNR
    }
    acc_rrcn <- acc_rrcn / (n1 - m)
    TPR_rrcn <- TPR_rrcn / (n1 - m)
    TNR_rrcn <- TNR_rrcn / (n1 - m)

    result[t, 3] <- acc_rrcn
    result[t, 6] <- TPR_rrcn
  }

  return(result)
}

# Run analysis
result <- run_SSE180(Y1, A3, A12, m, n1)
write.table(result, file = "result.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
