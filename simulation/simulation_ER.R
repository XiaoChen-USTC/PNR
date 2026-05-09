# simulation_ER.R
# Simulation for Partial Network Recovery on Erdős-Rényi random graphs

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
source('../Y_theta.R')
source('../SMC1.R')
source('../SL1.R')
source('../main.R')
source('../get_A_ER.R')

# Compile C++ functions
Rcpp::sourceCpp('../main.cpp')

# Load required packages
required_packages <- c("parallel", "devtools", "lars", "MASS", "glmnet", "ncvreg", "pbmcapply")
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

# Simulation parameters
N <- 100
T <- 50
ran_num <- 100
par <- c(0.2, 2)


#' Simulation for Partial Network Recovery on Erdős-Rényi random graphs
#'
#' @param N Number of nodes
#' @param T Number of time periods
#' @param par Vector of parameters: c(alpha, beta)
#' @return List containing accuracy, TPR, and TNR for different methods
simulation_ER <- function(N, T, par) {
  # Generate adjacency matrix A using ER model
  A <- get_A_ER(N)

  # Run SMC algorithm for partial network recovery
  SMC <- SMC1(A, 1)
  A <- SMC$A
  m <- SMC$m
  A12 <- A[1:m, -c(1:m)]
  S <- SMC$S
  W <- A / rowSums(A)
  n_t <- rowSums(A)
  n_t <- n_t[-c(1:m)]
  nr1 <- nr2 <- nr3 <- SMC$nr
  acc_SMC <- SMC$acc
  TPR_SMC <- SMC$TPR
  TNR_SMC <- SMC$TNR

  # Generate response data
  M1 <- gen_y(N, T, par, W, S)
  Y <- M1$Y
  rhoEst <- M1$rhoEst
  sigmaHat2 <- M1$sigmaHat2

  # Initialize results for oracle 1
  acc_o1 <- TPR_o1 <- TNR_o1 <- 0

  # Oracle method 1: using true parameters and degree
  for (i in 1:(N - m)) {
    Y_o1 <- Y[(m + 1):N,] - par[2] * (1 / n_t[i]) * t(A12) %*% (Y[1:m, ]) - par[1]
    X_o1 <- t(par[2] * (1 / n_t[i]) * (Y[, ]))
    X_o1 <- X_o1[, -c(1:m, i + m)]
    SL_o1 <- SL(X_o1, Y_o1, i, N, m, A)

    acc_o1 <- acc_o1 + SL_o1$acc
    TPR_o1 <- TPR_o1 + SL_o1$TPR
    TNR_o1 <- TNR_o1 + SL_o1$TNR
  }
  acc_o1 <- acc_o1 / (N - m)
  TPR_o1 <- TPR_o1 / (N - m)
  TNR_o1 <- TNR_o1 / (N - m)

  # Iterative estimation methods
  for (k in 1:6) {
    acc_o2 <- TPR_o2 <- TNR_o2 <- 0
    acc_e <- TPR_e <- TNR_e <- 0
    acc_l <- TPR_l <- TNR_l <- 0
    A_o2 <- matrix(0, N - m, N - m)
    A_e <- matrix(0, N - m, N - m)
    A_l <- matrix(0, N - m, N - m)

    for (i in 1:(N - m)) {
      # Oracle 2: using true parameters and estimated degree
      Y_o2 <- Y[(m + 1):N,] - par[2] * (1 / nr1[i]) * t(A12) %*% (Y[1:m, ]) - par[1]
      X_o2 <- t(par[2] * (1 / nr1[i]) * (Y[, ]))
      X_o2 <- X_o2[, -c(1:m, i + m)]
      SL_o2 <- SL(X_o2, Y_o2, i, N, m, A)

      # Estimation: using estimated parameters and degree
      Y_e <- Y[(m + 1):N,] - rhoEst[2] * (1 / nr2[i]) * t(A12) %*% (Y[1:m, ]) - rhoEst[1]
      X_e <- t(rhoEst[2] * (1 / nr2[i]) * (Y[, ]))
      X_e <- X_e[, -c(1:m, i + m)]
      SL_e <- SL(X_e, Y_e, i, N, m, A)

      # Lasso: using true parameters and estimated degree
      Y_l <- Y[(m + 1):N,] - par[2] * (1 / nr3[i]) * t(A12) %*% (Y[1:m, ]) - par[1]
      X_l <- t(par[2] * (1 / nr3[i]) * (Y[, ]))
      X_l <- X_l[, -c(1:m, i + m)]
      SL_l <- SL0(X_l, Y_l, i, N, m, A)

      A_o2[i, ] <- SL_o2$beta
      acc_o2 <- acc_o2 + SL_o2$acc
      TPR_o2 <- TPR_o2 + SL_o2$TPR
      TNR_o2 <- TNR_o2 + SL_o2$TNR

      A_e[i, ] <- SL_e$beta
      acc_e <- acc_e + SL_e$acc
      TPR_e <- TPR_e + SL_e$TPR
      TNR_e <- TNR_e + SL_e$TNR

      A_l[i, ] <- SL_l$beta
      acc_l <- acc_l + SL_l$acc
      TPR_l <- TPR_l + SL_l$TPR
      TNR_l <- TNR_l + SL_l$TNR
    }
    acc_o2 <- acc_o2 / (N - m)
    TPR_o2 <- TPR_o2 / (N - m)
    TNR_o2 <- TNR_o2 / (N - m)
    acc_e <- acc_e / (N - m)
    TPR_e <- TPR_e / (N - m)
    TNR_e <- TNR_e / (N - m)
    acc_l <- acc_l / (N - m)
    TPR_l <- TPR_l / (N - m)
    TNR_l <- TNR_l / (N - m)

    # Update estimated degrees
    nr1 <- rowSums(A_o2) + rowSums(t(A12))
    nr2 <- rowSums(A_e) + rowSums(t(A12))
    nr3 <- rowSums(A_l) + rowSums(t(A12))
    nr1[nr1 == 0] <- (5 / N + 1 / (3 * N^(0.45))) * (N - m) / 2
    nr2[nr2 == 0] <- (5 / N + 1 / (3 * N^(0.45))) * (N - m) / 2
    nr3[nr3 == 0] <- (5 / N + 1 / (3 * N^(0.45))) * (N - m) / 2
  }

  return(list(
    acc_o1 = acc_o1, TPR_o1 = TPR_o1, TNR_o1 = TNR_o1,
    acc_o2 = acc_o2, TPR_o2 = TPR_o2, TNR_o2 = TNR_o2,
    acc_e = acc_e, TPR_e = TPR_e, TNR_e = TNR_e,
    acc_l = acc_l, TPR_l = TPR_l, TNR_l = TNR_l,
    acc_SMC = acc_SMC, TPR_SMC = TPR_SMC, TNR_SMC = TNR_SMC
  ))
}

# Run simulation with parallel processing
cat("Simulation start time:", as.character(Sys.time()), "\n")

total <- ran_num

# Create progress log file
progress_log <- "progress_ER.txt"
if (file.exists(progress_log)) file.remove(progress_log)

result <- pbmclapply(1:ran_num, function(l) {
  res <- simulation_ER(N, T, par)
  cat("1", file = progress_log, append = TRUE)
  done_count <- file.info(progress_log)$size
  message(paste0("Done: index ", l, " | Completed: ", done_count, "/", ran_num, " | ", Sys.time()))
  return(res)
}, mc.cores = 20, mc.preschedule = FALSE)

# Remove progress log file
if (file.exists(progress_log)) file.remove(progress_log)

# Aggregate results
result1 <- 0
result2 <- matrix(0, ran_num, 15)
for (i in 1:ran_num) {
  result1 <- result1 + as.numeric(result[[i]])
  result2[i, ] <- as.numeric(result[[i]])
}
result1 <- result1 / ran_num
write.table(result1,
  file = paste0("result_ER_N=", N, "_T=", T, ".txt"),
  sep = "\t", row.names = FALSE, col.names = FALSE
)






