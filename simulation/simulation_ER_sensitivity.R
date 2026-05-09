# simulation_ER_sensitivity.R
# Sensitivity analysis: impact of n_i estimation error on Partial Network Recovery

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
T <- 100
ran_num <- 100
par <- c(0.2, 2)

# Sensitivity analysis: different levels of error in n_i estimation
error_levels <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)


#' Sensitivity analysis: impact of n_i estimation error on PNR
#'
#' @param N Number of nodes
#' @param T Number of time periods
#' @param par Vector of parameters: c(alpha, beta)
#' @param n_i_error Standard deviation of multiplicative noise
#' @return List containing accuracy, TPR, TNR for different methods
simulation_ER_sensitivity <- function(N, T, par, n_i_error = 0) {
  A <- get_A_ER(N)
  SMC <- SMC1(A, 1)
  A <- SMC$A
  m <- SMC$m
  A12 <- A[1:m, -c(1:m)]
  S <- SMC$S
  W <- A / rowSums(A)
  n_t <- rowSums(A)
  n_t <- n_t[-c(1:m)]

  # Apply error to n_i (multiplicative noise with normal distribution)
  if (n_i_error == 0) {
    nr <- n_t
  } else {
    nr_errors <- rnorm(length(n_t), mean = 0, sd = n_i_error)
    nr <- n_t * (1 + nr_errors)
    nr[nr <= 0] <- 1
    nr <- pmax(1, round(nr))
  }

  M1 <- gen_y(N, T, par, W, S)
  Y <- M1$Y

  # Compute with given n_i (either true or with error)
  acc <- TPR <- TNR <- 0
  for (i in 1:(N - m)) {
    Y_fit <- Y[(m + 1):N, ] - par[2] * (1 / nr[i]) * t(A12) %*% (Y[1:m, ]) - par[1]
    X_fit <- t(par[2] * (1 / nr[i]) * (Y[, ]))
    X_fit <- X_fit[, -c(1:m, i + m)]
    SL_fit <- SL(X_fit, Y_fit, i, N, m, A)

    acc <- acc + SL_fit$acc
    TPR <- TPR + SL_fit$TPR
    TNR <- TNR + SL_fit$TNR
  }
  acc <- acc / (N - m)
  TPR <- TPR / (N - m)
  TNR <- TNR / (N - m)

  return(list(acc = acc, TPR = TPR, TNR = TNR,
    acc_SMC = SMC$acc, TPR_SMC = SMC$TPR, TNR_SMC = SMC$TNR))
}


# Run sensitivity analysis
cat("Sensitivity Analysis: n_i estimation error impact on PNR\n")
cat("==================================================\n\n")

results_sensitivity <- data.frame()

for (err_level in error_levels) {
  cat("Testing n_i error level:", err_level * 100, "%\n")
  cat("Start time:", as.character(Sys.time()), "\n")

  progress_log <- paste0("progress_", err_level * 100, ".txt")
  if (file.exists(progress_log)) file.remove(progress_log)

  result <- pbmclapply(1:ran_num, function(l) {
    res <- simulation_ER_sensitivity(N, T, par, n_i_error = err_level)
    cat("1", file = progress_log, append = TRUE)
    return(res)
  }, mc.cores = 20, mc.preschedule = FALSE)

  if (file.exists(progress_log)) file.remove(progress_log)

  # Aggregate results
  result_avg <- 0
  for (i in 1:ran_num) {
    result_avg <- result_avg + as.numeric(result[[i]])
  }
  result_avg <- result_avg / ran_num

  results_sensitivity <- rbind(results_sensitivity,
    data.frame(error_level = err_level * 100,
      acc = result_avg[1],
      TPR = result_avg[2],
      TNR = result_avg[3]))

  cat("Error", err_level * 100, "% - Acc:", round(result_avg[1], 4),
    "TPR:", round(result_avg[2], 4), "TNR:", round(result_avg[3], 4), "\n\n")
}

# Save results
output_file <- paste0("sensitivity_ni_error_N", N, "_T", T, ".txt")
write.table(results_sensitivity, file = output_file, sep = "\t", row.names = FALSE)

cat("Results saved to:", output_file, "\n")