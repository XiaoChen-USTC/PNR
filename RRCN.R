RRCN <- function(Y2, i, N, n1, m, A_true, threshold = 1e-3) {
  # Y2: n1 x T data matrix (nodes x time)
  # i: index of target node (1 to N, where N = n1-m)
  # N: number of unknown nodes (= n1-m)
  # n1: total number of nodes
  # m: number of known nodes
  # A_true: true complete network adjacency matrix (n1 x n1)
  # threshold: threshold for determining connection existence

  # Phi = t(Y2): T x n1 (time x nodes)
  Phi = as.matrix(t(Y2))

  # Y_i: time series of target node i (length T)
  Y_i = as.numeric(Y2[i, ])
  T = length(Y_i)

  # Fit Lasso model with cross-validation (no intercept)
  cvfit0 = cv.glmnet(x = Phi, y = Y_i, nfolds = 5, intercept = FALSE, alpha = 1, type.measure = "mse")

  # Extract coefficients at optimal lambda (length n1)
  beta0 = as.numeric(coef(cvfit0, s = "lambda.min"))
  if(length(beta0) != n1) {
    beta0 = beta0[-1]  # Remove intercept term
  }

  # Set diagonal (self-connection) to 0
  beta0[i] = 0

  # Binarize based on threshold
  beta_binary = ifelse(abs(beta0) > threshold, 1, 0)

  # Calculate fitting loss
  pred = as.numeric(Phi %*% beta0)
  L_l = mean((Y_i - pred)^2)

  # Calculate evaluation metrics (only compare last n1-m edges, i.e., columns m+1 to n1)
  A_target = A_true[i, (m+1):n1]
  A_target[i - m] = 0

  beta_compare = beta_binary[(m+1):n1]

  M_compare = beta_compare + A_target

  # Accuracy
  acc = sum(beta_compare == A_target) / N

  # Avoid NA when denominator is 0
  P_true = sum(A_target == 1)
  N_true = sum(A_target == 0)

  TPR = ifelse(P_true > 0, sum(M_compare == 2) / P_true, 1)
  TNR = ifelse(N_true > 0, sum(M_compare == 0) / N_true, 1)

  return(list(beta = beta_binary, acc = acc, TPR = TPR, TNR = TNR, L_l = L_l))
}