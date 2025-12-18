# Install and load required packages
if (!require("vars")) {
  install.packages("vars")
  library(vars)
}

if (!require("BVAR")) {
  install.packages("BVAR")
  library(BVAR)
}

if (!require("MCMCpack")) {
  install.packages("MCMCpack")
  library(MCMCpack)
}

# Generate synthetic bivariate time series data
set.seed(123)
n <- nrow(X_2D_NETWPROP_TIME_1)  # Number of observations

# Create two time series with Granger causality relationship
# X1 -> X2 (X1 causes X2)
X1 <- as.ts(as.vector(as.matrix(X_2D_NETWPROP_TIME_1$Assortativity)))
X2 <- as.ts(as.vector(as.matrix(X_2D_NETWPROP_TIME_1$Mean.PD)))

# Initialize first observations
#X1[1:2] <- rnorm(2, 0, 1)
#X2[1:2] <- rnorm(2, 0, 1)

# Generate time series with causality
#for (t in 3:n) {
#  X1[t] <- 0.6 * X1[t-1] - 0.2 * X1[t-2] + rnorm(1, 0, 0.5)
#  X2[t] <- 0.4 * X2[t-1] + 0.3 * X1[t-1] + 0.2 * X1[t-2] + rnorm(1, 0, 0.5)
#}

# Combine into data frame
data <- data.frame(X1 = X1, X2 = X2)
ts_data <- ts(data)

# Method 1: Using BVAR package (corrected)
cat("\n=== METHOD 1: Using BVAR package ===\n")

# Prepare data for BVAR package
bvar_data <- as.matrix(data)

# Set priors correctly
priors <- bv_priors(
  hyper = c("lambda", "alpha", "psi"),
  mn = bv_minnesota(lambda = bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5),
                    alpha = bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3),
                    psi = bv_psi(scale = 0.004, shape = 0.004, mode = "auto", min = "auto", max = "auto"))
)

# Estimate model
bvar_model <- bvar(bvar_data, lags = 2, priors = priors,
                   n_draw = 5000, n_burn = 2000, verbose = FALSE)

# Print summary
print(summary(bvar_model))

# Calculate Granger causality using forecast error variance decomposition
fevd_results <- fevd(bvar_model, horizon = 10)
print(fevd_results)

# Method 2: Using vars package with Bayesian estimation
cat("\n=== METHOD 2: Using vars package ===\n")

# Estimate VAR model
var_model <- VAR(ts_data, p = 2, type = "const")

# Classical Granger causality test (for comparison)
gc_test <- causality(var_model, cause = "X1")
print(gc_test)

gc_test2 <- causality(var_model, cause = "X2")
print(gc_test2)

# Method 3: Manual Bayesian Granger Causality with MCMCpack
cat("\n=== METHOD 3: Bayesian VAR with MCMCpack ===\n")

# Prepare data matrices
Y <- as.matrix(data[3:n, ])
X <- matrix(1, nrow = n-2, ncol = 1)  # constant term

# Add lagged variables
for (i in 1:2) {
  X <- cbind(X, as.matrix(data[(3-i):(n-i), ]))
}

# Bayesian estimation of VAR using MCMCpack
bayesian_var <- function(Y, X, n_sim = 5000, n_burn = 1000) {

  # Priors
  b0 <- rep(0, ncol(X))
  B0 <- diag(1, ncol(X)) * 0.1  # Prior precision

  # Storage for posterior draws
  n_coef <- ncol(X)
  beta_draws <- array(NA, c(n_sim, n_coef, ncol(Y)))
  sigma_draws <- array(NA, c(n_sim, ncol(Y), ncol(Y)))

  # Initialize
  beta <- matrix(0, n_coef, ncol(Y))
  Sigma <- diag(ncol(Y))

  for (i in 1:(n_sim + n_burn)) {

    # Sample beta (equation by equation)
    for (j in 1:ncol(Y)) {
      V_beta <- solve(B0 + (1/Sigma[j,j]) * t(X) %*% X)
      mu_beta <- V_beta %*% (B0 %*% b0 + (1/Sigma[j,j]) * t(X) %*% Y[,j])
      beta[,j] <- as.vector(mvrnorm(1, mu_beta, V_beta))
    }

    # Sample Sigma
    S <- t(Y - X %*% beta) %*% (Y - X %*% beta)
    Sigma <- riwish(nrow(Y) + ncol(Y) + 1, S + diag(ncol(Y)))

    # Store draws after burn-in
    if (i > n_burn) {
      beta_draws[i - n_burn, , ] <- beta
      sigma_draws[i - n_burn, , ] <- Sigma
    }
  }

  return(list(beta = beta_draws, sigma = sigma_draws))
}

# Load required function from MASS
if (!require("MASS")) {
  install.packages("MASS")
  library(MASS)
}

# Run Bayesian estimation
bayesian_results <- bayesian_var(Y, X, n_sim = 3000, n_burn = 1000)

# Calculate Granger causality probabilities
calculate_granger_probability <- function(beta_draws, cause_idx, effect_idx) {
  # For our setup: columns 2,4 are lags of X1; columns 3,5 are lags of X2

  n_draws <- dim(beta_draws)[1]
  causality_count <- 0

  for (i in 1:n_draws) {
    # Test if coefficients of cause variable in effect equation are jointly zero
    if (cause_idx == 1 && effect_idx == 2) {  # X1 -> X2
      coefs <- beta_draws[i, c(2, 4), 2]  # X1 lags in X2 equation
    } else if (cause_idx == 2 && effect_idx == 1) {  # X2 -> X1
      coefs <- beta_draws[i, c(3, 5), 1]  # X2 lags in X1 equation
    }

    # Check if any coefficient is significantly different from zero
    if (any(abs(coefs) > 0.1)) {  # Simple threshold rule
      causality_count <- causality_count + 1
    }
  }

  return(causality_count / n_draws)
}

# Calculate probabilities
prob_x1_to_x2 <- calculate_granger_probability(bayesian_results$beta, 1, 2)
prob_x2_to_x1 <- calculate_granger_probability(bayesian_results$beta, 2, 1)

cat("\nBayesian Granger Causality Probabilities:\n")
cat("P(X1 -> X2):", round(prob_x1_to_x2, 4), "\n")
cat("P(X2 -> X1):", round(prob_x2_to_x1, 4), "\n")

# Method 4: Simple Bayes Factor approach
cat("\n=== METHOD 4: Bayes Factor Approach ===\n")

bayes_factor_granger <- function(y1, y2, lags = 2) {

  n <- length(y1)

  # Prepare data for restricted and unrestricted models
  Y <- y2[(lags+1):n]

  # Unrestricted model: include lags of both y1 and y2
  X_unrest <- matrix(1, nrow = length(Y), ncol = 1)  # constant
  for (i in 1:lags) {
    X_unrest <- cbind(X_unrest,
                      y1[(lags+1-i):(n-i)],    # y1 lags
                      y2[(lags+1-i):(n-i)])    # y2 lags
  }

  # Restricted model: include only lags of y2
  X_rest <- matrix(1, nrow = length(Y), ncol = 1)   # constant
  for (i in 1:lags) {
    X_rest <- cbind(X_rest, y2[(lags+1-i):(n-i)])   # only y2 lags
  }

  # Calculate marginal likelihoods using conjugate normal-gamma prior
  calc_marginal_likelihood <- function(Y, X) {
    n <- length(Y)
    k <- ncol(X)

    # Prior parameters
    b_0 <- rep(0, k)
    V_0 <- diag(k) * 100  # diffuse prior
    a_0 <- 0.001
    d_0 <- 0.001

    # Posterior parameters
    V_n <- solve(solve(V_0) + t(X) %*% X)
    b_n <- V_n %*% (solve(V_0) %*% b_0 + t(X) %*% Y)
    a_n <- a_0 + n/2
    d_n <- d_0 + 0.5 * (t(Y) %*% Y + t(b_0) %*% solve(V_0) %*% b_0 - t(b_n) %*% solve(V_n) %*% b_n)

    # Log marginal likelihood
    log_ml <- lgamma(a_n) - lgamma(a_0) + a_0 * log(d_0) - a_n * log(d_n) +
      0.5 * (log(det(V_n)) - log(det(V_0))) - (n/2) * log(2*pi)

    return(as.numeric(log_ml))
  }

  # Calculate Bayes factor
  log_ml_unrest <- calc_marginal_likelihood(Y, X_unrest)
  log_ml_rest <- calc_marginal_likelihood(Y, X_rest)

  log_bf <- log_ml_unrest - log_ml_rest
  bf <- exp(log_bf)

  # Posterior probability (assuming equal prior probabilities)
  post_prob <- bf / (1 + bf)

  return(list(
    bayes_factor = bf,
    log_bayes_factor = log_bf,
    posterior_probability = post_prob
  ))
}

# Test both directions
bf_x1_to_x2 <- bayes_factor_granger(data$X1, data$X2, lags = 2)
bf_x2_to_x1 <- bayes_factor_granger(data$X2, data$X1, lags = 2)

cat("Bayes Factor Results:\n")
cat("X1 -> X2: BF =", round(bf_x1_to_x2$bayes_factor, 4),
    ", P =", round(bf_x1_to_x2$posterior_probability, 4), "\n")
cat("X2 -> X1: BF =", round(bf_x2_to_x1$bayes_factor, 4),
    ", P =", round(bf_x2_to_x1$posterior_probability, 4), "\n")

# Visualization
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(data$X1, type = "l", main = "Time Series X1", ylab = "X1", xlab = "Time", col = "blue")
plot(data$X2, type = "l", main = "Time Series X2", ylab = "X2", xlab = "Time", col = "red")
par(mfrow = c(1, 1))

cat("\n=== INTERPRETATION ===\n")
cat("Bayes Factor interpretation:\n")
cat("BF > 10: Strong evidence for causality\n")
cat("BF 3-10: Moderate evidence\n")
cat("BF 1-3: Weak evidence\n")
cat("BF < 1: Evidence against causality\n")
