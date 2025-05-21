# functions.R - Core Functions for Impulse Response Comparison


# Load required packages
library(vars)
library(lpirfs)
library(ggplot2)

#--------------------------------------------------------------
# Data Generation Functions
#--------------------------------------------------------------

# Generate VAR(1) data given parameter matrix
generate_var_data <- function(Phi, sample_size, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  K <- nrow(Phi)
  mu <- c(0.4, 0.4)  # From the paper
  y <- matrix(0, nrow = sample_size, ncol = K)
  
  # Initial values
  y[1, ] <- mu + rnorm(K)
  
  # Generate time series
  for (t in 2:sample_size) {
    y[t, ] <- mu + Phi %*% (y[t-1, ] - mu) + rnorm(K)
  }
  
  colnames(y) <- paste0("y", 1:K)
  return(y)
}

#--------------------------------------------------------------
# Impulse Response Functions
#--------------------------------------------------------------

# Calculate true impulse responses
# Modify calculate_true_irf to include horizon 0
calculate_true_irf <- function(Phi, horizons) {
  K <- nrow(Phi)
  # Include horizon 0 by adding +1 to matrix size
  true_irf_e1 <- matrix(0, nrow = horizons + 1, ncol = K)
  true_irf_e2 <- matrix(0, nrow = horizons + 1, ncol = K)
  
  # For τ = 0, the response is a unit shock
  true_irf_e1[1, ] <- c(1, 0)  # Shock e1 = (1,0)'
  true_irf_e2[1, ] <- c(0, 1)  # Shock e2 = (0,1)'
  
  # Calculate for τ = 1, 2, ..., horizons
  Phi_power <- diag(K)
  for (tau in 1:horizons) {
    Phi_power <- Phi_power %*% Phi
    true_irf_e1[tau+1, ] <- Phi_power %*% c(1, 0)
    true_irf_e2[tau+1, ] <- Phi_power %*% c(0, 1)
  }
  
  return(list(e1 = true_irf_e1, e2 = true_irf_e2))
}

# Estimate VAR-based IRFs
estimate_var_irf <- function(data, horizons, var_lag = 1) {
  # Extract data dimensions
  n <- nrow(data)
  k <- ncol(data)
  
  # Create matrices for regression
  Y <- data[(var_lag+1):n, ]
  X <- cbind(1, data[var_lag:(n-1), ])
  
  # Estimate VAR coefficients
  B <- solve(t(X) %*% X) %*% t(X) %*% Y
  const <- B[1, ]
  A <- t(B[-1, ])  # VAR coefficient matrix (transpose to get k×k)
  
  # Check stability - print warning if unstable
  eigvals <- eigen(A)$values
  if (any(abs(eigvals) >= 0.999)) {
    warning("Estimated VAR may be unstable: max eigenvalue = ", max(abs(eigvals)))
  }
  
  # Initialize impulse responses (including horizon 0)
  irf_e1 <- matrix(0, nrow = horizons + 1, ncol = k)
  irf_e2 <- matrix(0, nrow = horizons + 1, ncol = k)
  
  # For τ = 0, response to shock is 1 for same variable, 0 for others
  irf_e1[1, ] <- c(1, 0)
  irf_e2[1, ] <- c(0, 1)
  
  # For τ = 1, 2, ..., horizons, calculate I*A^τ
  A_power <- diag(k)
  for (tau in 1:horizons) {
    A_power <- A_power %*% A
    irf_e1[tau+1, ] <- A_power %*% c(1, 0)
    irf_e2[tau+1, ] <- A_power %*% c(0, 1)
  }
  
  return(list(e1 = irf_e1, e2 = irf_e2))
}

# Estimate LP-based IRFs
estimate_lp_irf <- function(data, horizons, lp_lag = 1) {
  n <- nrow(data)
  k <- ncol(data)
  
  # Initialize impulse responses (including horizon 0)
  irf_e1 <- matrix(0, nrow = horizons + 1, ncol = k)
  irf_e2 <- matrix(0, nrow = horizons + 1, ncol = k)
  
  # For τ = 0, response to shock is 1 for same variable, 0 for others
  irf_e1[1, ] <- c(1, 0)
  irf_e2[1, ] <- c(0, 1)
  
  # For each horizon τ = 1, 2, ..., horizons
  for (tau in 1:horizons) {
    # Create dependent and independent variables
    if (n < lp_lag + tau + 1) {
      # Not enough data points
      warning("Sample too small for horizon ", tau)
      # Copy previous horizon values as fallback
      if (tau > 1) {
        irf_e1[tau+1, ] <- irf_e1[tau, ]
        irf_e2[tau+1, ] <- irf_e2[tau, ]
      }
      next
    }
    
    # Create y_{t+τ}
    Y_future <- data[(lp_lag+tau):n, ]
    
    # Create lagged predictors including constant
    X <- cbind(1, data[lp_lag:(n-tau), ])
    
    # Ensure enough observations
    if (nrow(Y_future) < 10) {
      warning("Too few observations for horizon ", tau)
      next
    }
    
    # Separate regressions for each dependent variable
    for (i in 1:k) {
      # Use safe regression
      fit <- tryCatch({
        lm(Y_future[, i] ~ X - 1)
      }, error = function(e) {
        warning("LP regression failed for horizon ", tau, ", variable ", i, ": ", e$message)
        NULL
      })
      
      if (!is.null(fit)) {
        # Extract coefficients (skip intercept)
        irf_e1[tau+1, i] <- coef(fit)[2]  # Effect of y1 on variable i
        irf_e2[tau+1, i] <- coef(fit)[3]  # Effect of y2 on variable i
      }
    }
  }
  
  return(list(e1 = irf_e1, e2 = irf_e2))
}

#--------------------------------------------------------------
# Accuracy Metrics
#--------------------------------------------------------------

# Calculate accuracy metrics (MSE, MAE, RMSE, Bias)
calculate_accuracy <- function(true_vals, estimated_vals) {
  # MSE (Mean Squared Error)
  mse <- mean((true_vals - estimated_vals)^2)
  # MAE (Mean Absolute Error)
  mae <- mean(abs(true_vals - estimated_vals))
  # RMSE (Root Mean Squared Error)
  rmse <- sqrt(mse)
  # Bias (μεροληψία)
  bias <- mean(estimated_vals - true_vals)
  
  return(list(MSE = mse, MAE = mae, RMSE = rmse, Bias = bias))
}

#--------------------------------------------------------------
# Plotting Functions
#--------------------------------------------------------------

# Create dataframe for plotting
create_plot_data <- function(true_irf, var_irf, lp_irf, T_value, dgp_name, e_idx, y_idx) {
  # Select correct vector for shock e_idx and response y_idx
  if (e_idx == 1) {
    true_vals <- true_irf$e1[, y_idx]
    var_vals <- var_irf$e1[, y_idx]
    lp_vals <- lp_irf$e1[, y_idx]
    shock_name <- "e1"
  } else {
    true_vals <- true_irf$e2[, y_idx]
    var_vals <- var_irf$e2[, y_idx]
    lp_vals <- lp_irf$e2[, y_idx]
    shock_name <- "e2"
  }
  
  response_name <- paste0("y", y_idx)
  
  # Create dataframe
  plot_data <- data.frame(
    Horizon = rep(0:(length(true_vals)-1), 3),
    Method = c(rep("True IRF", length(true_vals)), 
               rep("VAR", length(true_vals)), 
               rep("LP", length(true_vals))),
    IRF = c(true_vals, var_vals, lp_vals),
    T_value = rep(T_value, 3 * length(true_vals)),
    DGP = rep(dgp_name, 3 * length(true_vals)),
    Shock = rep(shock_name, 3 * length(true_vals)),
    Response = rep(response_name, 3 * length(true_vals))
  )
  
  return(plot_data)
}

# Plot bias comparison
plot_bias_comparison <- function(bias_data_var, bias_data_lp, dgp_name, sample_sizes) {
  # Combine bias data
  plot_data <- rbind(
    data.frame(Horizon = rep(1:length(bias_data_var), length(sample_sizes)),
               Method = "VAR",
               Bias = unlist(bias_data_var),
               Sample = rep(paste0("T=", sample_sizes), each=length(bias_data_var[1]))),
    data.frame(Horizon = rep(1:length(bias_data_lp), length(sample_sizes)),
               Method = "LP",
               Bias = unlist(bias_data_lp),
               Sample = rep(paste0("T=", sample_sizes), each=length(bias_data_lp[1])))
  )
  
  # Create plot
  ggplot(plot_data, aes(x = Horizon, y = Bias, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    facet_wrap(~ Sample, scales = "free_y") +
    labs(title = paste0("Bias Comparison: ", dgp_name),
         x = "Horizon (τ)",
         y = "Bias") +
    theme_minimal() +
    scale_color_manual(values = c("VAR" = "blue", "LP" = "red"))
}