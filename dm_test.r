# dm_test.R - Diebold-Mariano Test Implementation

#--------------------------------------------------------------
# Diebold-Mariano Test for Equal Predictive Accuracy
#--------------------------------------------------------------

# Calculate the modified Diebold-Mariano test statistic as in Harvey et al. (1997)
# Returns the test statistic and p-value
dm_test <- function(true_response, var_response, lp_response) {
  # Calculate the squared errors
  se_var <- (true_response - var_response)^2
  se_lp <- (true_response - lp_response)^2
  
  # Calculate the differentials (LP - VAR)
  d <- se_lp - se_var
  
  # Sample size
  n <- length(d)
  
  # Mean of differentials
  d_bar <- mean(d)
  
  # Estimate asymptotic variance (autocovariance up to tau-1 lags)
  tau <- min(n-1, 15)  # Using horizon as maximum lag, capped at n-1
  
  # Autocovariance calculation
  gamma0 <- sum((d - d_bar)^2) / n
  
  # Autocovariance for lags 1 to tau-1
  gamma <- numeric(tau-1)
  if (tau > 1) {
    for (k in 1:(tau-1)) {
      gamma[k] <- sum((d[(k+1):n] - d_bar) * (d[1:(n-k)] - d_bar)) / n
    }
  }
  
  # Asymptotic variance estimate
  omega <- gamma0 + 2 * sum(gamma)
  
  # Harvey et al. (1997) small-sample correction
  small_sample_factor <- (n + 1 - 2*tau + n^(-1)*tau*(tau-1)) / n
  
  # Modified Diebold-Mariano test statistic
  if (omega > 0) {
    dm_stat <- small_sample_factor * d_bar / sqrt(omega/n)
    
    # p-value for one-sided test (VAR better than LP if d_bar > 0)
    p_value <- pt(dm_stat, df = n-1, lower.tail = FALSE)
  } else {
    # Handle case where variance estimate is not positive
    dm_stat <- NA
    p_value <- NA
  }
  
  return(list(statistic = dm_stat, p_value = p_value, differential = d_bar))
}

# Function to compare VAR and LP accuracy across horizons
compare_accuracy <- function(true_irf, var_irf, lp_irf, shock_idx, response_idx) {
  # Number of horizons
  horizons <- nrow(true_irf$e1)
  
  # Extract true IRF values
  if (shock_idx == 1) {
    true_vals <- true_irf$e1[, response_idx]
    var_vals <- var_irf$e1[, response_idx]
    lp_vals <- lp_irf$e1[, response_idx]
  } else {
    true_vals <- true_irf$e2[, response_idx]
    var_vals <- var_irf$e2[, response_idx]
    lp_vals <- lp_irf$e2[, response_idx]
  }
  
  # Results container
  results <- data.frame(
    Horizon = 1:horizons,
    DM_Statistic = numeric(horizons),
    P_Value = numeric(horizons),
    VAR_Better = logical(horizons)
  )
  
  # Run test for each horizon
  for (tau in 1:horizons) {
    test_result <- dm_test(true_vals[tau], var_vals[tau], lp_vals[tau])
    results$DM_Statistic[tau] <- test_result$statistic
    results$P_Value[tau] <- test_result$p_value
    results$VAR_Better[tau] <- !is.na(test_result$p_value) && test_result$p_value < 0.05
  }
  
  return(results)
}

# Generate differentials for recursive estimation (as in the paper)
generate_dm_differentials <- function(true_irf, T, T0, dgp_fn, Phi, horizons, e_idx, y_idx) {
  N <- T - T0 + 1
  differentials <- matrix(NA, nrow = N, ncol = horizons)
  
  for (s in 1:N) {
    # Data length for this iteration
    current_T <- T0 + s - 1
    
    # Generate data with fixed seed for reproducibility
    data <- dgp_fn(Phi, current_T, seed = 123 + s)
    
    # Estimate VAR and LP IRFs
    var_irf <- estimate_var_irf(data, horizons)
    lp_irf <- estimate_lp_irf(data, horizons)
    
    # Extract true, VAR, and LP values
    if (e_idx == 1) {
      true_vals <- true_irf$e1[, y_idx]
      var_vals <- var_irf$e1[, y_idx]
      lp_vals <- lp_irf$e1[, y_idx]
    } else {
      true_vals <- true_irf$e2[, y_idx]
      var_vals <- var_irf$e2[, y_idx]
      lp_vals <- lp_irf$e2[, y_idx]
    }
    
    # Calculate squared error differentials
    se_var <- (true_vals - var_vals)^2
    se_lp <- (true_vals - lp_vals)^2
    
    # Store differentials (LP - VAR)
    differentials[s, ] <- se_lp - se_var
  }
  
  return(differentials)
}

# Perform Diebold-Mariano test with recursive differentials
dm_test_recursive <- function(differentials, horizon) {
  # Extract differentials for the specific horizon
  d <- differentials[, horizon]
  d <- d[!is.na(d)]  # Remove NAs if any
  
  # Sample size
  n <- length(d)
  if (n < 2) return(list(statistic = NA, p_value = NA))
  
  # Mean of differentials
  d_bar <- mean(d)
  
  # Estimate asymptotic variance (autocovariance up to tau-1 lags)
  tau <- min(n-1, horizon)
  
  # Autocovariance calculation
  gamma0 <- sum((d - d_bar)^2) / n
  
  # Autocovariance for lags 1 to tau-1
  gamma <- numeric(tau)
  if (tau > 0) {
    for (k in 1:tau) {
      if (k < n) {
        gamma[k] <- sum((d[(k+1):n] - d_bar) * (d[1:(n-k)] - d_bar)) / n
      }
    }
  }
  
  # Asymptotic variance estimate
  omega <- gamma0 + 2 * sum(gamma)
  
  # Harvey et al. (1997) small-sample correction
  small_sample_factor <- (n + 1 - 2*tau + n^(-1)*tau*(tau-1)) / n
  
  # Modified Diebold-Mariano test statistic
  if (omega > 0) {
    dm_stat <- small_sample_factor * d_bar / sqrt(omega/n)
    
    # p-value for one-sided test (VAR better than LP if d_bar > 0)
    p_value <- pt(dm_stat, df = n-1, lower.tail = FALSE)
  } else {
    # Handle case where variance estimate is not positive
    dm_stat <- NA
    p_value <- NA
  }
  
  return(list(statistic = dm_stat, p_value = p_value))
}