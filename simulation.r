  # simulation.R - Monte Carlo Simulations
  # Based on Psaradakis et al. (2024)
  
  #--------------------------------------------------------------
  # Bias Simulation
  #--------------------------------------------------------------

  run_bias_simulation <- function(Phi, T_values, max_horizon, n_simulations, dgp_name) {
    # Calculate true impulse responses
    true_irf <- calculate_true_irf(Phi, max_horizon)
    
    # Store bias and RMSE results for each sample size
    bias_results <- list()
    rmse_results <- list()
    
    # Loop over sample sizes
    for (T_value in T_values) {
      cat("Running bias and RMSE simulation for", dgp_name, "with T =", T_value, "\n")
      
      # Initialize bias accumulators
      var_bias_R11 <- numeric(max_horizon + 1)
      lp_bias_R11 <- numeric(max_horizon + 1)
      var_bias_R21 <- numeric(max_horizon + 1)
      lp_bias_R21 <- numeric(max_horizon + 1)
      var_bias_R12 <- numeric(max_horizon + 1)
      lp_bias_R12 <- numeric(max_horizon + 1)
      var_bias_R22 <- numeric(max_horizon + 1)
      lp_bias_R22 <- numeric(max_horizon + 1)
      
      # Initialize RMSE accumulators
      var_rmse_R11 <- numeric(max_horizon + 1)
      lp_rmse_R11 <- numeric(max_horizon + 1)
      var_rmse_R21 <- numeric(max_horizon + 1)
      lp_rmse_R21 <- numeric(max_horizon + 1)
      var_rmse_R12 <- numeric(max_horizon + 1)
      lp_rmse_R12 <- numeric(max_horizon + 1)
      var_rmse_R22 <- numeric(max_horizon + 1)
      lp_rmse_R22 <- numeric(max_horizon + 1)
      
      # Track squared errors for RMSE calculation
      var_se_R11 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      lp_se_R11 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      var_se_R21 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      lp_se_R21 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      var_se_R12 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      lp_se_R12 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      var_se_R22 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      lp_se_R22 <- matrix(0, nrow = n_simulations, ncol = max_horizon + 1)
      
      # Valid simulation counter
      valid_sims <- 0
      
      # Run simulations
      for (sim in 1:n_simulations) {
        # Generate data
        data <- generate_var_data(Phi, T_value, seed = 1000 + sim)
        
        # Estimate VAR and LP IRFs
        var_ok <- TRUE
        lp_ok <- TRUE
        
        tryCatch({
          var_irf <- estimate_var_irf(data, max_horizon)
        }, error = function(e) {
          warning("VAR estimation failed in simulation ", sim, ": ", e$message)
          var_ok <<- FALSE
        })
        
        tryCatch({
          lp_irf <- estimate_lp_irf(data, max_horizon)
        }, error = function(e) {
          warning("LP estimation failed in simulation ", sim, ": ", e$message)
          lp_ok <<- FALSE
        })
        
        # If VAR and LP estimation successful
        if (var_ok && lp_ok) {
          # Accumulate bias
          var_bias_R11 <- var_bias_R11 + (var_irf$e1[, 1] - true_irf$e1[, 1])
          lp_bias_R11 <- lp_bias_R11 + (lp_irf$e1[, 1] - true_irf$e1[, 1])
          var_bias_R21 <- var_bias_R21 + (var_irf$e1[, 2] - true_irf$e1[, 2])
          lp_bias_R21 <- lp_bias_R21 + (lp_irf$e1[, 2] - true_irf$e1[, 2])
          var_bias_R12 <- var_bias_R12 + (var_irf$e2[, 1] - true_irf$e2[, 1])
          lp_bias_R12 <- lp_bias_R12 + (lp_irf$e2[, 1] - true_irf$e2[, 1])
          var_bias_R22 <- var_bias_R22 + (var_irf$e2[, 2] - true_irf$e2[, 2])
          lp_bias_R22 <- lp_bias_R22 + (lp_irf$e2[, 2] - true_irf$e2[, 2])
          
          # Store squared errors for RMSE
          var_se_R11[valid_sims + 1, ] <- (var_irf$e1[, 1] - true_irf$e1[, 1])^2
          lp_se_R11[valid_sims + 1, ] <- (lp_irf$e1[, 1] - true_irf$e1[, 1])^2
          var_se_R21[valid_sims + 1, ] <- (var_irf$e1[, 2] - true_irf$e1[, 2])^2
          lp_se_R21[valid_sims + 1, ] <- (lp_irf$e1[, 2] - true_irf$e1[, 2])^2
          var_se_R12[valid_sims + 1, ] <- (var_irf$e2[, 1] - true_irf$e2[, 1])^2
          lp_se_R12[valid_sims + 1, ] <- (lp_irf$e2[, 1] - true_irf$e2[, 1])^2
          var_se_R22[valid_sims + 1, ] <- (var_irf$e2[, 2] - true_irf$e2[, 2])^2
          lp_se_R22[valid_sims + 1, ] <- (lp_irf$e2[, 2] - true_irf$e2[, 2])^2
          
          valid_sims <- valid_sims + 1
        }
        
        # Progress reporting
        if (sim %% max(1, floor(n_simulations/10)) == 0) {
          cat(sprintf("  Completed simulation %d/%d (%.1f%%) for T = %d\n", 
                      sim, n_simulations, 100*sim/n_simulations, T_value))
        }
      }
      
      if (valid_sims > 0) {
        # Calculate average bias
        var_bias_R11 <- var_bias_R11 / valid_sims
        lp_bias_R11 <- lp_bias_R11 / valid_sims
        var_bias_R21 <- var_bias_R21 / valid_sims
        lp_bias_R21 <- lp_bias_R21 / valid_sims
        var_bias_R12 <- var_bias_R12 / valid_sims
        lp_bias_R12 <- lp_bias_R12 / valid_sims
        var_bias_R22 <- var_bias_R22 / valid_sims
        lp_bias_R22 <- lp_bias_R22 / valid_sims
        
        # Calculate RMSE from squared errors
        var_rmse_R11 <- sqrt(colMeans(var_se_R11[1:valid_sims, ]))
        lp_rmse_R11 <- sqrt(colMeans(lp_se_R11[1:valid_sims, ]))
        var_rmse_R21 <- sqrt(colMeans(var_se_R21[1:valid_sims, ]))
        lp_rmse_R21 <- sqrt(colMeans(lp_se_R21[1:valid_sims, ]))
        var_rmse_R12 <- sqrt(colMeans(var_se_R12[1:valid_sims, ]))
        lp_rmse_R12 <- sqrt(colMeans(lp_se_R12[1:valid_sims, ]))
        var_rmse_R22 <- sqrt(colMeans(var_se_R22[1:valid_sims, ]))
        lp_rmse_R22 <- sqrt(colMeans(lp_se_R22[1:valid_sims, ]))
        
        # Store bias results
        bias_results[[paste0("T", T_value)]] <- list(
          var_bias_R11 = var_bias_R11,
          lp_bias_R11 = lp_bias_R11,
          var_bias_R21 = var_bias_R21,
          lp_bias_R21 = lp_bias_R21,
          var_bias_R12 = var_bias_R12,
          lp_bias_R12 = lp_bias_R12,
          var_bias_R22 = var_bias_R22,
          lp_bias_R22 = lp_bias_R22
        )
        
        # Store RMSE results
        rmse_results[[paste0("T", T_value)]] <- list(
          var_rmse_R11 = var_rmse_R11,
          lp_rmse_R11 = lp_rmse_R11,
          var_rmse_R21 = var_rmse_R21,
          lp_rmse_R21 = lp_rmse_R21,
          var_rmse_R12 = var_rmse_R12,
          lp_rmse_R12 = lp_rmse_R12,
          var_rmse_R22 = var_rmse_R22,
          lp_rmse_R22 = lp_rmse_R22
        )
        
        # Save intermediate results
        intermediate_results <- list(
          true_irf = true_irf,
          bias_results = bias_results,
          rmse_results = rmse_results
        )
        save(intermediate_results, file = paste0("output/bias_results_", dgp_name, "_T", T_value, "_intermediate.RData"))
        cat("  Saved intermediate results for T =", T_value, "\n")
      } else {
        warning("No valid simulations for T =", T_value)
      }
    }
    
    # Return true IRF, bias results, and RMSE results
    return(list(
      true_irf = true_irf,
      bias_results = bias_results,
      rmse_results = rmse_results
    ))
  }
  
  #--------------------------------------------------------------
  # Predictive Accuracy Simulation (Diebold-Mariano)
  #--------------------------------------------------------------
  
  # Generate differentials for recursive scheme
  generate_dm_differentials <- function(true_irf, T_value, T0, dgp_fn, Phi, horizons, e_idx, y_idx, step = 10) {
    # Calculate number of steps
    N <- floor((T_value - T0) / step) + 1
    differentials <- matrix(NA, nrow = N, ncol = horizons)
    
    for (s in 1:N) {
      # Data length for this iteration
      current_T <- T0 + (s-1) * step
      
      # Generate data with fixed seed for reproducibility
      data <- dgp_fn(Phi, sample_size = current_T, seed = 123 + s)
      
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
  
  # Main DM simulation function with fixed number of steps
  run_dm_simulation <- function(Phi, T_values, max_horizon, n_simulations, dgp_name, T0 = 60, 
                                fixed_steps = 15, use_fixed_steps = TRUE) {
    # Calculate true impulse responses (only once)
    true_irf <- calculate_true_irf(Phi, max_horizon)
    dm_results <- list()
    
    # Detect and create parallel clusters
    cl <- makeCluster(min(detectCores() - 1, 8))  # Limit to 8 cores max
    
    # Load required packages on all workers
    clusterEvalQ(cl, {
      library(vars)
      library(lpirfs)
      library(stats)
    })
    
    # Export necessary objects and functions to clusters
    clusterExport(cl, c("dm_test_recursive", "generate_var_data",
                        "calculate_true_irf", "estimate_var_irf", "estimate_lp_irf",
                        "Phi", "true_irf", "T0", "fixed_steps", "use_fixed_steps"), 
                  envir = environment())
    
    for (T_value in T_values) {
      cat("Running DM simulation for", dgp_name, "with T =", T_value, "\n")
      
      # Calculate step size based on whether we want fixed steps or not
      if (use_fixed_steps) {
        # Adjust step size to ensure consistent number of steps across all T values
        step <- max(1, floor((T_value - T0) / fixed_steps))
        steps <- fixed_steps
        cat("  Using fixed", fixed_steps, "steps with step size", step, "\n")
      } else {
        # Original behavior - step is fixed and steps varies with T
        step <- 10
        steps <- floor((T_value - T0) / step) + 1
        cat("  Using variable", steps, "steps with fixed step size", step, "\n")
      }
      
      # Export step to all worker nodes
      clusterExport(cl, "step", envir = environment())
      
      dm_sim_results <- parLapply(cl, 1:n_simulations, function(sim) {
        # Pre-generate data for all recursive steps with unique seeds
        all_data <- list()
        for (s in 1:steps) {
          # Calculate current sample size based on step size
          if (use_fixed_steps) {
            # Evenly distribute sample sizes between T0 and T_value
            current_T <- T0 + round((s-1) * (T_value - T0) / (steps - 1))
          } else {
            # Original calculation
            current_T <- T0 + (s-1) * step
          }
          
          # Ensure we don't exceed T_value
          current_T <- min(current_T, T_value)
          
          # Generate data with unique seed
          all_data[[s]] <- generate_var_data(Phi, sample_size = current_T, seed = 1000 + sim*100 + s)
        }
        
        # Initialize matrices for all differentials
        diff_R11 <- matrix(NA, nrow = steps, ncol = max_horizon)
        diff_R21 <- matrix(NA, nrow = steps, ncol = max_horizon)
        diff_R12 <- matrix(NA, nrow = steps, ncol = max_horizon)
        diff_R22 <- matrix(NA, nrow = steps, ncol = max_horizon)
        
        # Fill differential matrices for all recursive steps
        for (s in 1:steps) {
          # Get data for this step
          data <- all_data[[s]]
          
          # Estimate models only once per data set
          var_irf <- estimate_var_irf(data, max_horizon)
          lp_irf <- estimate_lp_irf(data, max_horizon)
          
          # Calculate all differentials from a single estimation
          for (horizon in 1:max_horizon) {
            # R11: y1 response to e1 shock
            true_val_R11 <- true_irf$e1[horizon, 1]
            var_val_R11 <- var_irf$e1[horizon, 1]
            lp_val_R11 <- lp_irf$e1[horizon, 1]
            diff_R11[s, horizon] <- (true_val_R11 - lp_val_R11)^2 - (true_val_R11 - var_val_R11)^2
            
            # R21: y2 response to e1 shock
            true_val_R21 <- true_irf$e1[horizon, 2]
            var_val_R21 <- var_irf$e1[horizon, 2]
            lp_val_R21 <- lp_irf$e1[horizon, 2]
            diff_R21[s, horizon] <- (true_val_R21 - lp_val_R21)^2 - (true_val_R21 - var_val_R21)^2
            
            # R12: y1 response to e2 shock
            true_val_R12 <- true_irf$e2[horizon, 1]
            var_val_R12 <- var_irf$e2[horizon, 1]
            lp_val_R12 <- lp_irf$e2[horizon, 1]
            diff_R12[s, horizon] <- (true_val_R12 - lp_val_R12)^2 - (true_val_R12 - var_val_R12)^2
            
            # R22: y2 response to e2 shock
            true_val_R22 <- true_irf$e2[horizon, 2]
            var_val_R22 <- var_irf$e2[horizon, 2]
            lp_val_R22 <- lp_irf$e2[horizon, 2]
            diff_R22[s, horizon] <- (true_val_R22 - lp_val_R22)^2 - (true_val_R22 - var_val_R22)^2
          }
        }
        
        # Clear data to free memory
        all_data <- NULL
        gc(verbose = FALSE)
        
        # Test differentials for each horizon
        rejection_counts <- matrix(0, nrow = 4, ncol = max_horizon)
        for (horizon in 1:max_horizon) {
          # Run the DM test for each response-shock combination
          test_R11 <- dm_test_recursive(diff_R11, horizon)
          rejection_counts[1, horizon] <- as.numeric(!is.na(test_R11$p_value) && test_R11$p_value < 0.05)
          
          test_R21 <- dm_test_recursive(diff_R21, horizon)
          rejection_counts[2, horizon] <- as.numeric(!is.na(test_R21$p_value) && test_R21$p_value < 0.05)
          
          test_R12 <- dm_test_recursive(diff_R12, horizon)
          rejection_counts[3, horizon] <- as.numeric(!is.na(test_R12$p_value) && test_R12$p_value < 0.05)
          
          test_R22 <- dm_test_recursive(diff_R22, horizon)
          rejection_counts[4, horizon] <- as.numeric(!is.na(test_R22$p_value) && test_R22$p_value < 0.05)
        }
        
        if (sim %% max(1, floor(n_simulations/10)) == 0) {
          cat(sprintf("  Completed simulation %d/%d (%.1f%%) for T = %d\n", 
                      sim, n_simulations, 100*sim/n_simulations, T_value))
        }
        
        return(rejection_counts)
      })
      
      # Calculate rejection rates (percentage)
      rejection_R11 <- Reduce('+', lapply(dm_sim_results, function(x) x[1,])) / n_simulations * 100
      rejection_R21 <- Reduce('+', lapply(dm_sim_results, function(x) x[2,])) / n_simulations * 100
      rejection_R12 <- Reduce('+', lapply(dm_sim_results, function(x) x[3,])) / n_simulations * 100
      rejection_R22 <- Reduce('+', lapply(dm_sim_results, function(x) x[4,])) / n_simulations * 100
      
      # Save results for this sample size
      dm_results[[paste0("T", T_value)]] <- list(
        rejection_R11 = rejection_R11,
        rejection_R21 = rejection_R21,
        rejection_R12 = rejection_R12,
        rejection_R22 = rejection_R22
      )
      
      # Save intermediate results after each T value (for safety)
      intermediate_dm_results <- list(
        true_irf = true_irf,
        dm_results = dm_results,
        settings = list(use_fixed_steps = use_fixed_steps, fixed_steps = fixed_steps)
      )
      save(intermediate_dm_results, file = paste0("output/dm_results_", dgp_name, "_T", T_value, "_intermediate.RData"))
      cat("  Saved intermediate results for T =", T_value, "\n")
    }
    
    # Clean up and return
    stopCluster(cl)
    cat("Completed all DM simulations for", dgp_name, "\n")
    
    return(list(true_irf = true_irf, dm_results = dm_results))
  }