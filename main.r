# main.R - Main Script to Run All Components
# Based on Psaradakis et al. (2024)

#--------------------------------------------------------------
# Load Required Packages
#--------------------------------------------------------------
# Check and install missing packages
required_packages <- c("vars","grid", "lpirfs", "ggplot2", "gridExtra", "reshape2", "dplyr", "parallel")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Script Parameters
#--------------------------------------------------------------
# Simulation parameters
Phi_DGP1 <- matrix(c(0.49, 0.36, 0.36, 0.49), nrow = 2, byrow = TRUE)
Phi_DGP2 <- matrix(c(0.80, 0.10, 0.75, 0.40), nrow = 2, byrow = TRUE)

# reduce data for example 
use_reduced_samples <- FALSE  # Set to TRUE to use only T=100, 400, 1600 (faster)

T_values <- if(use_reduced_samples) c(100, 400, 1600) else c(100, 200, 400, 800, 1600)  # Sample sizes

max_horizon <- 15  # Maximum impulse response horizon
n_simulations <- 1000  # Number of Monte Carlo iterations
T0 <- 60  # Initial sample size for recursive estimation

# Control which simulations to run
run_bias_sim <- TRUE     # Set to FALSE to skip bias simulations
run_dm_sim <- TRUE       # Set to FALSE to skip DM test simulations
use_small_sample <- FALSE # Set to TRUE for quicker testing (with fewer iterations)
use_saved_results <- FALSE # Set to TRUE to load saved results instead of running simulations

# Output directory
dir.create("output", showWarnings = FALSE)

#--------------------------------------------------------------
# Source Required Scripts
#--------------------------------------------------------------
cat("Loading core functions...\n")
source("functions.R")  # basic functions
source("dm_test.R")    # DM test functions 
source("analysis.R")   # analysis functions (moved up before usage)

# Load simulation.R which have all simulation functions
cat("Loading simulation functions...\n")
source("simulation.R")

#--------------------------------------------------------------
# Preparation of parameters according to user settings
#--------------------------------------------------------------
# Adjust sample sizes and iterations for quick testing
if (use_small_sample) {
  cat("Using reduced sample for quick testing\n")
  T_values <- c(100, 400)
  n_simulations <- 100
  
  # User notification about the reduced parameters
  cat("  Reduced Monte Carlo iterations to", n_simulations, "\n")
  cat("  Using only sample sizes:", paste(T_values, collapse = ", "), "\n")
}

#--------------------------------------------------------------
# Run Simulations or Load Saved Results
#--------------------------------------------------------------
# Create a list to store all the results
all_results <- list()

# Initialize results variables (in case they're not created in the simulation steps)
bias_results_DGP1 <- NULL
bias_results_DGP2 <- NULL
dm_results_DGP1 <- NULL
dm_results_DGP2 <- NULL

# Check if saved simulation results exist and load them if requested
if (use_saved_results) {
  cat("Checking for saved simulation results...\n")
  
  # Check for full results file
  sim_results_file <- "output/simulation_results.RData"
  if (file.exists(sim_results_file)) {
    cat("Loading saved simulation results from", sim_results_file, "\n")
    load(sim_results_file)
    cat("Successfully loaded all simulation results.\n")
  } else {
    # Check for individual results files if full results file doesn't exist
    bias_file <- "output/bias_results.RData"
    dm_file <- "output/dm_results.RData"
    
    if (run_bias_sim && file.exists(bias_file)) {
      cat("Loading saved bias simulation results from", bias_file, "\n")
      load(bias_file)
      cat("Successfully loaded bias simulation results.\n")
    }
    
    if (run_dm_sim && file.exists(dm_file)) {
      cat("Loading saved DM test simulation results from", dm_file, "\n")
      load(dm_file)
      cat("Successfully loaded DM test simulation results.\n")
    }
    
    # If neither individual files exist, try to load from intermediate files
    if (run_bias_sim && !exists("bias_results_DGP1")) {
      cat("Full bias results not found. Attempting to reconstruct from intermediate files...\n")
      bias_results_DGP1 <- list(true_irf = NULL, bias_results = list(), rmse_results = list())
      bias_results_DGP2 <- list(true_irf = NULL, bias_results = list(), rmse_results = list())
      
      # Look for intermediate bias results
      reconstructed_bias <- FALSE
      for (dgp_name in c("DGP-1", "DGP-2")) {
        for (T_value in T_values) {
          intermediate_file <- paste0("output/bias_results_", dgp_name, "_T", T_value, "_intermediate.RData")
          if (file.exists(intermediate_file)) {
            load(intermediate_file)
            if (dgp_name == "DGP-1") {
              if (is.null(bias_results_DGP1$true_irf)) bias_results_DGP1$true_irf <- intermediate_results$true_irf
              bias_results_DGP1$bias_results[[paste0("T", T_value)]] <- 
                intermediate_results$bias_results[[paste0("T", T_value)]]
              if (!is.null(intermediate_results$rmse_results)) {
                bias_results_DGP1$rmse_results[[paste0("T", T_value)]] <- 
                  intermediate_results$rmse_results[[paste0("T", T_value)]]
              }
              reconstructed_bias <- TRUE
            } else {
              if (is.null(bias_results_DGP2$true_irf)) bias_results_DGP2$true_irf <- intermediate_results$true_irf
              bias_results_DGP2$bias_results[[paste0("T", T_value)]] <- 
                intermediate_results$bias_results[[paste0("T", T_value)]]
              if (!is.null(intermediate_results$rmse_results)) {
                bias_results_DGP2$rmse_results[[paste0("T", T_value)]] <- 
                  intermediate_results$rmse_results[[paste0("T", T_value)]]
              }
              reconstructed_bias <- TRUE
            }
            cat("  Loaded intermediate bias results for", dgp_name, "T =", T_value, "\n")
          }
        }
      }
      
      if (reconstructed_bias) {
        cat("Successfully reconstructed bias results from intermediate files.\n")
      }
    }
    
    if (run_dm_sim && !exists("dm_results_DGP1")) {
      cat("Full DM test results not found. Attempting to reconstruct from intermediate files...\n")
      dm_results_DGP1 <- list(true_irf = NULL, dm_results = list())
      dm_results_DGP2 <- list(true_irf = NULL, dm_results = list())
      
      # Look for intermediate DM test results
      reconstructed_dm <- FALSE
      for (dgp_name in c("DGP-1", "DGP-2")) {
        for (T_value in T_values) {
          intermediate_file <- paste0("output/dm_results_", dgp_name, "_T", T_value, "_intermediate.RData")
          if (file.exists(intermediate_file)) {
            load(intermediate_file)
            if (dgp_name == "DGP-1") {
              if (is.null(dm_results_DGP1$true_irf)) dm_results_DGP1$true_irf <- intermediate_dm_results$true_irf
              dm_results_DGP1$dm_results[[paste0("T", T_value)]] <- 
                intermediate_dm_results$dm_results[[paste0("T", T_value)]]
              reconstructed_dm <- TRUE
            } else {
              if (is.null(dm_results_DGP2$true_irf)) dm_results_DGP2$true_irf <- intermediate_dm_results$true_irf
              dm_results_DGP2$dm_results[[paste0("T", T_value)]] <- 
                intermediate_dm_results$dm_results[[paste0("T", T_value)]]
              reconstructed_dm <- TRUE
            }
            cat("  Loaded intermediate DM test results for", dgp_name, "T =", T_value, "\n")
          }
        }
      }
      
      if (reconstructed_dm) {
        cat("Successfully reconstructed DM test results from intermediate files.\n")
      }
    }
  }
}

# Run bias simulations if required and not loaded from saved results
if (run_bias_sim && (is.null(bias_results_DGP1) || is.null(bias_results_DGP2))) {
  cat("Running bias simulations...\n")
  
  # DGP-1
  cat("- For DGP-1\n")
  bias_results_DGP1 <- run_bias_simulation(Phi_DGP1, T_values, max_horizon, n_simulations, "DGP-1")
  all_results$bias_results_DGP1 <- bias_results_DGP1
  
  # DGP-2
  cat("- For DGP-2\n")
  bias_results_DGP2 <- run_bias_simulation(Phi_DGP2, T_values, max_horizon, n_simulations, "DGP-2")
  all_results$bias_results_DGP2 <- bias_results_DGP2
  
  # Saving intermediate results
  save(Phi_DGP1, Phi_DGP2, T_values, max_horizon, n_simulations, T0,
       bias_results_DGP1, bias_results_DGP2,
       file = "output/bias_results.RData")
  cat("Bias simulation results saved to output/bias_results.RData\n")
}

# Run Diebold-Mariano test simulations if required and not loaded from saved results
if (run_dm_sim && (is.null(dm_results_DGP1) || is.null(dm_results_DGP2))) {
  cat("Running Diebold-Mariano test simulations...\n")
  
  # DGP-1
  cat("- For DGP-1\n")
  dm_results_DGP1 <- run_dm_simulation(Phi_DGP1, T_values, max_horizon, n_simulations, "DGP-1", T0)
  all_results$dm_results_DGP1 <- dm_results_DGP1
  
  # DGP-2
  cat("- For DGP-2\n")
  dm_results_DGP2 <- run_dm_simulation(Phi_DGP2, T_values, max_horizon, n_simulations, "DGP-2", T0)
  all_results$dm_results_DGP2 <- dm_results_DGP2
  
  # Saving intermediate results
  save(Phi_DGP1, Phi_DGP2, T_values, max_horizon, n_simulations, T0,
       dm_results_DGP1, dm_results_DGP2,
       file = "output/dm_results.RData")
  cat("DM test results saved to output/dm_results.RData\n")
}

# Save all simulation results
save(Phi_DGP1, Phi_DGP2, T_values, max_horizon, n_simulations, T0,
     bias_results_DGP1, bias_results_DGP2, 
     dm_results_DGP1, dm_results_DGP2, 
     file = "output/simulation_results.RData")

cat("All simulation results saved to output/simulation_results.RData\n")

#--------------------------------------------------------------
# Analysis and Visualization
#--------------------------------------------------------------
cat("Analyzing results and creating visualizations...\n")

# Create bias plots
if (run_bias_sim && !is.null(bias_results_DGP1) && !is.null(bias_results_DGP2)) {
  cat("Creating bias plots...\n")
  
  # DGP-1
  bias_plots_DGP1_R11 <- plot_bias(bias_results_DGP1, "DGP-1", 1, 1)
  bias_plots_DGP1_R21 <- plot_bias(bias_results_DGP1, "DGP-1", 2, 1)
  bias_plots_DGP1_R12 <- plot_bias(bias_results_DGP1, "DGP-1", 1, 2)
  bias_plots_DGP1_R22 <- plot_bias(bias_results_DGP1, "DGP-1", 2, 2)
  
  # DGP-2
  bias_plots_DGP2_R11 <- plot_bias(bias_results_DGP2, "DGP-2", 1, 1)
  bias_plots_DGP2_R21 <- plot_bias(bias_results_DGP2, "DGP-2", 2, 1)
  bias_plots_DGP2_R12 <- plot_bias(bias_results_DGP2, "DGP-2", 1, 2)
  bias_plots_DGP2_R22 <- plot_bias(bias_results_DGP2, "DGP-2", 2, 2)
  
  # Save plots
  cat("Saving bias plots...\n")
  if (!is.null(bias_plots_DGP1_R11$plot_small)) ggsave("output/bias_DGP1_R11_small.png", bias_plots_DGP1_R11$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP1_R21$plot_small)) ggsave("output/bias_DGP1_R21_small.png", bias_plots_DGP1_R21$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP1_R12$plot_small)) ggsave("output/bias_DGP1_R12_small.png", bias_plots_DGP1_R12$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP1_R22$plot_small)) ggsave("output/bias_DGP1_R22_small.png", bias_plots_DGP1_R22$plot_small, width = 10, height = 6)
  
  if (!is.null(bias_plots_DGP2_R11$plot_small)) ggsave("output/bias_DGP2_R11_small.png", bias_plots_DGP2_R11$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP2_R21$plot_small)) ggsave("output/bias_DGP2_R21_small.png", bias_plots_DGP2_R21$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP2_R12$plot_small)) ggsave("output/bias_DGP2_R12_small.png", bias_plots_DGP2_R12$plot_small, width = 10, height = 6)
  if (!is.null(bias_plots_DGP2_R22$plot_small)) ggsave("output/bias_DGP2_R22_small.png", bias_plots_DGP2_R22$plot_small, width = 10, height = 6)
  
  # Save bias plots for large samples if they exist
  if (max(T_values) > 400) {
    if (!is.null(bias_plots_DGP1_R11$plot_large)) ggsave("output/bias_DGP1_R11_large.png", bias_plots_DGP1_R11$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP1_R21$plot_large)) ggsave("output/bias_DGP1_R21_large.png", bias_plots_DGP1_R21$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP1_R12$plot_large)) ggsave("output/bias_DGP1_R12_large.png", bias_plots_DGP1_R12$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP1_R22$plot_large)) ggsave("output/bias_DGP1_R22_large.png", bias_plots_DGP1_R22$plot_large, width = 10, height = 6)
    
    if (!is.null(bias_plots_DGP2_R11$plot_large)) ggsave("output/bias_DGP2_R11_large.png", bias_plots_DGP2_R11$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP2_R21$plot_large)) ggsave("output/bias_DGP2_R21_large.png", bias_plots_DGP2_R21$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP2_R12$plot_large)) ggsave("output/bias_DGP2_R12_large.png", bias_plots_DGP2_R12$plot_large, width = 10, height = 6)
    if (!is.null(bias_plots_DGP2_R22$plot_large)) ggsave("output/bias_DGP2_R22_large.png", bias_plots_DGP2_R22$plot_large, width = 10, height = 6)
  }
  
  # Create combined bias figure (similar to Figure 1 in the paper)
  cat("Creating combined bias figure...\n")
  # Select three representative sample sizes if available
  if (length(T_values) >= 3) {
    sample_sizes_to_show <- c(min(T_values), median(T_values), max(T_values))
  } else {
    sample_sizes_to_show <- T_values
  }
  
  combined_figure <- create_combined_bias_figure(bias_results_DGP1, bias_results_DGP2, sample_sizes_to_show)
  if (!is.null(combined_figure)) ggsave("output/combined_bias_figure.png", combined_figure, width = 12, height = 16)
}

# Create RMSE tables
if (run_bias_sim && !is.null(bias_results_DGP1) && !is.null(bias_results_DGP2)) {
  cat("Creating RMSE tables...\n")
  
  # Select representative horizons
  horizons_to_show <- c(1, 3, 5, 10, 15)
  
  # Create tables comparing VAR and LP RMSE values
  rmse_tables <- create_rmse_tables(bias_results_DGP1, bias_results_DGP2, horizons_to_show)
  
  # Format for LaTeX
  latex_tables <- format_tables_for_latex(rmse_tables)
  
  # Save individual tables
  cat("Saving RMSE tables...\n")
  for (name in names(latex_tables)) {
    file_name <- paste0("output/rmse_table_", gsub("[(),]", "", name), ".tex")
    cat(latex_tables[[name]], file = file_name)
  }
  
  # Create combined tables
  combined_table_DGP1 <- create_combined_rmse_table(rmse_tables, "DGP-1")
  combined_table_DGP2 <- create_combined_rmse_table(rmse_tables, "DGP-2")
  
  # Save combined tables
  cat(combined_table_DGP1, file = "output/rmse_table_combined_DGP1.tex")
  cat(combined_table_DGP2, file = "output/rmse_table_combined_DGP2.tex")
  
  cat("RMSE tables created and saved to output directory.\n")
}

# Create DM test rejection plots and tables
if (run_dm_sim && !is.null(dm_results_DGP1) && !is.null(dm_results_DGP2)) {
  cat("Creating DM test rejection plots and tables...\n")
  
  # Create plots
  dm_plots_DGP1 <- plot_dm_rejection(dm_results_DGP1, "DGP-1")
  dm_plots_DGP2 <- plot_dm_rejection(dm_results_DGP2, "DGP-2")
  
  # Save plots
  cat("Saving DM test plots...\n")
  for (T_value in T_values) {
    T_name <- paste0("T", T_value)
    if (!is.null(dm_plots_DGP1$plots) && T_name %in% names(dm_plots_DGP1$plots)) {
      ggsave(paste0("output/dm_DGP1_T", T_value, ".png"), dm_plots_DGP1$plots[[T_name]], 
             width = 10, height = 6)
    }
    if (!is.null(dm_plots_DGP2$plots) && T_name %in% names(dm_plots_DGP2$plots)) {
      ggsave(paste0("output/dm_DGP2_T", T_value, ".png"), dm_plots_DGP2$plots[[T_name]], 
             width = 10, height = 6)
    }
  }
  
  # Create rejection tables (similar to Tables 1 and 2 in the paper)
  cat("Creating rejection tables...\n")
  rejection_tables_DGP1 <- create_rejection_table(dm_results_DGP1, "DGP-1")
  rejection_tables_DGP2 <- create_rejection_table(dm_results_DGP2, "DGP-2")
  
  # Save tables
  cat("Saving rejection tables...\n")
  write.csv(rejection_tables_DGP1$R11, "output/rejection_DGP1_R11.csv")
  write.csv(rejection_tables_DGP1$R21, "output/rejection_DGP1_R21.csv")
  write.csv(rejection_tables_DGP1$R12, "output/rejection_DGP1_R12.csv")
  write.csv(rejection_tables_DGP1$R22, "output/rejection_DGP1_R22.csv")
  
  write.csv(rejection_tables_DGP2$R11, "output/rejection_DGP2_R11.csv")
  write.csv(rejection_tables_DGP2$R21, "output/rejection_DGP2_R21.csv")
  write.csv(rejection_tables_DGP2$R12, "output/rejection_DGP2_R12.csv")
  write.csv(rejection_tables_DGP2$R22, "output/rejection_DGP2_R22.csv")
}

# Save tables as images
cat("Saving rejection tables as images...\n")
save_rejection_tables_as_images(rejection_tables_DGP1, "DGP-1")
save_rejection_tables_as_images(rejection_tables_DGP2, "DGP-2")

cat("Analysis complete. Results saved to output directory.\n")
cat("Project execution completed successfully.\n")