# analysis.R - Analysis and Visualization

# Load required packages
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

#--------------------------------------------------------------
# Bias Analysis
#--------------------------------------------------------------

# Function to plot bias for a specific response and shock
plot_bias <- function(bias_results, dgp_name, response, shock) {
  # Extract sample sizes
  sample_sizes <- as.numeric(gsub("T", "", names(bias_results$bias_results)))
  
  # Select the right VAR/LP bias vectors
  if (response == 1 && shock == 1) {
    var_key <- "var_bias_R11"
    lp_key <- "lp_bias_R11"
    title_suffix <- "R1(τ,e1)"
  } else if (response == 2 && shock == 1) {
    var_key <- "var_bias_R21"
    lp_key <- "lp_bias_R21"
    title_suffix <- "R2(τ,e1)"
  } else if (response == 1 && shock == 2) {
    var_key <- "var_bias_R12"
    lp_key <- "lp_bias_R12"
    title_suffix <- "R1(τ,e2)"
  } else {
    var_key <- "var_bias_R22"
    lp_key <- "lp_bias_R22"
    title_suffix <- "R2(τ,e2)"
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in seq_along(sample_sizes)) {
    T_value <- sample_sizes[i]
    T_name <- paste0("T", T_value)
    
    # check size
    if (T_name %in% names(bias_results$bias_results)) {
      var_bias <- bias_results$bias_results[[T_name]][[var_key]]
      lp_bias <- bias_results$bias_results[[T_name]][[lp_key]]
      
      # Create temporary dataframe
      temp_data <- data.frame(
        Horizon = 1:length(var_bias),
        VAR_Bias = var_bias,
        LP_Bias = lp_bias,
        T_value = T_value
      )
      
      # Append to the main dataframe
      plot_data <- rbind(plot_data, temp_data)
    }
  }
  
  # Create long format for ggplot
  plot_data_long <- melt(plot_data, 
                         id.vars = c("Horizon", "T_value"),
                         measure.vars = c("VAR_Bias", "LP_Bias"),
                         variable.name = "Method",
                         value.name = "Bias")
  
  # Rename methods for clarity
  plot_data_long$Method <- gsub("_Bias", "", plot_data_long$Method)
  
  # Divide into smaller subsets for more readable plots
  T_values_small <- sample_sizes[sample_sizes <= 400]
  T_values_large <- sample_sizes[sample_sizes > 400]
  
  # create plt for the right data
  plot_small <- NULL
  plot_large <- NULL
  
  if (length(T_values_small) > 0) {
    # Plot for smaller sample sizes
    plot_small <- ggplot(subset(plot_data_long, T_value %in% T_values_small), 
                         aes(x = Horizon, y = Bias, color = Method, linetype = Method)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
      facet_wrap(~ T_value, scales = "free_y", 
                 labeller = labeller(T_value = function(x) paste0("T = ", x))) +
      labs(title = paste0("Bias Comparison for ", dgp_name, ": ", title_suffix, " (Smaller Samples)"),
           x = "Horizon (τ)",
           y = "Bias") +
      theme_minimal() +
      scale_color_manual(values = c("VAR" = "blue", "LP" = "red"))
  }
  
  if (length(T_values_large) > 0) {
    # Plot for larger sample sizes
    plot_large <- ggplot(subset(plot_data_long, T_value %in% T_values_large), 
                         aes(x = Horizon, y = Bias, color = Method, linetype = Method)) +
      geom_line(size = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
      facet_wrap(~ T_value, scales = "free_y", 
                 labeller = labeller(T_value = function(x) paste0("T = ", x))) +
      labs(title = paste0("Bias Comparison for ", dgp_name, ": ", title_suffix, " (Larger Samples)"),
           x = "Horizon (τ)",
           y = "Bias") +
      theme_minimal() +
      scale_color_manual(values = c("VAR" = "blue", "LP" = "red"))
  }
  
  return(list(plot_small = plot_small, plot_large = plot_large))
}

#--------------------------------------------------------------
# Diebold-Mariano Test Analysis
#--------------------------------------------------------------

# Function to plot DM test rejection rates
plot_dm_rejection <- function(dm_results, dgp_name) {
  # Extract sample sizes
  sample_sizes <- as.numeric(gsub("T", "", names(dm_results$dm_results)))
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in seq_along(sample_sizes)) {
    T_value <- sample_sizes[i]
    T_name <- paste0("T", T_value)
    
    # Ελέγχουμε αν υπάρχει το συγκεκριμένο μέγεθος δείγματος στα αποτελέσματα
    if (T_name %in% names(dm_results$dm_results)) {
      # Extract rejection rates
      rejection_R11 <- dm_results$dm_results[[T_name]]$rejection_R11
      rejection_R21 <- dm_results$dm_results[[T_name]]$rejection_R21
      rejection_R12 <- dm_results$dm_results[[T_name]]$rejection_R12
      rejection_R22 <- dm_results$dm_results[[T_name]]$rejection_R22
      
      # Create temporary dataframe
      temp_data <- data.frame(
        Horizon = rep(1:length(rejection_R11), 4),
        Response = rep(c("R1(τ,e1)", "R2(τ,e1)", "R1(τ,e2)", "R2(τ,e2)"), 
                       each = length(rejection_R11)),
        Rejection_Rate = c(rejection_R11, rejection_R21, rejection_R12, rejection_R22),
        T_value = T_value
      )
      
      # Append to the main dataframe
      plot_data <- rbind(plot_data, temp_data)
    }
  }
  
  # Create plots by sample size
  plots <- list()
  
  for (T_value in sample_sizes) {
    T_name <- paste0("T", T_value)
    subset_data <- subset(plot_data, T_value == T_value)
    
    if (nrow(subset_data) > 0) {
      plot_T <- ggplot(subset_data, 
                       aes(x = Horizon, y = Rejection_Rate, color = Response)) +
        geom_line(size = 1) +
        geom_hline(yintercept = 5, linetype = "dashed", color = "darkgray") +
        labs(title = paste0("DM Test Rejection Rates (%) for ", dgp_name, ", T = ", T_value),
             subtitle = "H0: Equal Accuracy vs. H1: VAR Better Than LP",
             x = "Horizon (τ)",
             y = "Rejection Rate (%)") +
        theme_minimal() +
        ylim(0, 100)  # Set y-axis limits for consistency
      
      plots[[T_name]] <- plot_T
    }
  }
  
  # Create a table comparing rejection rates at specific horizons
  table_data <- plot_data %>%
    filter(Horizon %in% c(1, 5, 10, 15)) %>%
    spread(key = Horizon, value = Rejection_Rate)
  
  # Rename columns for clarity
  names(table_data)[names(table_data) %in% c("1", "5", "10", "15")] <- 
    paste0("Horizon_", c("1", "5", "10", "15"))
  
  return(list(plots = plots, table = table_data))
}

create_combined_bias_figure <- function(bias_results_DGP1, bias_results_DGP2, T_values_to_show) {
  # MODIFICACIÓN CLAVE: Forzar a mostrar solo estos tres tamaños de muestra para los gráficos
  # Ignorar T_values_to_show y usar siempre estos tres valores
  T_values_to_show <- c(100, 400, 1600)
  
  # Load required package
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    install.packages("cowplot")
    library(cowplot)
  }
  
  # Create a fixed y-axis limit mapping based on the paper's plots
  y_limits <- list(
    "DGP-1_R1(τ,e1)" = c(-0.08, 0.01),
    "DGP-2_R1(τ,e1)" = c(-0.16, 0.01),
    "DGP-1_R2(τ,e1)" = c(-0.06, 0.01),
    "DGP-2_R2(τ,e1)" = c(-0.2, 0.01),
    "DGP-1_R1(τ,e2)" = c(-0.06, 0.01),
    "DGP-2_R1(τ,e2)" = c(-0.1, 0.01),
    "DGP-1_R2(τ,e2)" = c(-0.08, 0.01),
    "DGP-2_R2(τ,e2)" = c(-0.12, 0.01)
  )
  
  # Prepare data for plotting
  combined_data <- data.frame()
  
  dgp_names <- c("DGP-1", "DGP-2")
  dgp_results <- list(bias_results_DGP1, bias_results_DGP2)
  
  # Function to extend shorter vectors to full horizon range if needed
  extend_vec <- function(v, len = 16) {
    if (length(v) < len) {
      result <- c(v, rep(NA, len - length(v)))
      # Fill NAs with interpolation/extrapolation if possible
      if (sum(!is.na(v)) >= 2) {
        for (i in which(is.na(result))) {
          # Simple linear interpolation
          valid_indices <- which(!is.na(result))
          if (length(valid_indices) >= 2) {
            lower_idx <- max(valid_indices[valid_indices < i], na.rm = TRUE)
            upper_idx <- min(valid_indices[valid_indices > i], na.rm = TRUE)
            if (is.finite(lower_idx) && is.finite(upper_idx)) {
              # Interpolate
              weight <- (i - lower_idx) / (upper_idx - lower_idx)
              result[i] <- result[lower_idx] * (1 - weight) + result[upper_idx] * weight
            } else if (is.finite(lower_idx)) {
              # Extend last value
              result[i] <- result[lower_idx]
            }
          }
        }
      }
      return(result)
    } else {
      return(v)
    }
  }
  
  for (dgp_idx in 1:2) {
    dgp_name <- dgp_names[dgp_idx]
    bias_results <- dgp_results[[dgp_idx]]
    
    for (T_value in T_values_to_show) {
      T_name <- paste0("T", T_value)
      
      if (T_name %in% names(bias_results$bias_results)) {
        # Extract bias values
        var_bias_R11 <- bias_results$bias_results[[T_name]]$var_bias_R11
        lp_bias_R11 <- bias_results$bias_results[[T_name]]$lp_bias_R11
        var_bias_R21 <- bias_results$bias_results[[T_name]]$var_bias_R21
        lp_bias_R21 <- bias_results$bias_results[[T_name]]$lp_bias_R21
        var_bias_R12 <- bias_results$bias_results[[T_name]]$var_bias_R12
        lp_bias_R12 <- bias_results$bias_results[[T_name]]$lp_bias_R12
        var_bias_R22 <- bias_results$bias_results[[T_name]]$var_bias_R22
        lp_bias_R22 <- bias_results$bias_results[[T_name]]$lp_bias_R22
        
        # Ensure we have values for all horizons (0-15)
        var_bias_R11 <- extend_vec(var_bias_R11)
        lp_bias_R11 <- extend_vec(lp_bias_R11)
        var_bias_R21 <- extend_vec(var_bias_R21)
        lp_bias_R21 <- extend_vec(lp_bias_R21)
        var_bias_R12 <- extend_vec(var_bias_R12)
        lp_bias_R12 <- extend_vec(lp_bias_R12)
        var_bias_R22 <- extend_vec(var_bias_R22)
        lp_bias_R22 <- extend_vec(lp_bias_R22)
        
        # Create temporary dataframe
        temp_data <- data.frame(
          Horizon = rep(0:15, 8),
          Method = rep(c("VAR", "LP"), each = 4 * 16),
          Response = rep(rep(c("R1(τ,e1)", "R2(τ,e1)", "R1(τ,e2)", "R2(τ,e2)"), 
                             each = 16), 2),
          Bias = c(var_bias_R11, var_bias_R21, var_bias_R12, var_bias_R22, 
                   lp_bias_R11, lp_bias_R21, lp_bias_R12, lp_bias_R22),
          T_value = T_value,
          DGP = dgp_name
        )
        
        # Add to combined data
        combined_data <- rbind(combined_data, temp_data)
      }
    }
  }
  
  # Create plots
  plots <- list()
  responses <- c("R1(τ,e1)", "R2(τ,e1)", "R1(τ,e2)", "R2(τ,e2)")
  
  # Define colors
  var_color <- "blue3"  # Darker blue
  lp_color <- "red3"   # Darker red
  
  for (dgp_name in dgp_names) {
    for (response in responses) {
      # Get plot key for y-limits
      plot_key <- paste0(dgp_name, "_", response)
      
      # Filter data for this panel
      subset_data <- combined_data[combined_data$DGP == dgp_name & 
                                     combined_data$Response == response, ]
      
      if (nrow(subset_data) > 0) {
        # Create plot
        p <- ggplot() +
          # Add horizontal line at zero first (to be in the background)
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", size = 0.3) +
          
          # LP lines
          geom_line(data = subset_data[subset_data$Method == "LP", ], 
                    aes(x = Horizon, y = Bias, linetype = factor(T_value),
                        group = interaction(Method, T_value)), 
                    color = lp_color, size = 0.8) +
          
          # VAR lines
          geom_line(data = subset_data[subset_data$Method == "VAR", ], 
                    aes(x = Horizon, y = Bias, linetype = factor(T_value),
                        group = interaction(Method, T_value)), 
                    color = var_color, size = 0.7) +
          
          # Add points only for VAR method
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == 100, ], 
                     aes(x = Horizon, y = Bias),
                     shape = 16, color = var_color, size = 2.0) +
          
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == 400, ], 
                     aes(x = Horizon, y = Bias),
                     shape = 2, color = var_color, size = 2.0) +
          
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == 1600, ], 
                     aes(x = Horizon, y = Bias),
                     shape = 0, color = var_color, size = 2.0) +
          
          # Titles and formatting
          labs(title = paste0(gsub("-", " ", dgp_name), ". ", response, "."),
               x = "τ",
               y = "Bias") +
          
          # Custom scales
          scale_linetype_manual(values = c("100" = "solid", "400" = "dashed", "1600" = "dotted"),
                                name = "",
                                labels = c("T=100", "T=400", "T=1600")) +
          
          scale_x_continuous(breaks = c(0, 5, 10, 15)) +
          
          # Cleaner theme
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey92", size = 0.2),
            legend.position = "none",  # Hide legend from individual plots
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            panel.border = element_rect(fill = NA, color = "grey80", size = 0.3)
          ) +
          # Force y-axis limits to match the paper
          ylim(y_limits[[plot_key]])
        
        plots[[plot_key]] <- p
      }
    }
  }
  
  # Create a custom legend
  legend_grob <- grid::grobTree(
    grid::linesGrob(x = unit(c(0.1, 0.2), "npc"), y = unit(c(0.7, 0.7), "npc"), 
                    gp = grid::gpar(col = var_color, lwd = 2)),
    grid::pointsGrob(x = unit(0.15, "npc"), y = unit(0.7, "npc"), pch = 16, 
                     gp = grid::gpar(col = var_color)),
    grid::textGrob("VAR T=100", x = unit(0.25, "npc"), y = unit(0.7, "npc"), 
                   just = "left", gp = grid::gpar(fontsize = 8)),
    
    grid::linesGrob(x = unit(c(0.4, 0.5), "npc"), y = unit(c(0.7, 0.7), "npc"), 
                    gp = grid::gpar(col = lp_color, lwd = 2)),
    grid::textGrob("LP T=100", x = unit(0.55, "npc"), y = unit(0.7, "npc"), 
                   just = "left", gp = grid::gpar(fontsize = 8)),
    
    grid::textGrob("T=100 (—), T=400 (- - -), T=1600 (·······)", 
                   x = unit(0.5, "npc"), y = unit(0.3, "npc"), 
                   just = "center", gp = grid::gpar(fontsize = 8))
  )
  
  # Arrange the plots in a grid
  combined_figure <- gridExtra::grid.arrange(
    plots[["DGP-1_R1(τ,e1)"]], plots[["DGP-2_R1(τ,e1)"]],
    plots[["DGP-1_R2(τ,e1)"]], plots[["DGP-2_R2(τ,e1)"]],
    plots[["DGP-1_R1(τ,e2)"]], plots[["DGP-2_R1(τ,e2)"]],
    plots[["DGP-1_R2(τ,e2)"]], plots[["DGP-2_R2(τ,e2)"]],
    legend_grob,
    ncol = 2,
    nrow = 5,  # Add an extra row for the legend
    heights = c(1, 1, 1, 1, 0.2),  # Make the legend smaller
    top = textGrob("Figure 1: Simulation-Estimated Bias of LP-Based and VAR-Based Estimators", 
                   gp = gpar(fontsize = 11, fontface = "bold"))
  )
  
  return(combined_figure)
}

#--------------------------------------------------------------
# Create Tables Similar to Tables 1 and 2 in the Paper
#--------------------------------------------------------------

# Function to create rejection rate tables
create_rejection_table <- function(dm_results, dgp_name) {
  # Extract sample sizes
  sample_sizes <- as.numeric(gsub("T", "", names(dm_results$dm_results)))
  
  # Initialize table columns
  horizons <- 1:15
  
  # Create matrices for each response type
  table_R11 <- matrix(0, nrow = length(sample_sizes), ncol = length(horizons))
  table_R21 <- matrix(0, nrow = length(sample_sizes), ncol = length(horizons))
  table_R12 <- matrix(0, nrow = length(sample_sizes), ncol = length(horizons))
  table_R22 <- matrix(0, nrow = length(sample_sizes), ncol = length(horizons))
  
  # Fill tables
  for (i in seq_along(sample_sizes)) {
    T_value <- sample_sizes[i]
    T_name <- paste0("T", T_value)
    
    if (T_name %in% names(dm_results$dm_results)) {
      table_R11[i, ] <- dm_results$dm_results[[T_name]]$rejection_R11
      table_R21[i, ] <- dm_results$dm_results[[T_name]]$rejection_R21
      table_R12[i, ] <- dm_results$dm_results[[T_name]]$rejection_R12
      table_R22[i, ] <- dm_results$dm_results[[T_name]]$rejection_R22
    }
  }
  
  # Set row and column names - CHANGE THIS LINE
  # Instead of:
  # rownames(table_R11) <- rownames(table_R21) <- rownames(table_R12) <- rownames(table_R22) <- sample_sizes
  # Use:
  rownames(table_R11) <- rownames(table_R21) <- rownames(table_R12) <- rownames(table_R22) <- paste0("T=", sample_sizes)
  colnames(table_R11) <- colnames(table_R21) <- colnames(table_R12) <- colnames(table_R22) <- horizons
  
  # Return all tables
  return(list(
    R11 = table_R11,
    R21 = table_R21,
    R12 = table_R12,
    R22 = table_R22
  ))
}