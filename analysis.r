# analysis.R - Analysis and Visualization

# Load required packages
library(ggplot2)
library(gridExtra)
library(reshape2)
library(tidyr)

#--------------------------------------------------------------
# Bias Analysis
#--------------------------------------------------------------

# Function to create RMSE comparison tables
create_rmse_tables <- function(results_DGP1, results_DGP2, horizons_to_show = c(1, 5, 10, 15)) {
  # Define response-shock combinations
  combinations <- list(
    "R1(τ,e1)" = c(1, 1),
    "R2(τ,e1)" = c(2, 1),
    "R1(τ,e2)" = c(1, 2),
    "R2(τ,e2)" = c(2, 2)
  )
  
  # Check if rmse_results exists in both results
  if (!("rmse_results" %in% names(results_DGP1)) || 
      !("rmse_results" %in% names(results_DGP2)) ||
      length(results_DGP1$rmse_results) == 0 || 
      length(results_DGP2$rmse_results) == 0) {
    warning("RMSE results not found in one or both DGP results")
    return(list())
  }
  
  # Sample sizes - try to extract from both and use the common ones
  sample_sizes_DGP1 <- as.numeric(gsub("T", "", names(results_DGP1$rmse_results)))
  sample_sizes_DGP2 <- as.numeric(gsub("T", "", names(results_DGP2$rmse_results)))
  sample_sizes <- sort(intersect(sample_sizes_DGP1, sample_sizes_DGP2))
  
  
  # Check if we have any common sample sizes
  if (length(sample_sizes) == 0) {
    warning("No common sample sizes found in RMSE results")
    return(list())
  }
  
  # Create tables for each DGP
  tables <- list()
  
  for (dgp_idx in 1:2) {
    dgp_name <- paste0("DGP-", dgp_idx)
    results <- if (dgp_idx == 1) results_DGP1 else results_DGP2
    
    for (combo_name in names(combinations)) {
      response <- combinations[[combo_name]][1]
      shock <- combinations[[combo_name]][2]
      
      # Determine the right RMSE keys
      if (response == 1 && shock == 1) {
        var_key <- "var_rmse_R11"
        lp_key <- "lp_rmse_R11"
      } else if (response == 2 && shock == 1) {
        var_key <- "var_rmse_R21"
        lp_key <- "lp_rmse_R21"
      } else if (response == 1 && shock == 2) {
        var_key <- "var_rmse_R12"
        lp_key <- "lp_rmse_R12"
      } else {
        var_key <- "var_rmse_R22"
        lp_key <- "lp_rmse_R22"
      }
      
      # Initialize table
      table_data <- matrix(NA, nrow = length(sample_sizes), ncol = 2 * length(horizons_to_show))
      rownames(table_data) <- paste0("T=", sample_sizes)
      
      # Create column names
      col_names <- c()
      for (h in horizons_to_show) {
        col_names <- c(col_names, paste0("VAR τ=", h), paste0("LP τ=", h))
      }
      colnames(table_data) <- col_names
      
      # Fill the table
      valid_data_found <- FALSE
      for (i in seq_along(sample_sizes)) {
        T_value <- sample_sizes[i]
        T_name <- paste0("T", T_value)
        
        if (T_name %in% names(results$rmse_results)) {
          rm_results <- results$rmse_results[[T_name]]
          
          # Check if both keys exist
          if (var_key %in% names(rm_results) && lp_key %in% names(rm_results)) {
            var_rmse <- rm_results[[var_key]]
            lp_rmse <- rm_results[[lp_key]]
            
            # Check for valid data
            if (length(var_rmse) > 0 && length(lp_rmse) > 0) {
              valid_data_found <- TRUE
              
              # For each selected horizon
              for (j in seq_along(horizons_to_show)) {
                h <- horizons_to_show[j]
                if (h < length(var_rmse) && h < length(lp_rmse)) {
                  # VAR RMSE at this horizon
                  table_data[i, 2*j-1] <- var_rmse[h+1]  # +1 because horizon 0 is at index 1
                  # LP RMSE at this horizon
                  table_data[i, 2*j] <- lp_rmse[h+1]
                }
              }
            }
          }
        }
      }
      
      # Only store the table if we found valid data
      if (valid_data_found) {
        tables[[paste0(dgp_name, "_", combo_name)]] <- table_data
      } else {
        # Create a dummy table with a message
        dummy_table <- matrix("No data", nrow = 1, ncol = 1)
        tables[[paste0(dgp_name, "_", combo_name)]] <- dummy_table
      }
    }
  }
  
  return(tables)
}

# Create LaTeX-formatted tables
format_tables_for_latex <- function(tables) {
  latex_tables <- list()
  
  if (length(tables) == 0) {
    # Return a dummy table if input is empty
    dummy_latex <- "\\begin{table}[htbp]\n\\centering\n\\caption{RMSE Comparison (No Data Available)}\n\\label{tab:rmse_dummy}\n\\begin{tabular}{c}\n\\hline\nNo RMSE data available. Please run simulations with RMSE calculations enabled. \\\\\n\\hline\n\\end{tabular}\n\\end{table}\n"
    latex_tables[["dummy"]] <- dummy_latex
    return(latex_tables)
  }
  
  for (name in names(tables)) {
    table_data <- tables[[name]]
    
    # Check if it's a dummy table
    if (identical(table_data, matrix("No data", nrow = 1, ncol = 1))) {
      latex <- "\\begin{table}[htbp]\n"
      latex <- paste0(latex, "\\centering\n")
      latex <- paste0(latex, "\\caption{RMSE Comparison for ", name, " (No Data Available)}\n")
      latex <- paste0(latex, "\\label{tab:rmse_", gsub("[(),]", "", name), "}\n")
      latex <- paste0(latex, "\\begin{tabular}{c}\n\\hline\nNo RMSE data available for this combination. \\\\\n\\hline\n\\end{tabular}\n\\end{table}\n")
      
      latex_tables[[name]] <- latex
      next
    }
    
    # For regular tables
    latex <- "\\begin{table}[htbp]\n"
    latex <- paste0(latex, "\\centering\n")
    latex <- paste0(latex, "\\caption{RMSE Comparison for ", name, "}\n")
    latex <- paste0(latex, "\\label{tab:rmse_", gsub("[(),]", "", name), "}\n")
    latex <- paste0(latex, "\\begin{tabular}{l", paste(rep("cc", ncol(table_data)/2), collapse=""), "}\n")
    latex <- paste0(latex, "\\hline\n")
    
    # Column headers
    latex <- paste0(latex, "Sample Size & ")
    col_headers <- colnames(table_data)
    for (i in seq(1, length(col_headers), by=2)) {
      if (i > 1) latex <- paste0(latex, " & ")
      horizon <- gsub("VAR τ=", "", col_headers[i])
      latex <- paste0(latex, "\\multicolumn{2}{c}{τ=", horizon, "}")
    }
    latex <- paste0(latex, " \\\\\n")
    
    # Subheaders
    latex <- paste0(latex, " & ")
    for (i in 1:(ncol(table_data)/2)) {
      if (i > 1) latex <- paste0(latex, " & ")
      latex <- paste0(latex, "VAR & LP")
    }
    latex <- paste0(latex, " \\\\\n")
    latex <- paste0(latex, "\\hline\n")
    
    # Data rows
    for (i in 1:nrow(table_data)) {
      latex <- paste0(latex, rownames(table_data)[i], " & ")
      for (j in 1:ncol(table_data)) {
        if (j > 1) latex <- paste0(latex, " & ")
        value <- table_data[i, j]
        if (!is.na(value)) {
          if (is.numeric(value)) {
            latex <- paste0(latex, sprintf("%.3f", value))
          } else {
            latex <- paste0(latex, as.character(value))
          }
        } else {
          latex <- paste0(latex, "---")
        }
      }
      latex <- paste0(latex, " \\\\\n")
    }
    
    # End LaTeX table
    latex <- paste0(latex, "\\hline\n")
    latex <- paste0(latex, "\\end{tabular}\n")
    latex <- paste0(latex, "\\end{table}\n")
    
    latex_tables[[name]] <- latex
  }
  
  return(latex_tables)
}

# Create a combined table for all response-shock combinations for each DGP
create_combined_rmse_table <- function(tables, dgp) {
  # Check if we have any tables for this DGP
  dgp_tables <- tables[grep(paste0(dgp, "_"), names(tables))]
  
  if (length(dgp_tables) == 0) {
    # No tables found for this DGP
    warning(paste0("No RMSE tables found for ", dgp))
    
    # Return a dummy table with a message
    return(paste0("\\begin{table}[htbp]
\\centering
\\caption{RMSE Comparison for ", dgp, " (No Data Available)}
\\label{tab:rmse_combined_", gsub("-", "", dgp), "}
\\begin{tabular}{c}
\\hline
No RMSE data available for this DGP \\\\
\\hline
\\end{tabular}
\\end{table}"))
  }
  
  # Find a valid table to extract dimensions
  valid_table_name <- NULL
  for (name in names(dgp_tables)) {
    if (is.matrix(dgp_tables[[name]]) && !identical(dgp_tables[[name]], matrix("No data", 1, 1))) {
      valid_table_name <- name
      break
    }
  }
  
  if (is.null(valid_table_name)) {
    # No valid tables found
    warning(paste0("No valid RMSE tables found for ", dgp))
    
    # Return a dummy table with a message
    return(paste0("\\begin{table}[htbp]
\\centering
\\caption{RMSE Comparison for ", dgp, " (No Valid Data Available)}
\\label{tab:rmse_combined_", gsub("-", "", dgp), "}
\\begin{tabular}{c}
\\hline
No valid RMSE data available for this DGP \\\\
\\hline
\\end{tabular}
\\end{table}"))
  }
  
  # Extract sample sizes and horizons
  table_name <- valid_table_name
  sample_sizes <- rownames(dgp_tables[[table_name]])
  column_headers <- colnames(dgp_tables[[table_name]])
  horizons <- unique(gsub("VAR τ=|LP τ=", "", column_headers))
  
  # Create LaTeX table
  latex <- "\\begin{table}[htbp]\n"
  latex <- paste0(latex, "\\centering\n")
  latex <- paste0(latex, "\\caption{RMSE Comparison for ", dgp, "}\n")
  latex <- paste0(latex, "\\label{tab:rmse_combined_", gsub("-", "", dgp), "}\n")
  latex <- paste0(latex, "\\begin{tabular}{ll", paste(rep("cc", length(horizons)), collapse=""), "}\n")
  latex <- paste0(latex, "\\hline\n")
  
  # Column headers
  latex <- paste0(latex, "Response & Sample Size & ")
  for (i in seq_along(horizons)) {
    if (i > 1) latex <- paste0(latex, " & ")
    latex <- paste0(latex, "\\multicolumn{2}{c}{τ=", horizons[i], "}")
  }
  latex <- paste0(latex, " \\\\\n")
  
  # Subheaders
  latex <- paste0(latex, " & & ")
  for (i in 1:length(horizons)) {
    if (i > 1) latex <- paste0(latex, " & ")
    latex <- paste0(latex, "VAR & LP")
  }
  latex <- paste0(latex, " \\\\\n")
  latex <- paste0(latex, "\\hline\n")
  
  # Data rows
  responses <- c("R1(τ,e1)", "R2(τ,e1)", "R1(τ,e2)", "R2(τ,e2)")
  for (response in responses) {
    response_key <- paste0(dgp, "_", response)
    
    if (response_key %in% names(dgp_tables) && 
        is.matrix(dgp_tables[[response_key]]) && 
        !identical(dgp_tables[[response_key]], matrix("No data", 1, 1))) {
      
      table_data <- dgp_tables[[response_key]]
      
      # Add multirow for response
      latex <- paste0(latex, "\\multirow{", length(sample_sizes), "}{*}{", response, "} & ")
      
      for (i in 1:length(sample_sizes)) {
        if (i > 1) latex <- paste0(latex, " & ")
        latex <- paste0(latex, sample_sizes[i], " & ")
        
        for (h_idx in seq_along(horizons)) {
          h <- as.numeric(horizons[h_idx])
          var_col <- 2*h_idx - 1
          lp_col <- 2*h_idx
          
          if (var_col <= ncol(table_data) && lp_col <= ncol(table_data)) {
            var_value <- table_data[i, var_col]
            lp_value <- table_data[i, lp_col]
            
            if (!is.na(var_value)) {
              if (is.numeric(var_value)) {
                latex <- paste0(latex, sprintf("%.3f", var_value))
              } else {
                latex <- paste0(latex, as.character(var_value))
              }
            } else {
              latex <- paste0(latex, "---")
            }
            
            latex <- paste0(latex, " & ")
            
            if (!is.na(lp_value)) {
              if (is.numeric(lp_value)) {
                latex <- paste0(latex, sprintf("%.3f", lp_value))
              } else {
                latex <- paste0(latex, as.character(lp_value))
              }
            } else {
              latex <- paste0(latex, "---")
            }
            
            if (h_idx < length(horizons)) latex <- paste0(latex, " & ")
          }
        }
        
        latex <- paste0(latex, " \\\\\n")
        if (i < length(sample_sizes)) latex <- paste0(latex, " & ")
      }
    } else {
      # No data for this response
      latex <- paste0(latex, response, " & \\multicolumn{", 1 + 2*length(horizons), "}{l}{No data available} \\\\\n")
    }
    
    if (response != responses[length(responses)]) {
      latex <- paste0(latex, "\\hline\n")
    }
  }
  
  # End LaTeX table
  latex <- paste0(latex, "\\hline\n")
  latex <- paste0(latex, "\\end{tabular}\n")
  latex <- paste0(latex, "\\end{table}\n")
  
  return(latex)
}

# Function to create combined RMSE figure
create_combined_rmse_figure <- function(results_DGP1, results_DGP2, T_values_to_show) {
  # Check if RMSE results exist
  if (!("rmse_results" %in% names(results_DGP1)) || 
      !("rmse_results" %in% names(results_DGP2)) ||
      length(results_DGP1$rmse_results) == 0 || 
      length(results_DGP2$rmse_results) == 0) {
    warning("RMSE results not found in one or both DGP results")
    return(NULL)
  }
  
  # Load required package
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    install.packages("cowplot")
    library(cowplot)
  }
  
  # Create a fixed y-axis limit mapping based on the data
  y_limits <- list(
    "DGP-1_R1(τ,e1)" = c(0, 0.2),
    "DGP-2_R1(τ,e1)" = c(0, 0.3),
    "DGP-1_R2(τ,e1)" = c(0, 0.2),
    "DGP-2_R2(τ,e1)" = c(0, 0.3),
    "DGP-1_R1(τ,e2)" = c(0, 0.2),
    "DGP-2_R1(τ,e2)" = c(0, 0.2),
    "DGP-1_R2(τ,e2)" = c(0, 0.2),
    "DGP-2_R2(τ,e2)" = c(0, 0.2)
  )
  
  # Prepare data for plotting
  combined_data <- data.frame()
  
  dgp_names <- c("DGP-1", "DGP-2")
  dgp_results <- list(results_DGP1, results_DGP2)
  
  for (dgp_idx in 1:2) {
    dgp_name <- dgp_names[dgp_idx]
    results <- dgp_results[[dgp_idx]]
    
    for (T_value in T_values_to_show) {
      T_name <- paste0("T", T_value)
      
      if (T_name %in% names(results$rmse_results)) {
        # Extract RMSE values
        rmse_values <- results$rmse_results[[T_name]]
        
        # Check if we have all required RMSE values
        required_keys <- c("var_rmse_R11", "lp_rmse_R11", "var_rmse_R21", "lp_rmse_R21",
                           "var_rmse_R12", "lp_rmse_R12", "var_rmse_R22", "lp_rmse_R22")
        has_all_keys <- all(required_keys %in% names(rmse_values))
        
        if (has_all_keys) {
          var_rmse_R11 <- rmse_values$var_rmse_R11
          lp_rmse_R11 <- rmse_values$lp_rmse_R11
          var_rmse_R21 <- rmse_values$var_rmse_R21
          lp_rmse_R21 <- rmse_values$lp_rmse_R21
          var_rmse_R12 <- rmse_values$var_rmse_R12
          lp_rmse_R12 <- rmse_values$lp_rmse_R12
          var_rmse_R22 <- rmse_values$var_rmse_R22
          lp_rmse_R22 <- rmse_values$lp_rmse_R22
          
          # Check lengths
          if (length(var_rmse_R11) > 0 && length(lp_rmse_R11) > 0) {
            # Create temporary dataframe
            temp_data <- data.frame(
              Horizon = rep(0:(length(var_rmse_R11)-1), 8),
              Method = rep(c("VAR", "LP"), each = 4 * length(var_rmse_R11)),
              Response = rep(rep(c("R1(τ,e1)", "R2(τ,e1)", "R1(τ,e2)", "R2(τ,e2)"), 
                                 each = length(var_rmse_R11)), 2),
              RMSE = c(var_rmse_R11, var_rmse_R21, var_rmse_R12, var_rmse_R22, 
                       lp_rmse_R11, lp_rmse_R21, lp_rmse_R12, lp_rmse_R22),
              T_value = T_value,
              DGP = dgp_name
            )
            
            # Add to combined data
            combined_data <- rbind(combined_data, temp_data)
          }
        }
      }
    }
  }
  
  # Check if we have data to plot
  if (nrow(combined_data) == 0) {
    warning("No valid RMSE data for plotting")
    return(NULL)
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
          # Add horizontal grid lines
          geom_hline(yintercept = seq(0, max(subset_data$RMSE, na.rm = TRUE) * 1.1, 
                                      length.out = 5), 
                     linetype = "dotted", color = "gray80", size = 0.3) +
          
          # LP lines
          geom_line(data = subset_data[subset_data$Method == "LP", ], 
                    aes(x = Horizon, y = RMSE, linetype = factor(T_value),
                        group = interaction(Method, T_value)), 
                    color = lp_color, size = 0.8) +
          
          # VAR lines
          geom_line(data = subset_data[subset_data$Method == "VAR", ], 
                    aes(x = Horizon, y = RMSE, linetype = factor(T_value),
                        group = interaction(Method, T_value)), 
                    color = var_color, size = 0.7) +
          
          # Add points for identification
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == min(T_values_to_show), ], 
                     aes(x = Horizon, y = RMSE),
                     shape = 16, color = var_color, size = 2.0) +
          
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == median(T_values_to_show), ], 
                     aes(x = Horizon, y = RMSE),
                     shape = 2, color = var_color, size = 2.0) +
          
          geom_point(data = subset_data[subset_data$Method == "VAR" & subset_data$T_value == max(T_values_to_show), ], 
                     aes(x = Horizon, y = RMSE),
                     shape = 0, color = var_color, size = 2.0) +
          
          # Titles and formatting
          labs(title = paste0(gsub("-", " ", dgp_name), ". ", response, "."),
               x = "τ",
               y = "RMSE") +
          
          # Custom scales
          scale_linetype_manual(values = c(
            "100" = "solid", 
            "200" = "longdash",
            "400" = "dashed", 
            "800" = "dotdash",
            "1600" = "dotted"
          ),
          name = "",
          labels = paste0("T=", T_values_to_show)) +
          
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
          # Force y-axis limits if available
          ylim(c(0, max(subset_data$RMSE, na.rm = TRUE) * 1.1))
        
        plots[[plot_key]] <- p
      } else {
        # Create a placeholder plot for missing data
        p <- ggplot() +
          annotate("text", x = 7.5, y = 0.5, label = "No RMSE data available") +
          labs(title = paste0(gsub("-", " ", dgp_name), ". ", response, "."),
               x = "τ", y = "RMSE") +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey92", size = 0.2),
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            panel.border = element_rect(fill = NA, color = "grey80", size = 0.3)
          ) +
          xlim(0, 15) + ylim(0, 1)
        
        plots[[plot_key]] <- p
      }
    }
  }
  
  # Check if we have all needed plots
  if (length(plots) < 8) {
    warning("Not all required plots could be created")
    return(NULL)
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
    top = grid::textGrob("RMSE of LP-Based and VAR-Based Estimators", 
                         gp = grid::gpar(fontsize = 11, fontface = "bold"))
  )
  
  return(combined_figure)
}

# Function to save rejection tables as images
save_rejection_tables_as_images <- function(rejection_tables, dgp_name) {
  if (is.null(rejection_tables)) {
    warning("No rejection tables provided for ", dgp_name)
    return(NULL)
  }
  
  library(ggplot2)
  library(reshape2)
  
  # Function to create a table plot
  create_table_plot <- function(data, table_number, response) {
    # Prepare data for plotting
    # Add row names as a column
    data_plot <- as.data.frame(data)
    data_plot$Sample <- rownames(data_plot)
    
    # Melt data for ggplot (convert to long format)
    melted_data <- reshape2::melt(data_plot, id.vars = "Sample", 
                                  variable.name = "Horizon", value.name = "Value")
    
    # Convert Horizon from factor to numeric for proper ordering
    melted_data$Horizon <- as.numeric(as.character(melted_data$Horizon))
    
    # Create the plot
    p <- ggplot(melted_data, aes(x = Horizon, y = Sample)) +
      # Add colored cells
      geom_tile(fill = "#f0f0f0", color = "gray90") +
      # Add text values
      geom_text(aes(label = sprintf("%.1f", Value)), size = 3) +
      # Add column labels (Horizons)
      scale_x_continuous(breaks = 1:15, 
                         labels = paste0("τ=", 1:15),
                         expand = c(0, 0)) +
      # Formatting
      labs(
        title = paste0("Table ", table_number, ": Percentage Rejections for ", dgp_name, ", ", response),
        x = "Horizon",
        y = ""
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
        axis.text.y = element_text(face = "bold", size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
      )
    
    return(p)
  }
  
  # Check if key tables exist
  if (!all(c("R11", "R21", "R12", "R22") %in% names(rejection_tables))) {
    warning("Missing tables for ", dgp_name)
    return(NULL)
  }
  
  # Create each table
  table_R11 <- create_table_plot(
    rejection_tables$R11, 
    ifelse(dgp_name == "DGP-1", "1a", "2a"), 
    "R1(τ,e1)"
  )
  
  table_R21 <- create_table_plot(
    rejection_tables$R21, 
    ifelse(dgp_name == "DGP-1", "1b", "2b"),
    "R2(τ,e1)"
  )
  
  table_R12 <- create_table_plot(
    rejection_tables$R12, 
    ifelse(dgp_name == "DGP-1", "1c", "2c"),
    "R1(τ,e2)"
  )
  
  table_R22 <- create_table_plot(
    rejection_tables$R22, 
    ifelse(dgp_name == "DGP-1", "1d", "2d"),
    "R2(τ,e2)"
  )
  
  # Save individual tables
  ggsave(paste0("output/table_", dgp_name, "_R11.png"), table_R11, width = 10, height = 3)
  ggsave(paste0("output/table_", dgp_name, "_R21.png"), table_R21, width = 10, height = 3)
  ggsave(paste0("output/table_", dgp_name, "_R12.png"), table_R12, width = 10, height = 3)
  ggsave(paste0("output/table_", dgp_name, "_R22.png"), table_R22, width = 10, height = 3)
  
  # Create combined plot
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    
    # Add main title
    title <- grid::textGrob(paste0(dgp_name, " Rejection Rates"), 
                            gp = grid::gpar(fontsize = 14, fontface = "bold"))
    
    # Arrange plots
    combined <- gridExtra::grid.arrange(
      title,
      table_R11, table_R21, table_R12, table_R22,
      ncol = 1,
      heights = c(0.5, 3, 3, 3, 3)
    )
    
    # Save combined image
    ggsave(paste0("output/tables_combined_", dgp_name, ".png"), combined, width = 10, height = 14)
  }
  
  cat("Saved rejection tables as images for", dgp_name, "\n")
  return(TRUE)
}

# Function to plot bias for a specific response and shock
plot_bias <- function(bias_results, dgp_name, response, shock) {
  # Check if bias_results has required structure
  if (!("bias_results" %in% names(bias_results))) {
    warning("Invalid bias_results structure for ", dgp_name)
    return(list(plot_small = NULL, plot_large = NULL))
  }
  
  # Extract sample sizes
  sample_sizes <- as.numeric(gsub("T", "", names(bias_results$bias_results)))
  
  if (length(sample_sizes) == 0) {
    warning("No sample sizes found in bias_results for ", dgp_name)
    return(list(plot_small = NULL, plot_large = NULL))
  }
  
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
      bias_res <- bias_results$bias_results[[T_name]]
      
      if (var_key %in% names(bias_res) && lp_key %in% names(bias_res)) {
        var_bias <- bias_res[[var_key]]
        lp_bias <- bias_res[[lp_key]]
        
        if (length(var_bias) > 0 && length(lp_bias) > 0) {
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
    }
  }
  
  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    warning("No valid bias data for ", dgp_name, ", response=", response, ", shock=", shock)
    return(list(plot_small = NULL, plot_large = NULL))
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
  # Check if dm_results has required structure
  if (!("dm_results" %in% names(dm_results))) {
    warning("Invalid dm_results structure for ", dgp_name)
    return(list(plots = list(), table = NULL))
  }
  
  # Extract sample sizes
  sample_sizes <- as.numeric(gsub("T", "", names(dm_results$dm_results)))
  
  if (length(sample_sizes) == 0) {
    warning("No sample sizes found in dm_results for ", dgp_name)
    return(list(plots = list(), table = NULL))
  }
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (i in seq_along(sample_sizes)) {
    T_value <- sample_sizes[i]
    T_name <- paste0("T", T_value)
    
    # Check if this sample size exists in the results
    if (T_name %in% names(dm_results$dm_results)) {
      dm_res <- dm_results$dm_results[[T_name]]
      
      # Check if all required keys exist
      required_keys <- c("rejection_R11", "rejection_R21", "rejection_R12", "rejection_R22")
      
      if (all(required_keys %in% names(dm_res))) {
        # Extract rejection rates
        rejection_R11 <- dm_res$rejection_R11
        rejection_R21 <- dm_res$rejection_R21
        rejection_R12 <- dm_res$rejection_R12
        rejection_R22 <- dm_res$rejection_R22
        
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
  }
  
  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    warning("No valid rejection data for ", dgp_name)
    return(list(plots = list(), table = NULL))
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
  if (nrow(plot_data) > 0) {
    table_data <- plot_data %>%
      filter(Horizon %in% c(1, 5, 10, 15)) %>%
      spread(key = Horizon, value = Rejection_Rate)
    
    # Rename columns for clarity
    names(table_data)[names(table_data) %in% c("1", "5", "10", "15")] <- 
      paste0("Horizon_", c("1", "5", "10", "15"))
  } else {
    table_data <- NULL
  }
  
  return(list(plots = plots, table = table_data))
}

# Function to create combined bias figure
create_combined_bias_figure <- function(bias_results_DGP1, bias_results_DGP2, T_values_to_show) {
  # Force to show only these three sample sizes for the plots
  T_values_to_show <- c(100, 400, 1600)
  
  # Check if bias_results have required structure
  if (!all(c("bias_results" %in% names(bias_results_DGP1), "bias_results" %in% names(bias_results_DGP2)))) {
    warning("Invalid bias_results structure")
    return(NULL)
  }
  
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
  
  # Prepare data for plotting
  combined_data <- data.frame()
  
  dgp_names <- c("DGP-1", "DGP-2")
  dgp_results <- list(bias_results_DGP1, bias_results_DGP2)
  
  for (dgp_idx in 1:2) {
    dgp_name <- dgp_names[dgp_idx]
    bias_results <- dgp_results[[dgp_idx]]
    
    for (T_value in T_values_to_show) {
      T_name <- paste0("T", T_value)
      
      if (T_name %in% names(bias_results$bias_results)) {
        # Extract bias values
        bias_res <- bias_results$bias_results[[T_name]]
        
        # Check if all required keys exist
        required_keys <- c("var_bias_R11", "lp_bias_R11", "var_bias_R21", "lp_bias_R21",
                           "var_bias_R12", "lp_bias_R12", "var_bias_R22", "lp_bias_R22")
        
        if (all(required_keys %in% names(bias_res))) {
          var_bias_R11 <- bias_res$var_bias_R11
          lp_bias_R11 <- bias_res$lp_bias_R11
          var_bias_R21 <- bias_res$var_bias_R21
          lp_bias_R21 <- bias_res$lp_bias_R21
          var_bias_R12 <- bias_res$var_bias_R12
          lp_bias_R12 <- bias_res$lp_bias_R12
          var_bias_R22 <- bias_res$var_bias_R22
          lp_bias_R22 <- bias_res$lp_bias_R22
          
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
  }
  
  # Check if we have data to plot
  if (nrow(combined_data) == 0) {
    warning("No valid bias data for combined figure")
    return(NULL)
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
      } else {
        # Create a placeholder plot for missing data
        p <- ggplot() +
          annotate("text", x = 7.5, y = 0, label = "No bias data available") +
          labs(title = paste0(gsub("-", " ", dgp_name), ". ", response, "."),
               x = "τ", y = "Bias") +
          theme_minimal() +
          theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "grey92", size = 0.2),
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.title = element_text(size = 9),
            panel.border = element_rect(fill = NA, color = "grey80", size = 0.3)
          ) +
          xlim(0, 15) + ylim(y_limits[[plot_key]])
        
        plots[[plot_key]] <- p
      }
    }
  }
  
  # Check if we have all needed plots
  if (length(plots) < 8) {
    warning("Not all required plots could be created for combined bias figure")
    return(NULL)
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
    top = grid::textGrob("Figure 1: Simulation-Estimated Bias of LP-Based and VAR-Based Estimators", 
                         gp = grid::gpar(fontsize = 11, fontface = "bold"))
  )
  
  return(combined_figure)
}

#--------------------------------------------------------------
# Create Tables Similar to Tables 1 and 2 in the Paper
#--------------------------------------------------------------

# Function to create rejection rate tables
create_rejection_table <- function(dm_results, dgp_name) {
  # Check if dm_results has required structure
  if (!("dm_results" %in% names(dm_results))) {
    warning("Invalid dm_results structure for ", dgp_name)
    return(NULL)
  }
  
  # Extract sample sizes

  sample_sizes <- as.numeric(gsub("T", "", names(dm_results$dm_results)))
  # decreasing False
  sample_sizes <- sort(sample_sizes)
  
  
  if (length(sample_sizes) == 0) {
    warning("No sample sizes found in dm_results for ", dgp_name)
    return(NULL)
  }
  
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
      dm_res <- dm_results$dm_results[[T_name]]
      
      if (all(c("rejection_R11", "rejection_R21", "rejection_R12", "rejection_R22") %in% names(dm_res))) {
        table_R11[i, ] <- dm_res$rejection_R11
        table_R21[i, ] <- dm_res$rejection_R21
        table_R12[i, ] <- dm_res$rejection_R12
        table_R22[i, ] <- dm_res$rejection_R22
      }
    }
  }
  
  # Set row and column names
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