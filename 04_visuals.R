# Direct analysis of MCMC chains for CJS model results
# This approach bypasses the parameter extraction issues by directly analyzing the MCMC chains
rm(list = ls())
gc()

# Load necessary libraries
library(tidyverse)
library(coda)     # For MCMC analysis
library(ggplot2)  # For visualization
library(viridis)  # For color scales

# Set theme for consistent visualization
theme_set(theme_classic())

#========================================================================
# 1. Load Results
#========================================================================
# List all files in the directory that contain "min_captures"
min_capture_files <- list.files(path = "Data/",
                                pattern = ".*min_capture.*\\.rds$", 
                                full.names = TRUE)

# Create a data frame with file information
file_info <- data.frame(
  file_path = min_capture_files,
  file_name = basename(min_capture_files),
  stringsAsFactors = FALSE
)

# Extract bird code and sample size from filenames
file_info <- file_info %>%
  mutate(
    bird_code = str_extract(file_name, "^[A-Z]{4}"),
    sample_size = as.numeric(str_extract(file_name, "\\d+(?=\\.rds$)"))
  )

# Print the files found
cat("Found", nrow(file_info), "files matching the pattern:\n")
print(file_info)

# Now load each file into a list
results_list <- list()

for (i in 1:nrow(file_info)) {
  bird_code <- file_info$bird_code[i]
  file_path <- file_info$file_path[i]
  
  # Skip files with missing bird code
  if (is.na(bird_code)) next
  
  # Try to load the file
  tryCatch({
    results_list[[bird_code]] <- readRDS(file_path)
    cat("Loaded file for", bird_code, "\n")
  }, error = function(e) {
    cat("Error loading file for", bird_code, ":", e$message, "\n")
  })
}

#========================================================================
# 2. Direct analysis of MCMC chains
#========================================================================
# Function to extract summary statistics directly from MCMC chains
extract_mcmc_summary <- function(results_list) {
  all_params <- data.frame()
  
  for (sp in names(results_list)) {
    cat("\nExtracting parameters for", sp, "\n")
    
    # Get the MCMC data
    mcmc_data <- results_list[[sp]]
    
    # First try to handle the format you showed (with $chain1)
    if (!is.null(mcmc_data$chain1) && is.matrix(mcmc_data$chain1)) {
      cat("  Found chain1 matrix format\n")
      
      # Extract chain data
      chain_data <- mcmc_data$chain1
      
      # Get parameter names
      param_names <- colnames(chain_data)
      
      # Calculate summary statistics
      param_means <- colMeans(chain_data)
      param_sds <- apply(chain_data, 2, sd)
      param_q025 <- apply(chain_data, 2, quantile, probs = 0.025)
      param_q975 <- apply(chain_data, 2, quantile, probs = 0.975)
      
      # Create a data frame with results
      sp_params <- data.frame(
        parameter = param_names,
        mean = param_means,
        sd = param_sds,
        lower_ci = param_q025,
        upper_ci = param_q975,
        species = sp,
        sample_size = file_info$sample_size[file_info$bird_code == sp]
      )
      
      # Add to combined results
      all_params <- rbind(all_params, sp_params)
      
    } else if (inherits(mcmc_data, "mcmc") || inherits(mcmc_data, "mcmc.list")) {
      # Handle standard mcmc objects
      cat("  Found standard mcmc or mcmc.list format\n")
      
      # Get summary statistics using coda
      mcmc_summary <- summary(mcmc_data)
      
      # Extract parameter statistics
      param_stats <- as.data.frame(mcmc_summary$statistics)
      param_quant <- as.data.frame(mcmc_summary$quantiles)
      
      # Create parameter data frame
      sp_params <- data.frame(
        parameter = rownames(param_stats),
        mean = param_stats$Mean,
        sd = param_stats$SD,
        lower_ci = param_quant[,"2.5%"],
        upper_ci = param_quant[,"97.5%"],
        species = sp,
        sample_size = file_info$sample_size[file_info$bird_code == sp]
      )
      
      # Add to combined results
      all_params <- rbind(all_params, sp_params)
      
    } else {
      cat("  Unrecognized format - skipping\n")
      next
    }
  }
  
  cat("\nTotal parameters extracted:", nrow(all_params), "\n")
  return(all_params)
}

# Extract parameter summaries
mcmc_summaries <- extract_mcmc_summary(results_list)

#========================================================================
# 3. Categorize parameters and prepare for visualization
#========================================================================
# Function to identify parameter types from names
categorize_parameters <- function(param_df) {
  # Make a copy to avoid modifying the original
  result_df <- param_df
  
  # Add columns for parameter type, station, and time
  result_df$type <- "Unknown"
  result_df$station <- 1
  result_df$time <- 1
  
  # Process each parameter name
  for (i in 1:nrow(result_df)) {
    param_name <- result_df$parameter[i]
    
    # Skip if parameter name is missing
    if (is.na(param_name) || param_name == "") next
    
    # Identify parameter type
    if (grepl("^phi", param_name)) {
      result_df$type[i] <- "Survival"
    } else if (grepl("^p\\[", param_name) || grepl("^p\\.", param_name)) {
      result_df$type[i] <- "Detection"
    }
    
    # Try to extract station and time if the parameter has brackets
    if (grepl("\\[", param_name)) {
      # Extract indices from brackets: format like p[time, species, station]
      indices <- gsub(".*\\[(.*)\\].*", "\\1", param_name)
      index_parts <- strsplit(indices, ",")[[1]]
      
      # Clean up whitespace and convert to numbers
      index_parts <- trimws(index_parts)
      
      # Extract time (first index)
      if (length(index_parts) >= 1) {
        result_df$time[i] <- as.numeric(index_parts[1])
      }
      
      # Extract station (last index, or second-to-last if more than 2 indices)
      if (length(index_parts) >= 3) {
        # Format is likely [time, species, station]
        result_df$station[i] <- as.numeric(index_parts[3])
      } else if (length(index_parts) == 2) {
        # Format might be [time, station]
        result_df$station[i] <- as.numeric(index_parts[2])
      }
    }
  }
  
  return(result_df)
}

# Categorize parameters
categorized_params <- categorize_parameters(mcmc_summaries)

# Print summary of categorized parameters
cat("\nParameter Categorization Results:\n")
cat("Parameter types found:", paste(unique(categorized_params$type), collapse = ", "), "\n")
cat("Number of stations:", length(unique(categorized_params$station)), "\n")
cat("Number of time periods:", length(unique(categorized_params$time)), "\n")

# Check number of parameters by type
type_counts <- table(categorized_params$type)
cat("\nParameter counts by type:\n")
print(type_counts)

#========================================================================
# 4. Visualize survival and detection probabilities
#========================================================================
# Basic parameter visualization function
plot_parameters <- function(param_data) {
  # Filter for known parameter types
  known_params <- param_data %>% 
    filter(type %in% c("Survival", "Detection"))
  
  if (nrow(known_params) == 0) {
    cat("No valid parameters found for plotting\n")
    return(NULL)
  }
  
  # Create plot
  param_plot <- ggplot(known_params, aes(x = interaction(species, station), 
                                         y = mean, 
                                         color = species)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
    facet_wrap(~type) +
    labs(title = "Parameter Estimates by Species and Station",
         x = "Species-Station", 
         y = "Parameter Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_viridis_d()
  
  return(param_plot)
}

# Generate and print the plot
parameter_plot <- plot_parameters(categorized_params)
if (!is.null(parameter_plot)) {
  print(parameter_plot)
}

#========================================================================
# 5. Visualize time series of parameter values
#========================================================================
# Plot parameters over time
plot_time_series <- function(param_data) {
  # Check if we have time information
  if (length(unique(param_data$time)) <= 1) {
    cat("No time series information available\n")
    return(NULL)
  }
  
  # Filter for known parameter types
  known_params <- param_data %>% 
    filter(type %in% c("Survival", "Detection"))
  
  if (nrow(known_params) == 0) {
    cat("No valid parameters found for time series plotting\n")
    return(NULL)
  }
  
  # Create time series plot by species and station
  time_plot <- ggplot(known_params, aes(x = time, y = mean, color = type)) +
    geom_line() +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = type), alpha = 0.2, color = NA) +
    facet_grid(species ~ station, labeller = labeller(
      station = function(x) paste("Station", x),
      species = function(x) paste("Species", x)
    )) +
    labs(title = "Parameter Estimates Over Time",
         x = "Time Period", 
         y = "Parameter Value") +
    scale_color_manual(values = c("Survival" = "blue", "Detection" = "red")) +
    scale_fill_manual(values = c("Survival" = "lightblue", "Detection" = "pink"))
  
  return(time_plot)
}

# Generate and print the time series plot
time_series_plot <- plot_time_series(categorized_params)
if (!is.null(time_series_plot)) {
  print(time_series_plot)
}

#========================================================================
# 6. Compare precision of estimates across species
#========================================================================
# Precision analysis
analyze_precision <- function(param_data) {
  # Calculate precision (width of credible interval)
  param_data$precision <- param_data$upper_ci - param_data$lower_ci
  
  # Summarize precision by species and type
  precision_summary <- param_data %>%
    filter(type %in% c("Survival", "Detection")) %>%
    group_by(species, type) %>%
    summarize(
      mean_precision = mean(precision, na.rm = TRUE),
      sd_precision = sd(precision, na.rm = TRUE),
      median_precision = median(precision, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create precision plot
  precision_plot <- ggplot(precision_summary, 
                           aes(x = species, y = mean_precision, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_precision - sd_precision, 
                      ymax = mean_precision + sd_precision),
                  position = position_dodge(width = 0.9), width = 0.2) +
    labs(title = "Precision of Parameter Estimates by Species",
         subtitle = "Lower values indicate higher precision (narrower credible intervals)",
         x = "Species",
         y = "Mean Width of 95% Credible Interval",
         fill = "Parameter Type") +
    scale_fill_manual(values = c("Survival" = "blue", "Detection" = "red")) +
    theme(legend.position = "bottom")
  
  return(list(summary = precision_summary, plot = precision_plot))
}

# Generate precision analysis
precision_results <- analyze_precision(categorized_params)
print(precision_results$plot)
print(precision_results$summary)

#========================================================================
# 7. Compare sample size vs. precision
#========================================================================
# This analysis can only be done if we have multiple sample sizes per species
plot_sample_size_vs_precision <- function(param_data) {
  # Check if we have sample size variation
  sample_sizes <- unique(param_data$sample_size)
  
  if (length(sample_sizes) <= 1) {
    cat("No variation in sample size - can't analyze relationship\n")
    return(NULL)
  }
  
  # Calculate precision by sample size, species, and parameter type
  sample_precision <- param_data %>%
    filter(type %in% c("Survival", "Detection")) %>%
    mutate(precision = upper_ci - lower_ci) %>%
    group_by(species, sample_size, type) %>%
    summarize(
      mean_precision = mean(precision, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create plot
  sample_plot <- ggplot(sample_precision, 
                        aes(x = sample_size, y = mean_precision, color = species)) +
    geom_line() +
    geom_point(size = 3) +
    facet_wrap(~type) +
    labs(title = "Effect of Sample Size on Parameter Precision",
         x = "Sample Size", 
         y = "Mean Width of 95% Credible Interval",
         color = "Species") +
    scale_color_viridis_d() +
    theme(legend.position = "bottom")
  
  return(sample_plot)
}

# Generate and print the sample size vs. precision plot
sample_size_plot <- plot_sample_size_vs_precision(categorized_params)
if (!is.null(sample_size_plot)) {
  print(sample_size_plot)
}

#========================================================================
# 8. Direct access to MCMC chains for advanced diagnostics
#========================================================================
# Function to visualize specific parameters from the raw MCMC chains
plot_mcmc_chains <- function(results_list, species_code, param_pattern = "phi") {
  # Get the MCMC data for this species
  mcmc_data <- results_list[[species_code]]
  
  if (is.null(mcmc_data)) {
    cat("No data found for species", species_code, "\n")
    return(NULL)
  }
  
  # Handle different formats
  if (!is.null(mcmc_data$chain1) && is.matrix(mcmc_data$chain1)) {
    # Format with separate chains
    chain1 <- mcmc_data$chain1
    
    # Find parameters matching the pattern
    matching_params <- grep(param_pattern, colnames(chain1), value = TRUE)
    
    if (length(matching_params) == 0) {
      cat("No parameters matching", param_pattern, "found for", species_code, "\n")
      return(NULL)
    }
    
    # Limit to first few parameters for readability
    display_params <- matching_params[1:min(6, length(matching_params))]
    
    # Create trace plots
    par(mfrow = c(length(display_params), 1))
    for (param in display_params) {
      plot(chain1[, param], type = "l", 
           main = paste(species_code, "-", param),
           xlab = "Iteration", ylab = "Value")
    }
    par(mfrow = c(1, 1))
    
    # Create density plots
    par(mfrow = c(2, 3))
    for (param in display_params) {
      density_obj <- density(chain1[, param])
      plot(density_obj, main = paste(species_code, "-", param), 
           xlab = "Value", ylab = "Density")
    }
    par(mfrow = c(1, 1))
    
    cat("Created trace and density plots for", species_code, "\n")
    
  } else if (inherits(mcmc_data, "mcmc") || inherits(mcmc_data, "mcmc.list")) {
    # Standard mcmc format
    # Use coda's plot methods
    matching_params <- grep(param_pattern, varnames(mcmc_data), value = TRUE)
    
    if (length(matching_params) == 0) {
      cat("No parameters matching", param_pattern, "found for", species_code, "\n")
      return(NULL)
    }
    
    # Limit to first few parameters for readability
    display_params <- matching_params[1:min(6, length(matching_params))]
    
    # Create trace and density plots
    plot(mcmc_data[, display_params])
    
    cat("Created trace and density plots for", species_code, "\n")
  } else {
    cat("Unrecognized MCMC format for", species_code, "\n")
    return(NULL)
  }
  
  return(TRUE)  # Indicate success
}

# Direct visualization of MCMC chains for each species
for (sp in names(results_list)) {
  cat("\nCreating diagnostic plots for", sp, "\n")
  
  # Create trace plots for survival parameters
  cat("Survival parameter diagnostics:\n")
  plot_mcmc_chains(results_list, sp, "phi")
  
  # Create trace plots for detection parameters
  cat("Detection parameter diagnostics:\n")
  plot_mcmc_chains(results_list, sp, "^p\\[")
}

#========================================================================
# 9. Create summary tables for reporting
#========================================================================
# Function to generate summary tables
create_summary_tables <- function(param_data) {
  # Summary by species and parameter type
  species_type_summary <- param_data %>%
    filter(type %in% c("Survival", "Detection")) %>%
    group_by(species, type) %>%
    summarize(
      Mean = mean(mean, na.rm = TRUE),
      SD = sd(mean, na.rm = TRUE),
      Min = min(mean, na.rm = TRUE),
      Max = max(mean, na.rm = TRUE),
      `CI Width` = mean(upper_ci - lower_ci, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(type, species)
  
  # Summary by species, station, and parameter type
  station_summary <- param_data %>%
    filter(type %in% c("Survival", "Detection")) %>%
    group_by(species, station, type) %>%
    summarize(
      Mean = mean(mean, na.rm = TRUE),
      SD = sd(mean, na.rm = TRUE),
      `CI Width` = mean(upper_ci - lower_ci, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(type, species, station)
  
  # Return both tables
  return(list(
    species_summary = species_type_summary,
    station_summary = station_summary
  ))
}

# Generate summary tables
summary_tables <- create_summary_tables(categorized_params)

# Print summary tables
cat("\nSummary by Species and Parameter Type:\n")
print(summary_tables$species_summary)

cat("\nSummary by Species, Station, and Parameter Type:\n")
print(summary_tables$station_summary)

#========================================================================
# 10. Minimum sample size analysis
#========================================================================
# This function determines the minimum sample size for reliable estimates
analyze_min_sample_size <- function(param_data, precision_threshold = 0.3) {
  # Calculate precision
  param_data$precision <- param_data$upper_ci - param_data$lower_ci
  
  # Calculate mean precision by species, sample size, and parameter type
  precision_by_sample <- param_data %>%
    filter(type %in% c("Survival", "Detection")) %>%
    group_by(species, sample_size, type) %>%
    summarize(
      mean_precision = mean(precision, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Find minimum sample size that meets precision threshold for each species
  min_samples <- precision_by_sample %>%
    filter(type == "Survival") %>%  # Focus on survival parameters
    filter(mean_precision <= precision_threshold) %>%
    group_by(species) %>%
    summarize(
      min_sample_size = min(sample_size, na.rm = TRUE),
      achieved_precision = mean_precision[which.min(sample_size)],
      .groups = "drop"
    ) %>%
    arrange(min_sample_size)
  
  # Create bar plot of minimum sample sizes
  min_sample_plot <- ggplot(min_samples, aes(x = reorder(species, min_sample_size), 
                                          y = min_sample_size, 
                                          fill = species)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = min_sample_size), vjust = -0.5, size = 4) +
    labs(title = "Minimum Birds Per Station Required for Reliable Estimates",
         subtitle = paste0("For precision threshold of ", precision_threshold),
         x = "Species", 
         y = "Minimum Sample Size") +
    theme(legend.position = "none") +
    scale_fill_viridis_d()
  
  return(list(min_samples = min_samples, plot = min_sample_plot))
}

# Only run this if we have sample size variation
if (length(unique(categorized_params$sample_size)) > 1) {
  min_sample_results <- analyze_min_sample_size(categorized_params)
  print(min_sample_results$plot)
  cat("\nMinimum Sample Size Requirements by Species:\n")
  print(min_sample_results$min_samples)
} else {
  cat("\nNo sample size variation - can't determine minimum sample size\n")
}

#========================================================================
# 11. Export results (optional)
#========================================================================
# Uncomment to save processed results to CSV
# write.csv(categorized_params, "cjs_processed_parameters.csv", row.names = FALSE)
# write.csv(summary_tables$species_summary, "cjs_species_summary.csv", row.names = FALSE)
# write.csv(summary_tables$station_summary, "cjs_station_summary.csv", row.names = FALSE)

# Uncomment to export plots to PDF
# pdf("cjs_visualization_results.pdf", width = 10, height = 8)
# # Re-run plot code here
# dev.off()

#========================================================================
# 12. Detailed example of working with one species
#========================================================================
# Function to do a focused analysis for a specific species
analyze_single_species <- function(results_list, species_code) {
  cat("\n======================================================\n")
  cat("Detailed Analysis for Species:", species_code, "\n")
  cat("======================================================\n")
  
  # Get the MCMC data for this species
  mcmc_data <- results_list[[species_code]]
  
  if (is.null(mcmc_data)) {
    cat("No data found for species", species_code, "\n")
    return(NULL)
  }
  
  # 1. Examine the data structure
  cat("Data structure:\n")
  if (!is.null(mcmc_data$chain1) && is.matrix(mcmc_data$chain1)) {
    cat("  Format: List with chain matrices\n")
    cat("  Chain1 dimensions:", nrow(mcmc_data$chain1), "x", ncol(mcmc_data$chain1), "\n")
    cat("  Parameter names (first 5):", paste(colnames(mcmc_data$chain1)[1:min(5, ncol(mcmc_data$chain1))], collapse = ", "), "...\n")
    
    # 2. Extract all phi parameters
    phi_params <- grep("^phi", colnames(mcmc_data$chain1), value = TRUE)
    cat("\nSurvival parameters found:", length(phi_params), "\n")
    if (length(phi_params) > 0) {
      cat("  Example names:", paste(phi_params[1:min(3, length(phi_params))], collapse = ", "), "...\n")
      
      # Calculate summary statistics for phi parameters
      phi_means <- colMeans(mcmc_data$chain1[, phi_params, drop = FALSE])
      phi_sds <- apply(mcmc_data$chain1[, phi_params, drop = FALSE], 2, sd)
      phi_q025 <- apply(mcmc_data$chain1[, phi_params, drop = FALSE], 2, quantile, probs = 0.025)
      phi_q975 <- apply(mcmc_data$chain1[, phi_params, drop = FALSE], 2, quantile, probs = 0.975)
      
      # Create summary data frame
      phi_summary <- data.frame(
        parameter = phi_params,
        mean = phi_means,
        sd = phi_sds,
        lower_ci = phi_q025,
        upper_ci = phi_q975
      )
      
      cat("\nSummary statistics for survival parameters:\n")
      print(head(phi_summary))
      
      # 3. Extract all p parameters
      p_params <- grep("^p\\[", colnames(mcmc_data$chain1), value = TRUE)
      cat("\nDetection parameters found:", length(p_params), "\n")
      if (length(p_params) > 0) {
        cat("  Example names:", paste(p_params[1:min(3, length(p_params))], collapse = ", "), "...\n")
        
        # Calculate summary statistics for p parameters
        p_means <- colMeans(mcmc_data$chain1[, p_params, drop = FALSE])
        p_sds <- apply(mcmc_data$chain1[, p_params, drop = FALSE], 2, sd)
        p_q025 <- apply(mcmc_data$chain1[, p_params, drop = FALSE], 2, quantile, probs = 0.025)
        p_q975 <- apply(mcmc_data$chain1[, p_params, drop = FALSE], 2, quantile, probs = 0.975)
        
        # Create summary data frame
        p_summary <- data.frame(
          parameter = p_params,
          mean = p_means,
          sd = p_sds,
          lower_ci = p_q025,
          upper_ci = p_q975
        )
        
        cat("\nSummary statistics for detection parameters:\n")
        print(head(p_summary))
        
        # 4. Create visualizations
        # Extract station information (if present in parameter names)
        extract_station <- function(param_name) {
          station_match <- regexpr(",(\\s*\\d+)(\\s*\\])$", param_name)
          if (station_match > 0) {
            station_text <- regmatches(param_name, station_match)
            return(as.numeric(gsub("[^0-9]", "", station_text)))
          } else {
            return(1)  # Default to station 1
          }
        }
        
        # Add station information
        phi_summary$station <- sapply(phi_summary$parameter, extract_station)
        p_summary$station <- sapply(p_summary$parameter, extract_station)
        
        # Add type information
        phi_summary$type <- "Survival"
        p_summary$type <- "Detection"
        
        # Combine for plotting
        combined_summary <- rbind(phi_summary, p_summary)
        
        # Create parameter plot by station
        param_plot <- ggplot(combined_summary, aes(x = as.factor(station), y = mean, color = type)) +
          geom_point(position = position_dodge(width = 0.3), size = 3) +
          geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                        position = position_dodge(width = 0.3), width = 0.2) +
          labs(title = paste(species_code, "Parameter Estimates by Station"),
               x = "Station", 
               y = "Parameter Estimate",
               color = "Parameter Type") +
          ylim(0, 1) +
          scale_color_manual(values = c("Survival" = "blue", "Detection" = "red"))
        
        print(param_plot)
        
        return(list(
          phi_summary = phi_summary,
          p_summary = p_summary,
          param_plot = param_plot
        ))
      }
    }
  } else if (inherits(mcmc_data, "mcmc") || inherits(mcmc_data, "mcmc.list")) {
    # Handle standard mcmc format
    cat("  Format: mcmc or mcmc.list object\n")
    
    # Use coda's summary method
    mcmc_summary <- summary(mcmc_data)
    
    cat("  Parameters:", length(varnames(mcmc_data)), "\n")
    cat("  Example names:", paste(varnames(mcmc_data)[1:min(5, length(varnames(mcmc_data)))], collapse = ", "), "...\n")
    
    # Further processing would be similar to above
    # ...
  } else {
    cat("  Unrecognized format\n")
  }
}

# Example: analyze a specific species (replace "BCCH" with your species code)
if ("BTNW" %in% names(results_list)) {
  analyze_single_species(results_list, "BTNW")
}

# Print session info for reproducibility
sessionInfo()