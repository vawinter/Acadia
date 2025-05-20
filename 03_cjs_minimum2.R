##########################################################################X
# CJS model to estimate survival from MAPS banding stations
#------------------------------------------------------------X
# What is the minimum number of birds per species necessary?
# Finding number by iterating through MAPS data and systematically
# removing samples each time.
# VAW
# 4/3/2025
##########################################################################X
# Clean env
rm(list=ls())
gc()
# reproducibility
set.seed(0235)

# Source dependencies
source("00_funs.R")

# Load banding data
banding_data <- read_excel('MAPS atlantic forest BCR banding data.xlsx', sheet = 'MAPS_BANDING_capture_data')
station_data <- read.csv('MAPS atlantic forest BCR station info total.csv', header = T)

############################################X
# Define NIMBLE Model ----
############################################X
# Init_z function adjusted for this specific CJS model
init_z <- function(y) {
  z <- array(0, dim = dim(y))
  for (i in 1:dim(y)[1]) {
    for (s in 1:dim(y)[3]) {
      for (l in 1:dim(y)[4]) {
        # Set first occasion to 1
        z[i, 1, s, l] <- 1
        
        # Initialize other occasions based on captures
        if (sum(y[i, , s, l]) > 0) {
          first <- which(y[i, , s, l] == 1)[1]
          last <- max(which(y[i, , s, l] == 1))
          
          if (first < dim(y)[2]) {
            for (j in first:last) {
              z[i, j, s, l] <- 1
            }
          }
        }
      }
    }
  }
  return(z)
}
test_min_recapture_per_station <- function(full_data, species, start_n = 10, step = 5, iter = 50000, burn = 10000, thin = 10) {
  library(nimble)
  library(coda)
  library(dplyr)
  library(MCMCvis)  # For MCMCsummary
  
  # Filter data for the specific species
  species_data <- full_data %>% filter(SPEC == species)
  stations <- unique(species_data$STATION)
  
  # Get max number of individuals per station
  station_counts <- species_data %>%
    group_by(STATION, year) %>%
    summarise(n = n_distinct(BAND), .groups = "drop") %>%
    group_by(STATION) %>%
    summarise(total = sum(n), .groups = "drop")
  
  max_n <- min(station_counts$total)  # Use the minimum across stations to ensure all can be sampled
  
  cat("Maximum available individuals at the most limited station:", max_n, "\n")
  
  # Create sequence of sample sizes to test
  sample_sizes <- seq(start_n, max_n, by = step)
  if(length(sample_sizes) == 0) {
    cat("No valid sample sizes to test. Adjust start_n or step size.\n")
    return(NULL)
  }
  
  results_list <- list()
  
  # Define CJS model code (assuming it's defined elsewhere)
  cjs_code <- nimbleCode({
    # Priors
    for (s in 1:sp) {
      for (t in 1:station) {
        phi.mean[s, t] ~ dbeta(1, 1)      # Survival probability
        p.mean[s, t] ~ dbeta(1, 1)        # Detection probability
      }
    }
    
    # Likelihood
    for (i in 1:N) {
      # First occasion
      for (s in 1:sp) {
        for (t in 1:station) {
          z[i, 1, s, t] <- y[i, 1, s, t]
          for (j in 2:time) {
            # State process
            z[i, j, s, t] ~ dbern(z[i, j-1, s, t] * phi.mean[s, t])
            # Observation process
            y[i, j, s, t] ~ dbern(z[i, j, s, t] * p.mean[s, t])
          }
        }
      }
    }
  })
  
  for (n_per_station in sample_sizes) {
    cat("\nTesting with", n_per_station, "individuals per station for", species, "\n")
    
    tryCatch({
      # Sample n_per_station individuals per station
      sampled_data <- data.frame()
      
      for (stn in stations) {
        station_data <- species_data %>% filter(STATION == stn)
        unique_bands <- unique(station_data$BAND)
        
        if (length(unique_bands) < n_per_station) {
          cat("Station", stn, "has fewer than", n_per_station, "individuals. Skipping this sample size.\n")
          next
        }
        
        # Randomly sample bands for this station
        sampled_bands <- sample(unique_bands, n_per_station)
        station_sampled <- station_data %>% filter(BAND %in% sampled_bands)
        sampled_data <- bind_rows(sampled_data, station_sampled)
      }
      
      if (nrow(sampled_data) == 0) {
        cat("No valid data after sampling. Skipping to next sample size.\n")
        next
      }
      
      # Process the filtered data
      filtered_data <- sampled_data %>%
        dplyr::select(STATION, year, SPEC, BAND) %>%
        mutate(
          station_id = as.integer(as.factor(STATION)),
          species_id = as.integer(as.factor(SPEC)),
          year_id = as.integer(as.factor(year))
        )
      
      # Setup dimensions
      N <- length(unique(filtered_data$BAND))
      sp <- 1
      time <- length(unique(filtered_data$year_id))
      station <- length(unique(filtered_data$STATION))
      
      # Empty capture history array
      capture_history <- array(0, dim = c(N, time, sp, station))
      
      # Fill capture history
      bird_ids <- unique(filtered_data$BAND)
      station_ids <- unique(filtered_data$station_id)
      
      for (i in 1:nrow(filtered_data)) {
        bird_idx <- which(bird_ids == filtered_data$BAND[i])
        yr_idx <- filtered_data$year_id[i]
        stn_idx <- which(station_ids == filtered_data$station_id[i])
        capture_history[bird_idx, yr_idx, 1, stn_idx] <- 1
      }
      
      # Check for sufficient recaptures
      recaptures <- apply(capture_history, c(1, 4), sum)
      recaptured_birds <- apply(recaptures > 1, 2, sum)
      
      cat("Birds with recaptures by station:\n")
      for (s in 1:length(station_ids)) {
        cat("Station", stations[s], ":", recaptured_birds[s], "birds with recaptures\n")
      }
      
      if (any(recaptured_birds < 5)) {
        cat("At least one station has fewer than 5 recaptured birds. This may lead to poor estimation.\n")
      }
      
      # Setup for NIMBLE
      constants <- list(
        N = dim(capture_history)[1],
        time = dim(capture_history)[2],
        sp = dim(capture_history)[3],
        station = dim(capture_history)[4]
      )
      
      data_list <- list(y = capture_history)
      z_inits <- init_z(capture_history)
      
      inits <- list(
        phi.mean = array(0.5, dim = c(sp, station)),
        p.mean = array(0.5, dim = c(sp, station)),
        z = z_inits
      )
      
      # Run NIMBLE model
      cjs_model <- suppressMessages(nimbleModel(cjs_code, data = data_list, constants = constants, inits = inits))
      compiled_model <- compileNimble(cjs_model)
      mcmc_conf <- configureMCMC(cjs_model, monitors = c("phi.mean", "p.mean"))
      mcmc <- buildMCMC(mcmc_conf)
      compiled_mcmc <- compileNimble(mcmc)
      
      samples <- runMCMC(
        compiled_mcmc,
        niter = iter,
        nburnin = burn,
        nchains = 2,
        thin = thin,
        samplesAsCodaMCMC = TRUE
      )
      
      # Analyze results
      summary_df <- MCMCsummary(samples)
      phi_values <- summary_df[grep("phi", rownames(summary_df)), ]
      p_values <- summary_df[grep("p", rownames(summary_df)), ]
      
      # Check convergence
      rhat_phi <- phi_values$Rhat
      ess_phi <- phi_values$n.eff
      
      cat("\nRhat for phi with", n_per_station, "per station:\n")
      print(rhat_phi)
      cat("\nESS for phi:\n")
      print(ess_phi)
      
      # Store results if convergence criteria met
      if (all(!is.na(rhat_phi)) && all(rhat_phi < 1.1) && all(ess_phi > 100)) {
        cat("✓ Convergence criteria met!\n")
        
        # Calculate precision (width of credible intervals)
        phi_precision <- phi_values$`97.5%` - phi_values$`2.5%`
        
        results_list[[as.character(n_per_station)]] <- list(
          n_per_station = n_per_station,
          station_counts = recaptured_birds,
          rhat = rhat_phi,
          ess = ess_phi,
          phi_summary = phi_values,
          p_summary = p_values,
          phi_precision = phi_precision,
          mean_precision = mean(phi_precision)
        )
        
        # If precision is good, we may have found our minimum
        if (mean(phi_precision) < 0.3) { # Adjust threshold as needed
          cat("Good precision achieved with", n_per_station, "individuals!\n")
        }
      } else {
        cat("✗ Convergence criteria not met.\n")
      }
      
    }, error = function(e) {
      cat("Error with", n_per_station, "individuals:", e$message, "\n")
    })
  }
  
  # Find minimum sample size that meets criteria
  if (length(results_list) > 0) {
    # Extract sample sizes that converged
    converged_sizes <- as.numeric(names(results_list))
    min_size <- min(converged_sizes)
    
    # Create summary dataframe
    summary_df <- data.frame(
      sample_size = converged_sizes,
      mean_precision = sapply(results_list, function(x) x$mean_precision),
      min_ess = sapply(results_list, function(x) min(x$ess)),
      max_rhat = sapply(results_list, function(x) max(x$rhat))
    )
    
    return(list(
      min_sample_size = min_size,
      all_results = results_list,
      summary = summary_df
    ))
  } else {
    cat("No sample sizes met convergence criteria.\n")
    return(NULL)
  }
}

############################################X
# Minimum captures? ----
############################################X
# First, get a random sample of stations
set.seed(125)  # For reproducibility
sampled_stations <- banding_data %>%
  filter(SPEC %in% c("BCCH", "BTNW", "HETH", "REVI")) %>%
  slice_sample(n = 10) %>%
  pull(STATION) %>%
  unique()

# Define study species and their respective MCMC settings
species_settings <- list(
  BCCH = list(iter = 350000, burn = 20000, thin = 10, start_n = 5),
  BTNW = list(iter = 130000, burn = 10000, thin = 10, start_n = 5),
  HETH = list(iter = 300000, burn = 20000, thin = 10, start_n = 5),
  REVI = list(iter = 200000, burn = 15000, thin = 10, start_n = 5)
)

# Species of interest
selected_species <- c("BCCH")  # Can be expanded to include other species

# Filter data to only include those stations and species
full_data <- banding_data %>%
  filter(SPEC %in% selected_species & STATION %in% sampled_stations) %>%
  select(STATION, year, SPEC, BAND)

# Run analysis for each species with appropriate settings
results_list <- list()

for (sp in selected_species) {
  cat("\n========================================\n")
  cat("Analyzing species:", sp, "\n")
  cat("========================================\n")
  
  settings <- species_settings[[sp]]
  
  results_list[[sp]] <- test_min_recapture_per_station(
    full_data = full_data,
    species = sp,
    start_n = settings$start_n,
    step = 2,  # Increase by 2 individuals each iteration
    iter = settings$iter,
    burn = settings$burn,
    thin = settings$thin
  )
  
  # Save results after each species (in case of crashes)
  saveRDS(results_list, file = paste0("min_capture_results_", Sys.Date(), ".rds"))
  
  cat("\nResults for", sp, ":\n")
  if (!is.null(results_list[[sp]])) {
    cat("Minimum sample size:", results_list[[sp]]$min_sample_size, "\n")
    print(results_list[[sp]]$summary)
  } else {
    cat("No convergent results found.\n")
  }
}

# Create summary across all species
species_summary <- data.frame(
  Species = character(),
  Min_Sample_Size = numeric(),
  Mean_Precision = numeric(),
  Max_Rhat = numeric(),
  Min_ESS = numeric(),
  stringsAsFactors = FALSE
)

for (sp in names(results_list)) {
  if (!is.null(results_list[[sp]])) {
    min_n <- results_list[[sp]]$min_sample_size
    min_result <- results_list[[sp]]$all_results[[as.character(min_n)]]
    
    species_summary <- rbind(species_summary, data.frame(
      Species = sp,
      Min_Sample_Size = min_n,
      Mean_Precision = min_result$mean_precision,
      Max_Rhat = max(min_result$rhat, na.rm = TRUE),
      Min_ESS = min(min_result$ess, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }
}

# Print final summary
cat("\n========================================\n")
cat("Final Summary Across Species\n")
cat("========================================\n")
print(species_summary)

# Save final results
saveRDS(results_list, file = "min_capture_results_final.rds")
write.csv(species_summary, file = "species_min_sample_sizes.csv", row.names = FALSE)