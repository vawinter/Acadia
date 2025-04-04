##########################################################################X
# CJS model to estimate survival from MAPS banding stations
#------------------------------------------------------------X
# What is the minimum number of birds per species necessary?
# VAW
# 4/3/2025
##########################################################################X
# Clean env
rm(list=ls())
gc()
# reproducibility
set.seed(0235)

# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(nimble)
library(MCMCvis)
library(coda)

# Function to create initial values for latent state z
init_z <- function(y_array) {
  # Dimensions: [N, time, sp, station]
  dims <- dim(y_array)
  N <- dims[1]
  time <- dims[2]
  sp <- dims[3]
  station <- dims[4]
  
  # Create an array for z
  z <- array(NA, dim = c(N, time, sp, station))
  
  for (s in 1:sp) {
    for (l in 1:station) {
      for (i in 1:N) {
        # Extract capture history for this individual
        y_vec <- y_array[i, , s, l]
        
        # First occasion - set initial state based on first capture (if any)
        first_capture <- which(y_vec == 1)[1]
        if (!is.na(first_capture)) {
          z[i, first_capture, s, l] <- 1  # Mark as alive at first capture
        }
        
        for (t in 2:time) {
          if (sum(y_vec[t:time]) > 0) {
            # If detected at or after t, must be alive at t
            z[i, t, s, l] <- 1
          } else if (sum(y_vec[1:(t-1)]) == 0) {
            # If never detected before t, assume not alive (conservative)
            z[i, t, s, l] <- 0
          } else {
            # Detected before t but not after, uncertain status - initialize as alive
            z[i, t, s, l] <- NA  # Uncertain if it was alive but not detected after
          }
        }
      }
    }
  }
  
  return(z)
}

# Load banding data
banding_data <- read_excel('MAPS atlantic forest BCR banding data.xlsx', sheet = 'MAPS_BANDING_capture_data')
station_data <- read.csv('MAPS atlantic forest BCR station info total.csv', header = T)

############################################X
# Define NIMBLE Model ----
############################################X
cjs_code <- nimbleCode({
  for (s in 1:sp) {  # Loop over species
    for (l in 1:station) {  # Loop over stations
      # Global parameters for this species-station combo
      phi.mean[s, l] ~ dunif(0, 1)  # Mean survival probability
      p.mean[s, l] ~ dunif(0, 1)    # Mean capture probability
      
      for (t in 1:(time-1)) {
        phi[t, s, l] <- phi.mean[s, l]  # Use same survival for all time periods
      }
      
      for (t in 1:time) {
        p[t, s, l] <- p.mean[s, l]  # Use same capture prob for all time periods
      }
      
      for (i in 1:N) {  # Loop over individuals
        # First occasion - all individuals alive at first occasion (simplified)
        z[i, 1, s, l] <- 1
        y[i, 1, s, l] ~ dbern(p[1, s, l] * z[i, 1, s, l])
        
        for (t in 2:time) {
          # State process: survival from t-1 to t
          z[i, t, s, l] ~ dbern(phi[t-1, s, l] * z[i, t-1, s, l])
          
          # Observation process: detection at time t
          y[i, t, s, l] ~ dbern(p[t, s, l] * z[i, t, s, l])
        }
      }
    }
  }
})

############################################X
# Minimum captures? ----
############################################X
# First, get a random sample of 6 stations
sampled_stations <- banding_data %>%
  distinct(STATION) %>%
  slice_sample(n = 6) %>%
  pull(STATION)

species_list <- c("BCCH", "BTNW", "HETH", "REVI")

# Then filter your data to only include those stations
full_data <- banding_data %>%
  filter(SPEC %in% species_list & STATION %in% sampled_stations) %>%
  select(STATION, year, SPEC, BAND)

# Function to test minimum sample size needed
test_min_recapture <- function(full_data, species, step = 2) {
  # Load dependencies
  library(nimble)
  library(coda)
  #----------------------X
  # Data prep ----
  #----------------------X
  # Filter species of interest
  species_data <- full_data %>% filter(SPEC == species)
  
  # Max captures per species
  captures <- n_distinct(species_data$BAND)

  for (i in 0:floor(captures)) {
    n <- captures - i * step
    cat("\nTesting with", n, "recaptures for", species, "\n")
    #----------------------X
    # Resample data ----
    sampled_data <- species_data %>% sample_n(n, replace = TRUE)
    #----------------------X
    # Prepare capture histories
    filtered_data <- sampled_data %>%
      dplyr::select(STATION, year, SPEC, BAND)
    
    # Assign numerical indices for categorical variables
    filtered_data <- filtered_data %>%
      mutate(
        station_id = as.integer(as.factor(STATION)),
        species_id = as.integer(as.factor(SPEC)),
        year_id = as.integer(as.factor(year))
      )
    
    # Get dataset dimensions
    N <- length(unique(filtered_data$BAND))  # Total unique birds
    sp <- 1  # Number of species
    time <- length(unique(filtered_data$year))  # Number of years
    station <- length(unique(filtered_data$STATION))  # Number of stations
    
    # Create empty capture history array
    capture_history <- array(0, dim = c(N, time, sp, station))
    
    # Fill in capture history: 1 if the bird was captured, 0 otherwise
    for (i in 1:nrow(filtered_data)) {
      bird <- 1
      yr <- filtered_data$year_id[i]
      stn <- filtered_data$station_id[i]
      spc <- filtered_data$species_id[i]
      
      capture_history[bird, yr, spc, stn] <- 1
    }
    
    
    #----------------------X
    # Prepare model specifications
    #----------------------X
    constants <- list(N = nrow(capture_history), 
                      time = ncol(capture_history), 
                      sp = 1, 
                      station = length(unique(sampled_data$STATION)))
    
    data_list <- list(y = capture_history)
    
    # Initializing z based on capture history [N x time x station x sp]
    z_inits <- init_z(capture_history)
    
    inits <- list(
      phi.mean = array(0.5, dim = c(sp, station)),
      p.mean = array(0.5, dim = c(sp, station)),
      z = z_inits
    )
    #----------------------X
    # Build model
    #----------------------X
    cjs_model <- nimbleModel(cjs_code, data = data_list, constants = constants, inits = inits)
    
    # Compile model
    compiled_model <- compileNimble(cjs_model)
    
    # Set up MCMC
    mcmc_conf <- configureMCMC(
      cjs_model, 
      monitors = c("phi", "p")
    )
    
    # Build and compile MCMC
    mcmc <- buildMCMC(mcmc_conf)
    compiled_mcmc <- compileNimble(mcmc)
    
    #----------------------X
    # Run MCMC sampling with error handling
    #----------------------X
    tryCatch({
      samples <- runMCMC(
        compiled_mcmc, 
        niter = 130000, 
        nburnin = 10000, 
        nchains = 2,
        thin = 10,
        samplesAsCodaMCMC = TRUE
      )
      
      # If successful, check convergence and ESS
      summary_df <- MCMCsummary(samples)
      phi_values <- summary_df[grep("phi", rownames(summary_df)), ]
      rhat <- phi_values$Rhat
      ess <- phi_values$n.eff
      
      # Convergence check
      cat("\nRhat values for", n, "recaptures:\n")
      print(rhat)
      
      # Convergence check
      cat("\nESS values for", n, "recaptures:\n")
      print(ess)
      
      if (any(rhat > 1.1)) {
        cat("Convergence issue for", n, "recaptures, skipping to next\n")
        next
      }
      
      if (min(ess) > 100) {
        best_n <- n
        min_rhat <- min(rhat)
        
        # Summarize how many individuals came from each station
        station_counts <- sampled_data %>%
          count(STATION, name = "n_individuals")
        
        # Store results
        results_list[[as.character(n)]] <- list(
          sample_size = n,
          station_counts = station_counts,
          phi_summary = summary_df[grep("phi", rownames(summary_df)), ],
          p_summary = summary_df[grep("p", rownames(summary_df)), ]
        )
      }
      
    }, error = function(e) {
      cat("Model broke with error: ", e$message, "\n")
      
      # Output phi and p values before the model broke
      if (exists("phi_values")) {
        cat("phi values before break:\n")
        print(phi_values)
      }
      if (exists("samples")) {
        p_values <- MCMCsummary(samples)[grep("p", rownames(MCMCsummary(samples))), ]
        cat("p values before break:\n")
        print(p_values)
      }
      
      # Optionally, store partial results before the break
      results_list[[n]] <- list(phi = phi_values, p = p_values)
    })
    
    if (exists("best_n") && n > best_n + step) {
      cat("\nOptimal number of recaptures found for", species, ":", best_n, "with Rhat:", min_rhat)
      break
    }
  }
  
  return(results_list)
}
#---------------------------------X
## Run this function for each species ----
#---------------------------------X
# results_list <- list()
# species <- c("BCCH", "BTNW", "HETH", "REVI")
# for (sp in species) {
#   results <- test_min_recapture(filtered_data, sp)
#   cat("\nResults for", sp, ":\n")
#   print(results)
# }

all_results <- list()
for (sp in species) {
  all_results[[sp]] <- test_min_recapture(full_data, sp)
}

# #---------------------------------------------------
# # Set up parallel processing
# #---------------------------------------------------
# library(parallel)
# library(doParallel)
# 
# # Calculate number of cores to use (leave 1 free for system)
# num_cores <- detectCores() - 1
# cat("Running on", num_cores, "cores\n")
# 
# # Create cluster
# cl <- makeCluster(num_cores)
# 
# # Register parallel backend
# registerDoParallel(cl)
# 
# # Export needed libraries and functions to all workers
# clusterEvalQ(cl, {
#   library(nimble)
#   library(coda)
#   library(dplyr)
# })
# 
# # Export data to all workers
# clusterExport(cl, c("filtered_data"))
# 
# #---------------------------------------------------
# # Run parallel analysis
# #---------------------------------------------------
# species_list <- c("BCCH", "BTNW", "HETH", "REVI")
# 
# # Run in parallel
# results <- parLapply(cl, species_list, function(sp) {
#   return(test_min_recapture(filtered_data, sp, step = 2))  # Using step=5 to test fewer sample sizes
# })
# 
# # Stop cluster
# stopCluster(cl)
# 
# #---------------------------------------------------
# # Process and save results
# #---------------------------------------------------
# names(results) <- species_list
# 
# # Save results to file
# saveRDS(results, "cjs_min_sample_results.rds")
# 
# # Print summary
# for (sp in species_list) {
#   sp_results <- results[[sp]]
#   cat("\n=== Results for", sp, "===\n")
#   
#   # Find minimum successful sample size
#   min_sample <- NULL
#   for (n in names(sp_results$results)) {
#     result <- sp_results$results[[n]]
#     if (!is.null(result$converged) && result$converged) {
#       min_sample <- result$n_recaptures
#       break
#     }
#   }
#   
#   if (!is.null(min_sample)) {
#     cat("Minimum required recaptures:", min_sample, "\n")
#   } else {
#     cat("No convergent sample size found\n")
#   }
# }
# 
