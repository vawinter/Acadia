#Functions
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

# Function to test minimum sample size needed
test_min_recapture <- function(full_data, species, step = 2, iter, burn, t) {
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
        niter = iter, 
        nburnin = burn, 
        nchains = 2,
        thin = t,
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
        
        # # Summarize how many individuals came from each station
        # station_counts <- sampled_data %>%
        #   count(STATION, name = "n_individuals")
        # 
        # # Store results
        # results_list[[as.character(n)]] <- list(
        #   sample_size = n,
        #   station_counts = station_counts,
        #   phi_summary = summary_df[grep("phi", rownames(summary_df)), ],
        #   p_summary = summary_df[grep("p", rownames(summary_df)), ]
        # )
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
  
  return(samples)
}
