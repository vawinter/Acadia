##########################################################################
# CJS model to estimate survival from MAPS banding stations
#------------------------------------------------------------X
# Building CJS model and fitting with simulated data
# VAW
# 4/3/2025
##########################################################################
# Clean env
rm(list=ls())
gc()
# reproducibility
set.seed(0235)

# Source dependencies
source("00_funs.R")

# Define a simpler CJS model for debugging
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

# Set up model parameters - using smaller values for testing
N <- 20      # Number of individuals per species per station
sp <- 2      # Number of species
time <- 3    # Number of capture occasions
station <- 2 # Number of stations

# True parameter values for simulation
phi_true <- 0.7  # Survival probability
p_true <- 0.5    # Capture probability

# Simulate data - simplified version
simulate_cjs_data <- function(N, time, sp, station, phi, p) {
  # Initialize arrays
  z <- array(NA, dim=c(N, time, sp, station))
  y <- array(NA, dim=c(N, time, sp, station))
  
  for (s in 1:sp) {
    for (l in 1:station) {
      for (i in 1:N) {
        # First occasion
        z[i, 1, s, l] <- 1  # All animals alive at first occasion
        y[i, 1, s, l] <- rbinom(1, 1, p)  # First detection
        
        # Subsequent occasions
        for (t in 2:time) {
          # Survival process
          if (z[i, t-1, s, l] == 1) {
            z[i, t, s, l] <- rbinom(1, 1, phi)
          } else {
            z[i, t, s, l] <- 0  # Dead stays dead
          }
          
          # Observation process
          if (z[i, t, s, l] == 1) {
            y[i, t, s, l] <- rbinom(1, 1, p)
          } else {
            y[i, t, s, l] <- 0  # Can't detect if dead
          }
        }
      }
    }
  }
  
  return(list(z=z, y=y))
}

# Generate simulated data
sim_data <- simulate_cjs_data(N, time, sp, station, phi_true, p_true)

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
        
        # First occasion - everyone starts as alive
        z[i, 1, s, l] <- 1
        
        for (t in 2:time) {
          if (sum(y_vec[t:time]) > 0) {
            # If detected at or after t, must be alive at t
            z[i, t, s, l] <- 1
          } else if (sum(y_vec[1:(t-1)]) == 0) {
            # If never detected before t, assume not alive (conservative)
            z[i, t, s, l] <- 0
          } else {
            # Detected before t but not after, uncertain status - initialize as alive
            z[i, t, s, l] <- 1
          }
        }
      }
    }
  }
  
  return(z)
}

# Prepare model components
constants <- list(
  N = N, 
  sp = sp, 
  time = time, 
  station = station
)

data <- list(
  y = sim_data$y
)

# Check dimensions
cat("\nDimensions check:\n")
cat("Y dimensions:", paste(dim(data$y), collapse=" x "), "\n")

# Initialize z states
z_inits <- init_z(data$y)

# Debug: Check if z initialization is consistent
check_z_consistency <- function(y, z) {
  inconsistent <- 0
  n_dims <- length(dim(y))
  
  # Loop through all elements
  for (i in 1:dim(y)[1]) {
    for (t in 1:dim(y)[2]) {
      for (s in 1:dim(y)[3]) {
        for (l in 1:dim(y)[4]) {
          # If y[i,t,s,l] = 1 but z[i,t,s,l] = 0, we have an inconsistency
          if (y[i,t,s,l] == 1 && z[i,t,s,l] == 0) {
            cat("Inconsistency at i=", i, ", t=", t, ", s=", s, ", l=", l, 
                ": y=", y[i,t,s,l], ", z=", z[i,t,s,l], "\n")
            inconsistent <- inconsistent + 1
          }
        }
      }
    }
  }
  
  return(inconsistent)
}

inconsistencies <- check_z_consistency(data$y, z_inits)
cat("\nFound", inconsistencies, "inconsistencies in z initialization\n")


# Initial values
inits <- list(
  z = z_inits,
  phi.mean = array(0.5, dim = c(sp, station)),
  p.mean = array(0.5, dim = c(sp, station))
)

# Build and compile the model
cjs_model <- nimbleModel(
  cjs_code, 
  constants = constants, 
  data = data, 
  inits = inits
)

# Check model
cjs_model$calculate()  # Should not return -Inf if model is properly specified

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

# Run MCMC sampling
samples <- runMCMC(
  compiled_mcmc, 
  niter = 10000, 
  nburnin = 2000, 
  nchains = 2,
  samplesAsCodaMCMC = TRUE
)

# Examine results for convergence
gelman.diag(samples)

# Summarize posterior distributions
mcmc_summary <- MCMCsummary(samples)
print(mcmc_summary)

# Plot results for first species, first station
# Compare estimates to true values
for (s in 1:sp) {
  for (l in 1:station) {
    # Extract phi estimates for this species and station
    phi_est <- mcmc_summary[grep(paste0("phi\\[.*,", s, ",", l, "\\]"), rownames(mcmc_summary)),]
    # Extract p estimates
    p_est <- mcmc_summary[grep(paste0("p\\[.*,", s, ",", l, "\\]"), rownames(mcmc_summary)),]
    
    # Print comparison
    cat("\nSpecies", s, "Station", l, "\n")
    cat("True phi:", phi_true, "\n")
    cat("Est. phi (mean):", phi_est$mean, "\n")
    cat("True p:", p_true, "\n")
    cat("Est. p (mean):", p_est$mean, "\n\n")
  }
}

# Summarize results
summary(samples)
MCMCvis::MCMCsummary(samples)

